"""
FASTA file merger. Take multiple alignments to the same reference sequence and produce a FASTA formatted string
containing all given sequences aligned to the same reference sequence.

Usage:
    fasta-merge <ref> (<file> ...)

Arguments:
    <ref>       name of tag in files for reference sequence. Reference sequences must be trivially
                alignable either by pairwise alignment or matching to longest sequence.

    <file>      the filenames of the sequences to align
"""
import copy
from itertools import chain, groupby, zip_longest
from functools import reduce, partial
from collections import namedtuple
from argopt import argopt
from pykxnr.string_tools import pad, split_on_indices, merge_strings
from pykxnr.utils import clamp
import re

# TODO: offset_marks

__version__ = "1.0"

mark = namedtuple('mark', ('start', 'len'))


class Sequence:
    '''
    Class to represent a fasta sequence with spacing marks added to control alignment of
    sequences. Enables offsetting of sequences for alignment between objects, access to
    raw sequence data, modification of marks in sequence, and rendering sequence.
    '''

    def __init__(self, label, data, offset=0, mark_char='-'):
        self.label = label
        self.data = data.replace(mark_char, '')
        self._offset = offset

        padding = re.search(rf'^{mark_char}*', data)  # match from beginning of string
        self.start_padding = padding.end() - padding.start()

        padding = re.search(rf'{mark_char}*$', data)  # match from beginning of string
        self.end_padding = padding.end() - padding.start()
        self._marks = matches_offset(data.strip(mark_char), rf'[{mark_char}]+', global_offset=offset)

    def marks(self, padding=False):
        if padding:
            return [mark(self.offset, self.start_padding + self.offset)] \
                   + self._marks \
                   + [mark(len(self) + self.offset, self.end_padding)]

        return self._marks

    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, offset):
        '''
        modify marks when offset is changed to match newly offset sequence so that
        accessing marks gives marks aligned to the same reference.

        :param offset: integer sequence offset
        :return: None
        '''
        self._marks = [mark(m.start - self._offset + offset, m.len) for m in self._marks]
        self._offset = offset

    def sequence(self, raw=True):
        if raw:
            return self.data
        else:
            return add_marks(self.data, self.marks(padding=True), offset=-self.offset)

    def add_marks(self, *marks: mark):
        '''
        add marks to existing marks sequence, combining marks that start at the same location additively

        :param marks: marks to add to existing marks
        :return: Sequence instance, modified in-place
        '''
        def additive_merge(marks):
            return clamp(sum(m.len for m in marks), 0, None)

        self._marks = merge_marks(self._marks, marks, mode=lambda k, m: mark(k, additive_merge(m)))
        return self

    def add_padding(self, start=0, end=0):
        self.start_padding = clamp(self.start_padding + start, 0, None)
        self.end_padding = clamp(self.end_padding + end, 0, None)
        return self

    def __len__(self):
        return len(self.data)

    def __str__(self):
        return self.sequence(raw=False)

    def __repr__(self):
        return f'<Sequence {self.label}: {str(self)}>'

    @classmethod
    def from_sequences(cls, tag, sequences):
        '''
        return a new sequence that represents the union of mulitple sequences, using very basic alignment to
        merge sequences. If longest sequence is a superset of other sequences, return longest sequence with
        merged marks. Otherwise, attempt to merge sequences pairwise and then find offsets before merging marks.

        :param tag: name to give to new sequence
        :param sequences: Sequence objects to merge
        :return: new Sequence object for merged sequence
        '''
        sequences = copy.deepcopy(sequences)
        merged_string = merge_strings([s.sequence(raw=True) for s in sequences])

        for s in sequences:
            s.offset = merged_string.find(s.sequence(raw=True))

        merged_marks = merge_marks(*(s.marks(padding=True) for s in sequences))
        print(merged_marks)

        return Sequence(tag, add_marks(merged_string, merged_marks)), sequences


def add_marks(sequence: Sequence, marks: list[mark], fill='-', offset=0):
    '''
    given a sequence, add 'fill' element to locations specified by gaps

    :param sequence: string to add fill element to
    :param marks: iterable of (start, length) tuples describing fill locations, indexed to unmodified sequence
    :param fill: item to repeat in fill
    :param offset: amount to move start positions when adding marks
    :return: new, modified string
    '''
    starts, lens = zip(*marks)
    split = split_on_indices(sequence, [clamp(s + offset, 0, len(sequence)) for s in starts])
    return ''.join([seq + fill * l for l, seq in zip_longest(lens, split, fillvalue=0)])


def matches_offset(seq: str, pattern=r'[-]+', global_offset=0):
    '''
    find all matches of compiled regex to string and adjust start
    indices to reflect position in string if all matches were removed

    :param seq: string to search
    :param pattern: compiled regex (re.compile)
    :param global_offset: an optional offest to add to the start of all matches
    :return: list of (start, length) tuples with adjusted start positions
    '''

    marks = [mark(m.start(), m.end() - m.start()) for m in re.finditer(pattern, seq)]

    local_offset = 0
    for i, m in enumerate(marks):
        marks[i] = mark(m.start + global_offset - local_offset, m.len)
        local_offset += m.len

    return marks


def merge_marks(*marks: list[mark], mode=partial(max, key=lambda x: x.len)):
    '''
    merge sequences of (start, length) marks for multiple sequences. Where
    marks share a start position, always take the longest.

    :param marks: lists of marks for different sequences, with each mark a (start, length) tuple
    :param mode: function that takes a collection of marks and returns the merged mark

    :return:  a single list of marks combining the marks from each sequence
    '''
    return [mode(g) for key, g in groupby(sorted(chain(*marks)), key=lambda x: x.start)]


def diff_marks(marks: list[mark], ref_marks: list[mark]):
    '''
    Given two sets of marks, return a new set of marks needed to turn the first into the second

    Note: it is an error to have multiple marks with the same start in the same sequence
    Note: this routine assumes that marks are strictly added between marks and ref_marks

    :param marks: list of marks to modify
    :param ref_marks: list of marks to diff against
    :return: new marks to add to marks to match ref marks
    '''
    # we strictly add marks, so we can just subtract within groups to find the changes
    merged = merge_marks(marks, ref_marks)
    return [reduce(lambda a, b: mark(key, b.len - a.len), g) for key, g in
            groupby(sorted(marks + merged), key=lambda x: x.start)]


def reindex_marks(to_adjust: list[mark], adjustment: list[mark]):
    '''
    given a collection of marks, add adjustment lengths to marks
    after adjustment positions, reindexing marks to a string with
    adjustment marks already present

    :param to_adjust: list of marks to adjust
    :param adjustment: list of marks already present
    :return: new list of marks with added offsets
    '''
    adjusted = to_adjust[:]

    i = 0
    start, adj = adjustment[i]
    total_adj = 0
    for j, val in enumerate(to_adjust):
        while val.start > start and i < len(adjustment) - 1:
            i += 1
            total_adj += adj
            start, adj = adjustment[i]

        adjusted[j] = (val.start + total_adj, val.len)

    return adjusted


def dict_to_fasta_str(data: dict):
    '''
    really simple fasta formatter, takes a dict of data and returns
    a fasta formatted string

    :param data: dict of {key: sequence} pairs
    :return: fasta formatted string
    '''
    max_len = len(max(data.values(), key=lambda x: len(x)))
    return "\n".join([f'>{k}\n{pad(v, 0, clamp(max_len - len(v), 0, max_len))}' for k, v in data.items()])


def faster_fasta_reader(fname: str):
    '''
    Really simple reader to process a fasta file into a dict of strings.
    Not suitable for very large files or cases where metadata is essential

    :param fname: filename of fasta file
    :return: dict of sequences
    '''

    # TODO: raise exception on invalid format

    genes = {}
    with open(fname, 'r') as f:
        for line in f.readlines():
            if line.startswith('>'):
                label = line.strip().split('|')[0][1:]
                genes[label] = ''
            else:
                genes[label] += line.strip()

    return genes

def main():
    parser = argopt(__doc__, version=__version__)
    args = parser.parse_args()

    tag = args.ref
    genes = [faster_fasta_reader(f) for f in args.file]

    seqs = [Sequence(tag, g[tag]) for g in genes]
    full_seq, seqs_offset = Sequence.from_sequences('reference', seqs)
    diffs = [reindex_marks(diff_marks(so.marks(padding=True), full_seq.marks(padding=False)), so.marks(padding=True))
             for so in seqs_offset]

    realigned = [Sequence(k, s, offset=so.offset).add_marks(d).add_padding(full_seq.start_padding * (so.offset > 0))
                 for so, d, g in zip(seqs_offset, diffs, genes) for k, s in g.items() if not k == tag]
    realigned.append(full_seq)

    print(dict_to_fasta_str({seq.label: seq.sequence(raw=False) for seq in realigned}))


if __name__ == '__main__':
    main()
