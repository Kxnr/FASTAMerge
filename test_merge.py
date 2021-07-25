from unittest import TestCase
from merge import *
from tempfile import NamedTemporaryFile
import re

SEQ_1 = '----AAAAAAA------AAAAAAAA----'
SEQ_2 = 'A-AA--AAA---AAAA----AAAAA'
MERGED = '----A-AA--AAA---A------AAA----AAAAA----'
FASTA = \
    f'''>1
----AA----
>2
---AAAAA-A'''
FASTA_DICT = {'1': '----AA----', '2': '---AAAAA-A'}


class TestStringMethods(TestCase):
    def test_split_on_indices(self):
        test_str = 'a' * 10
        self.assertEqual(split_on_indices(test_str, (5,)), ('a' * 5, 'a' * 5))
        self.assertEqual(split_on_indices(test_str, (2, 6)), ('a' * 2, 'a' * 4, 'a' * 4))
        self.assertEqual(split_on_indices(test_str, (2, 10)), ('a' * 2, 'a' * 8, ''))
        self.assertEqual(split_on_indices(test_str, (0, 10)), ('', test_str, ''))

    def test_pad(self):
        test_str = 'a' * 10
        offset = 3
        l = 20
        padded = pad(test_str, offset, l-offset-len(test_str), char='-')
        self.assertEqual(padded.find(test_str), offset)
        self.assertEqual(len(padded), l)

    def test_merge_strings(self):
        strings = ['abbccdd', 'bbcc', 'ccddeeff', 'aabbcc']
        self.assertEqual(merge_strings(strings), 'aabbccddeeff')


class TestMarkMethods(TestCase):
    def setUp(self):
        self.pattern = re.compile(r'(-+)')

    def test_merge_marks(self):
        marks_1 = [(0, 1), (2, 3)]
        marks_2 = [(0, 1), (2, 5), (7, 8)]
        self.assertEqual(merge_marks([marks_1, marks_1]), marks_1)
        self.assertEqual(merge_marks([marks_1, marks_2]), marks_2)

    def test_diff_marks(self):
        marks_1 = [(0, 1), (2, 3)]
        marks_2 = [(0, 1), (2, 4), (9, 3)]
        self.assertEqual(diff_marks(marks_1, marks_1), [(1, 0), (5, 0)])
        self.assertEqual(diff_marks(marks_1, marks_2), [(0, 0), (2, 1), (9, 3)])

    def test_reindex_marks(self):
        marks = [(0, 1), (2, 4), (9, 3)]
        adj = [(0, 2), (1, 2), (2, 6), (10, 11)]
        expected = [(0, 1), (6, 4), (19, 3)]
        self.assertEqual(reindex_marks(marks, adj), expected)

    def test_add_marks(self):
        marks = [(1, 1), (3, 2), (6, 3), (10, 4)]
        self.assertEqual(add_marks('A' * 15, marks), SEQ_2)

    def test_matches_offset(self):
        expected = [(1, 1), (3, 2), (6, 3), (10, 4)]
        self.assertEqual(matches_offset(SEQ_2), expected)


class TestFastaHandling(TestCase):
    def setUp(self):
        self.file = NamedTemporaryFile(mode='w+')
        self.file.write(FASTA)
        self.file.seek(0)

    def test_dict_to_fasta_str(self):
        self.assertEqual(dict_to_fasta_str(FASTA_DICT), FASTA)

    def test_faster_fasta_reader(self):
        self.assertEqual(faster_fasta_reader(self.file.name), FASTA_DICT)

    def tearDown(self) -> None:
        self.file.close()


class IntegrationTests(TestCase):
    def test_diff_add(self):
        diff = diff_marks(matches_offset(SEQ_1), matches_offset(MERGED))
        reindexed = reindex_marks(diff, matches_offset(SEQ_1))
        self.assertEqual(add_marks(SEQ_1, reindexed), MERGED)
