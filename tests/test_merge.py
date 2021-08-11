from unittest import TestCase
from fasta_merge.merge import *
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


class TestMarkMethods(TestCase):
    def setUp(self):
        self.pattern = re.compile(r'(-+)')

    def test_merge_marks(self):
        marks_1 = [mark(0, 1), mark(2, 3)]
        marks_2 = [mark(0, 1), mark(2, 5), mark(7, 8)]
        self.assertEqual(merge_marks(marks_1, marks_1), marks_1)
        self.assertEqual(merge_marks(marks_1, marks_2), marks_2)

    def test_diff_marks(self):
        marks_1 = [mark(0, 1), mark(2, 3)]
        marks_2 = [mark(0, 1), mark(2, 4), mark(9, 3)]
        self.assertEqual(diff_marks(marks_1, marks_1), [mark(0, 0), mark(2, 0)])
        self.assertEqual(diff_marks(marks_1, marks_2), [mark(0, 0), mark(2, 1), mark(9, 3)])

    def test_reindex_marks(self):
        marks = [mark(0, 1), mark(2, 4), mark(9, 3)]
        adj = [mark(0, 2), mark(1, 2), mark(2, 6), mark(10, 11)]
        expected = [mark(0, 1), mark(6, 4), mark(19, 3)]
        self.assertEqual(reindex_marks(marks, adj), expected)

    def test_add_marks(self):
        marks = [mark(1, 1), mark(3, 2), mark(6, 3), mark(10, 4)]
        self.assertEqual(add_marks('A' * 15, marks), SEQ_2)

    def test_matches_offset(self):
        expected = [mark(1, 1), mark(3, 2), mark(6, 3), mark(10, 4)]
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
