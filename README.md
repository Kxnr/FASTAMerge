# FASTA Merge

## Description

A tool to merge multiple FASTA alignments to the same reference sequence. Given multiple files with sequences
aligned to the same reference sequence, merge reference sequences and aligned sequences such that both alignment
and spacing are preserved. This tool may not be suitable for large files, as it reads the full contents of all
provided files into memory.

Note: This tool only does rudimentary alignment on the reference sequences. For best results, one of the reference
      sequences should be a superset of all other reference sequences. If this is not the case, this tool will fall
      back to aligning all reference sequences pairwise, which requires that all pairs of reference sequences overlap.

## Setup

### Dependencies

* python >= 3.6
* argopt
* pykxnr

### Installation

This project is available through pip as fasta-merge. To install, run:

`pip install fasta-merge`
 OR
`python -m pip install fasta-merge`

### Executing program

For basic usage, pass the tag name of the reference sequence and a list of files as input.
Output can be redirected to a file to save merged alignments.

```
fasta-merge <label> [<file>...] > output.fasta
```

## Help

```
fasta-merge -h
```

## Authors

* [Connor Keane](kxnr.me)

<!--- ## Version History --->
## License

This project is licensed under the GPL3 License - see the LICENSE.md file for details

## Acknowledgments

* Thanks to the CU Boulder Whiteley Lab

