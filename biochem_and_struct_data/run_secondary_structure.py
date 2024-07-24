#!/usr/bin/env/python
"""Runs secondary structure prediction on FASTA sequences."""

import argparse
from Bio import SeqIO
import collections
import logging
import os
import re
from seqfold import dg, fold
import sys

__author__ = "Ian Hoskins"
__credits__ = ["Ian Hoskins"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Ian Hoskins"
__email__ = "ianjameshoskins@utexas.edu"
__status__ = "Development"

__logger = logging.getLogger(__name__)

# Consider putting the following in a logging config file
__logger.setLevel(logging.DEBUG)
__fhandler = logging.FileHandler("stderr.log")
__fhandler.setLevel(logging.DEBUG)
__formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
__fhandler.setFormatter(__formatter)
__logger.addHandler(__fhandler)

FILE_NEWLINE = "\n"
FILE_DELIM = "\t"
NA_VALUE = "NA"


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    # Add new arguments for command line passing of files, options, etc; see argparse docs
    parser.add_argument("-f", "--fasta", type=str, help='Sequences in FASTA format.')

    parser.add_argument("-t", "--temp", type=int, default=RnaFolding.DEFAULT_FOLD_TEMP,
                        help='Integer temperature for folding. Default %i' % RnaFolding.DEFAULT_FOLD_TEMP)

    parser.add_argument("-m", "--max_length", type=int, default=RnaFolding.DEFAULT_MAX_LEN,
                        help='Maximum length to consider for folding. Default %i' % RnaFolding.DEFAULT_MAX_LEN)

    parser.add_argument("-o", "--outdir", type=str, default=".", help='Output directory. Default current directory.')

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def dna_to_rna(seq):
    """Converts Ts to Us.

    :param str seq: input sequence
    :return str: RNA sequence
    """

    dna_to_rna_trans = str.maketrans("Tt", "Uu")
    res = seq.translate(dna_to_rna_trans)
    return res


class RnaFolding(object):
    """Class for folding RNAs"""

    DEFAULT_FOLD_TEMP = 37
    STRUCT_STATS = collections.namedtuple(
        "STRUCT_STATS", "n_hairpins, n_bifurc, n_bulges, start_stem, max_stem_len, max_loop_len")

    STACK_RE = re.compile("STACK")
    HAIRPIN_RE = re.compile("HAIRPIN")
    BIFURC_RE = re.compile("BIFURCATION")
    BULGE_RE = re.compile("BULGE")

    DEFAULT_DG = 0
    DEFAULT_MAX_LEN = 200
    DEFAULT_MIN_STEM_LEN = 3
    DEFAULT_MAX_LOOP_LEN = 10

    def __init__(self, temperature=DEFAULT_FOLD_TEMP):
        """Constructor for RnaFolding

        :param float temperature: folding temperature in degrees Celsius
        """

        self.temperature = temperature

    def get_min_dG(self, query_seq):
        """Determines the delta delta G minimum free energy change.

        :param str query_seq: DNA or RNA sequence
        :return float: delta G minimum free energy change
        """

        res = dg(dna_to_rna(query_seq), temp=self.temperature)
        return res

    def get_min_ddG(self, query_seq, control_seq):
        """Determines the delta delta G minimum free energy change for two sequences.

        :param str query_seq: alternate DNA or RNA sequence
        :param str control_seq: reference DNA or RNA sequence
        :return float: delta delta G
        """

        res = self.get_min_dG(query_seq) - self.get_min_dG(control_seq)
        return res

    def has_structure(self, query_seq, dg_cutoff=DEFAULT_DG, min_stem_len=DEFAULT_MIN_STEM_LEN,
                      max_loop_len=DEFAULT_MAX_LOOP_LEN):
        """Determines if a sequence has structure and if so returns stats.

        :param str query_seq: sequence to fold
        :param float dg_cutoff: deltaG cutoff to consider for structure formation. Default 0.
        :param int min_stem_len: minimum number of paired bases to constitute a stem/hairpin
        :param int max_loop_len: max number of bases in a loop to constitute a hairpin
        :return None | collections.namedtuple: STRUCT_STATS
        """

        rna_seq = dna_to_rna(query_seq)
        dg = self.get_min_dG(rna_seq)

        if dg > dg_cutoff:
            return None

        folds = fold(rna_seq, temp=self.temperature)
        if len(folds) == 0:
            return None

        stem_lens = []
        loop_lens = []
        obs_stem_len = 0
        start_stem = 0  # should almost always be overwritten- it would be unusual to start a structure without a match
        n_hairpins = 0
        n_bifurc = 0
        n_bulges = 0

        for i, op in enumerate(folds):

            # Enumerate stem bases until we see a hairpin; allow bulges and bifurcations
            stack_match = self.STACK_RE.search(op.desc)
            hairpin_match = self.HAIRPIN_RE.search(op.desc)
            bifurc_match = self.BIFURC_RE.search(op.desc)
            bulge_match = self.BULGE_RE.search(op.desc)

            if i == 0 and stack_match:
                start_stem = op.ij[0][0] + 1

            if stack_match:
                obs_stem_len += 1
                continue

            # If we observe any bulges or bifurcations reset the counter
            if bifurc_match:
                n_bifurc += 1
                continue

            if bulge_match:
                n_bulges += 1
                continue

            if hairpin_match and obs_stem_len >= min_stem_len:
                # For some reason Struct reports the i of the hairpin as 1 based preceding it, so subtract 1
                obs_loop_len = op.ij[0][1] - op.ij[0][0] - 1
                if obs_loop_len <= max_loop_len:
                    n_hairpins += 1
                    stem_lens.append(obs_stem_len)
                    loop_lens.append(obs_loop_len)
                    # Reset the stem counter if we've seen a loop
                    obs_stem_len = 0

            stem_lens.append(obs_stem_len)

        max_stem_len = 0
        if len(stem_lens) != 0:
            max_stem_len = max(stem_lens)

        max_loop_len = 0
        if len(loop_lens) != 0:
            max_loop_len = max(loop_lens)

        struct_stats = self.STRUCT_STATS(
            n_hairpins=n_hairpins, n_bifurc=n_bifurc, n_bulges=n_bulges, start_stem=start_stem,
            max_stem_len=max_stem_len, max_loop_len=max_loop_len)

        return struct_stats


def workflow(fasta, temp=RnaFolding.DEFAULT_FOLD_TEMP, max_len=RnaFolding.DEFAULT_MAX_LEN, outdir="."):
    """Runs the secondary structure workflow.

    :param str fasta: path of the FASTA.
    :param int temp: temperature for fold. Default 37 C.
    :param int max_len: max length to consider a sequence for folding. Default 200.
    :param str outdir: Output directory to write results to.
    """

    rna_fold = RnaFolding(temperature=temp)
    outfile = os.path.join(outdir, os.path.basename(fasta) + ".sec.struct.txt")
    output_fields = ("seqname", "min_dG", "n_hairpins", "n_bifurc", "n_bulges", "start_stem", "max_stem_len", "max_loop_len",)

    with open(fasta, "r") as in_fh, \
            open(outfile, "w") as out_fh:

        outline = FILE_DELIM.join(output_fields) + FILE_NEWLINE
        out_fh.write(outline)

        for record in SeqIO.parse(in_fh, "fasta"):

            seqname = str(record.id)
            seq = str(record.seq)

            # Skip structure analysis on very long sequences
            if len(seq) > max_len:
                outline = FILE_DELIM.join([seqname] + [NA_VALUE] * 7) + FILE_NEWLINE
            else:
                min_dG = rna_fold.get_min_dG(seq)
                structs = rna_fold.has_structure(seq)

                if structs is None:
                    outline = FILE_DELIM.join([seqname, str(min_dG)] + [NA_VALUE] * 6) + FILE_NEWLINE
                else:
                    outline = FILE_DELIM.join([seqname, str(min_dG)] + list(
                        map(str, list(structs._asdict().values())))) + FILE_NEWLINE

            out_fh.write(outline)


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])
    workflow(fasta=parsed_args["fasta"], temp=parsed_args["temp"],
             max_len=parsed_args["max_length"], outdir=parsed_args["outdir"])


if __name__ == "__main__":
    __logger.info("Started %s" % sys.argv[0])
    main()
    __logger.info("Completed %s" % sys.argv[0])
