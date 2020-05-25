#!/usr/bin/env python
"""
Consensus Creator ver 2.01

This progam reads pre-aligned sequences as a single FASTA file, and creates
consensus sequences. The program code was written by Yongchan Lee, 2017,
and modified by Mizuki Takemoto.
Theoretical framework is based on a paper by Steipe B et al., 1994.

-Reference
1) Steipe, B., Schiller, B., Pluckthun, A. & Steinbacher, S.
Sequence statistics reliably predict stabilizing mutations in a protein
domain. J. Mol. Biol. 240, 188-92 (1994).
"""

from __future__ import division, print_function

import argparse
import datetime
import re
from copy import deepcopy

import numpy as np


class AminoAcidsCounter(object):
    # TODO: key error if unusual letter is contained
    all_aa_letter = [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
        "-",
        "X",
        "J",
    ]

    def __init__(self, num):
        self.id = num
        self.dict = {}
        for k in self.all_aa_letter:
            self.dict[k] = 0

    def calc_sum(self):
        n_total = sum(self.dict.values())
        self.n_gap = self.dict.pop("-")
        return n_total - self.n_gap

    def get_max_aa(self):
        return max([(v, k) for k, v in self.dict.items()])[1]


class MultiFasta(object):
    def __init__(self, filename):
        self.filename = filename
        self.fasta_ary = []
        self.find_all_fasta(filename)

    def __getitem__(self, i):
        return self.fasta_ary[i]

    def find_all_fasta(self, filename):
        with open(filename, "r") as f:
            line = f.readline().rstrip("\r\n")

            # the first line must be a header
            if not re.match("^>", line):
                raise

            header = re.sub("^>", "", line)

            line = f.readline().rstrip("\r\n")
            seq = ""

            while line:
                if re.match("^>", line):
                    self.fasta_ary.append(SingleFasta(header, seq))
                    header = re.sub("^>", "", line)
                    seq = ""
                else:
                    seq += line
                line = f.readline().rstrip("\r\n")

    def get_all(self, n):
        return [fasta.seq[n] for fasta in self.fasta_ary]


class SingleFasta(object):
    def __init__(self, header, seq):
        self.header = header
        self.seq = seq

    @property
    def length(self):
        return len(self.seq)


d = datetime.datetime.today()
str_run_datetime = d.strftime("%Y%m%d-%H-%M-%S")
print("Run on", str_run_datetime)


# arg parser #
str_description = """\
************************ Consensus Creator ver 2.01 ***********************
This progam reads pre-aligned sequences as a single FASTA file, and creates
consensus sequences. The program code was written by Yongchan Lee, 2017,
and modified by Mizuki Takemoto.
Theoretical framework is based on a paper by Steipe B et al., 1994.

-Reference
1) Steipe, B., Schiller, B., Pluckthun, A. & Steinbacher, S.
Sequence statistics reliably predict stabilizing mutations in a protein
domain. J. Mol. Biol. 240, 188-92 (1994).
***************************************************************************
"""

p = argparse.ArgumentParser(
    description=str_description, formatter_class=argparse.RawTextHelpFormatter
)
p.add_argument(
    "-f", "--filename", type=str, help="input filename",
)
default_threshold = [0, 0.5, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3]
p.add_argument(
    "-t",
    "--threshold",
    help="you can specify threshold values seperated by space",
    type=float,
    nargs="*",
    default=default_threshold,
)
args = p.parse_args()

print("Import file destination:")
print(args.filename)
print("Threshold values:")
print(", ".join([str(e) for e in args.threshold]))


# Load Multi-Fasta file #
multifasta = MultiFasta(args.filename)

n_length = multifasta[0].length
n_seq = len(multifasta.fasta_ary)

# raise Error if all sequences do not have same length
bool_ary = [fasta.length == n_length for fasta in multifasta.fasta_ary]
if not all(bool_ary):
    raise ValueError


# make counts for each amino acids
counts = {}
for i in range(n_length):
    counter = AminoAcidsCounter(i)
    for aa in multifasta.get_all(i):
        counter.dict[aa] += 1
    counts[i] = counter


# sum up over all sequences
all_count = {}
for i in range(n_length):
    all_count[i] = counts[i].calc_sum()


# header
print("AA_Pos    WT  Freq    Con Freq    FCon/FWT")

seq_id = 1
target_seq = []
consensus_seq = []
ratio_ary = []

for i in range(n_length):
    target_aa = multifasta[0].seq[i]

    # skip gap
    if target_aa == "-":
        continue

    target_freq = counts[i].dict[target_aa] / all_count[i] * 100
    max_aa = counts[i].get_max_aa()
    max_freq = counts[i].dict[max_aa] / all_count[i] * 100

    ratio = max_freq / target_freq
    ratio_ary.append(ratio)

    # save amino acid sequences
    target_seq.append(target_aa)
    consensus_seq.append(max_aa)

    output = (
        "{seq_id:<6d}{target_aa:>6s}  {target_freq:5.1f}    "
        "{max_aa:<3s} {max_freq:5.1f}    {ratio:4.2f}"
    ).format(**locals())
    # output = "%s, %5.2f, %s, %5.2f" % (target_aa, target_freq, max_aa, max_freq)
    print(output)

    seq_id += 1


# output sequences
title = multifasta[0].header
print()
print(">WT_" + title)
print("".join(target_seq))
print()
print(">Consensus_" + title)
print("".join(consensus_seq))
print()


# print suggested mutation and seq
target_seq = np.array(target_seq)
consensus_seq = np.array(consensus_seq)
ratio_ary = np.array(ratio_ary)

for t in args.threshold:
    bool_ary = np.log(ratio_ary) > t

    target = deepcopy(target_seq)

    # mutation suggestion
    base_seq = ["-" for i in range(seq_id - 1)]
    base_seq = np.array(base_seq)
    base_seq[bool_ary] = consensus_seq[bool_ary]

    n_mutation = sum(bool_ary)
    header = ">Threshold_value:{t} (suggested_mutation:{n_mutation})".format(**locals())
    print(header)
    print("".join(base_seq))
    print()

    # construct suggestion
    target[bool_ary] = consensus_seq[bool_ary]
    header = ">Construct_threshold_value:{}".format(t)
    print(header)
    print("".join(target))
    print()
