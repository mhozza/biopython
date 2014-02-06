# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parsing AlignACE output files
"""

from Bio.Motif import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


class Record(object):
    def __init__(self):
        self.motifs=[]
        self.parameters = None


def read(handle):
    """read(handle)"""
    record = Record()
    line = handle.next()
    record.version = line.strip()
    line = handle.next()
    record.command = line.strip()
    for line in handle:
        line = line.strip()
        if line=="":
            pass
        elif line[:4]=="Para":
            record.parameters={}
        elif line[0]=="#":
            seq_name = line.split("\t")[1]
            record.sequences.append(seq_name)
        elif "=" in line:
            par_name, par_value = line.split("=")
            par_name = par_name.strip()
            par_value = par_value.strip()
            record.parameters[par_name]=par_value
        elif line[:5]=="Input":
            record.sequences=[]
        elif line[:5]=="Motif":
            current_motif = Motif()
            current_motif.alphabet=IUPAC.unambiguous_dna
            record.motifs.append(current_motif)
        elif line[:3]=="MAP":
            current_motif.score = float(line.split()[-1])
        elif len(line.split("\t"))==4:
            seq = Seq(line.split("\t")[0],IUPAC.unambiguous_dna)
            current_motif.add_instance(seq)
        elif "*" in line:
            current_motif.set_mask(line.strip("\n\c"))
        else:
            raise ValueError(line)
    return record
