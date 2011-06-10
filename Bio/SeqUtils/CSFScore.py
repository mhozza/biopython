# Copyright 2011 by Andrew Sczesnak.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Implements the codon substitution frequency (CSF) score test put forth by
Manolis Kellis and Michael Lin in:

    http://www.ncbi.nlm.nih.gov/pubmed/17989253

"""

from __rfc_csf_shared import _MultizMafAccess, _RFC_CSF_Shared

class CSFScore(_RFC_CSF_Shared):
    def __init__(self, target_seqname):
        # TODO
        # pull these codon tables from within biopython, maybe based on
        # a provided translation table code
        self._all_codons = ["aaa", "aat", "aag", "aac", "ata", "att", "atg",
                            "atc", "aga", "agt", "agg", "agc", "aca", "act",
                            "acg", "acc", "taa", "tat", "tag", "tac", "tta",
                            "ttt", "ttg", "ttc", "tga", "tgt", "tgg", "tgc",
                            "tca", "tct", "tcg", "tcc", "gaa", "gat", "gag",
                            "gac", "gta", "gtt", "gtg", "gtc", "gga", "ggt",
                            "ggg", "ggc", "gca", "gct", "gcg", "gcc", "caa",
                            "cat", "cag", "cac", "cta", "ctt", "ctg", "ctc",
                            "cga", "cgt", "cgg", "cgc", "cca", "cct", "ccg",
                            "ccc"]

        self._stop_codons = ["tag", "taa", "tga"]
        self._target_seqname = target_seqname

    def __str__(self):
        from Bio.Seq import Seq

        output = []

        matrix = self.normal_matrix

        for seqid in matrix.iterkeys():
            output.append("\n>> %s" % (seqid,))

            header1 = "      "
            header2 = "      "

            for foreign_codon in self._all_codons:
                header1 += "  %s " % (str(Seq(foreign_codon).translate()),)
                header2 += " %s" % (foreign_codon,)

            output.append(header1)
            output.append(header2)

            for ref_codon in self._all_codons:
                ref_aa = str(Seq(ref_codon).translate())
                this_codon_list = "%s %s " % (ref_aa, ref_codon)

                for foreign_codon in self._all_codons:
                    foreign_aa = str(Seq(foreign_codon).translate())

                    try: value = matrix[seqid][ref_codon][foreign_codon]
                    except KeyError: value = "X"

                    if value in (None, "X"):
                        value = "X"
                    else:
                        value = int(round(value, 0))

                    if ref_aa == foreign_aa:
                        this_addt = " >%2s" % (value,)
                    else:
                        this_addt = "%4s" % (value,)

                    this_codon_list += this_addt

                output.append(this_codon_list)

        return "\n".join(output)

    def _split_into_codons(self, seq_by_id):
        codon_start = None
        num_bases = 0

        split_alignment = dict([(x, []) for x in seq_by_id.iterkeys()])

        for i in xrange(0, len(seq_by_id[self._target_seqname])):
            if seq_by_id[self._target_seqname][i] == "-":
                continue
            else:
                num_bases += 1

                if num_bases == 1:
                    codon_start = i
                elif num_bases == 3:
                    for seqid in seq_by_id.iterkeys():
                        split_alignment[seqid].append (seq_by_id[seqid][codon_start:i + 1])

                    num_bases = 0
                    codon_start = None

        return split_alignment

    @staticmethod
    def _fast_median(numlist):
        numlist_len = len(numlist)

        if numlist_len == 0:
            return 0
        elif numlist_len == 1:
            return numlist[0]
        else:
            sorted_list = sorted(numlist)

            target_key = numlist_len / 2

            if numlist_len % 2 == 0:
                median = float(numlist[target_key] + numlist[target_key - 1]) / float(2)
            else:
                median = float(numlist[target_key])

            return median

    def _score_split_alignment(self, split_alignment):
        # empty dict with key as codon number in the split alignment
        codon_scores = dict([(x, []) for x in xrange(0, len(split_alignment[self._target_seqname]))])

        for seqid in split_alignment.iterkeys():
            # don't score against the target_seqname
            if seqid == self._target_seqname: continue

            for codon_num, ref_codon, foreign_codon in zip(codon_scores.keys(), split_alignment[self._target_seqname], split_alignment[seqid]):
                # in this case, it won't matter if we remove gaps
                if ref_codon.count("-") != foreign_codon.count("-"):
                    score = None
                else:
                    # remove gaps that occur in the same position in both codons
                    (ungapped_ref_codon, ungapped_foreign_codon) = map("".join, zip(*filter(lambda x:x not in [("-", "-")], zip(ref_codon, foreign_codon))))

                    ref_codon = ungapped_ref_codon.lower()
                    foreign_codon = ungapped_foreign_codon.lower()

                    if ("-" in ref_codon or "-" in foreign_codon) or \
                       ("n" in ref_codon or "n" in foreign_codon) or \
                       (ref_codon == foreign_codon) or \
                       (ref_codon in self._stop_codons):
                        score = None
                    else:
                        score = self.normal_matrix[seqid][ref_codon][foreign_codon]

                if score != None:
                    codon_scores[codon_num].append(score)

        # take median score at each codon position, add to final score
        fast_median = self._fast_median

        return sum(map(fast_median, codon_scores.itervalues()))

    def _train_matrix (self, split_alignment, matrix):
        for seqid in split_alignment.iterkeys():
            for gapped_ref_codon, gapped_foreign_codon in zip(split_alignment[self._target_seqname], split_alignment[seqid]):
                (ref_codon, foreign_codon) = map("".join, zip(*filter(lambda x:x not in [("-", "-")], zip(gapped_ref_codon, gapped_foreign_codon))))

                ref_codon = ref_codon.lower()
                foreign_codon = foreign_codon.lower()

                if ("-" in ref_codon or "-" in foreign_codon) or \
                   ("n" in ref_codon or "n" in foreign_codon) or \
                   (ref_codon == foreign_codon) or \
                   (ref_codon in self._stop_codons):
                    continue

                matrix[seqid][ref_codon][foreign_codon] += 1

        return True

    def _normalize_matrix(self, matrix):
        for asm in self._seqnames:
            for i in self._all_codons:
                ref_codon_total = 0

                for j in self._all_codons:
                    ref_codon_total += matrix[asm][i][j]

                for j in self._all_codons:
                    if ref_codon_total > 0:
                        matrix[asm][i][j] = float(matrix[asm][i][j]) / float(ref_codon_total)
                    else:
                        matrix[asm][i][j] = 0

    def blank_matrix(self, seqnames):
        self._seqnames = seqnames

        self.coding_matrix = {}
        self.noncoding_matrix = {}
        self.normal_matrix = {}

        for matrix in [self.coding_matrix, self.noncoding_matrix, self.normal_matrix]:
            for asm in seqnames:
                matrix[asm] = {}
                for i in self._all_codons:
                    matrix[asm][i] = {}
                    for j in self._all_codons:
                        matrix[asm][i][j] = 0

    def crunch_matrix(self):
        from math import log as math_log

        coding = self.coding_matrix
        noncoding = self.noncoding_matrix
        normal = self.normal_matrix

        self._normalize_matrix(coding)
        self._normalize_matrix(noncoding)

        for asm in self._seqnames:
            for i in self._all_codons:
                for j in self._all_codons:
                    try:
                        normal[asm][i][j] = math_log(float(coding[asm][i][j]) / float(noncoding[asm][i][j]))
                    except (ZeroDivisionError, ValueError):
                        normal[asm][i][j] = None

        normal["_meta_data"] = {"check": True}

        return True

    def load_matrix(self, pickle_file):
        import cPickle, os

        if os.path.isfile(pickle_file):
            with open(pickle_file, "r") as fp:
                pickle_obj = cPickle.load(fp)

            if pickle_obj["_meta_data"]["check"] != True:
                raise ValueError("Error loading matrix - incomplete/corrupt: %s" % (pickle_file,))

            self.normal_matrix = pickle_obj
        else:
            raise ValueError("Error loading matrix - file not found: %s" % (pickle_file,))

    def save_matrix(self, pickle_file):
        import cPickle, os

        with open(pickle_file, "w") as fp:
            cPickle.dump(self.normal_matrix, fp)

        return True

    def train_coding(self, multiseq):
        return self._train_matrix(self._split_into_codons(self._multiseq_to_seq_by_id(multiseq)), self.coding_matrix)

    def train_noncoding(self, multiseq):
        return self._train_matrix(self._split_into_codons(self._multiseq_to_seq_by_id(multiseq)), self.noncoding_matrix)

    def _sliding_window_score(self, split_alignment):
        """Does the actual work"""

        # not implemented
        raise Warning("Not implemented")

    def sliding_window_score(self, multiseq):
        """Lets another class implement functionality on top of us."""

        return self._sliding_window_score(self._split_into_codons(self._multiseq_to_seq_by_id(multiseq)))

    def absolute_score(self, multiseq):
        return self._score_split_alignment(self._split_into_codons(self._multiseq_to_seq_by_id(multiseq)))

class CSFMultizMaf(CSFScore, _MultizMafAccess):
    def __init__(self, maf_dir, target_seqname):
        self._maf_dir = maf_dir
        self._mafindexers = {}

        CSFScore.__init__(self, target_seqname)

    def sliding_window_score(self, seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds):
        return self._sliding_window_score(self._split_into_codons(self._fetch_splice(seqname,
                                                                   strand,
                                                                   cdsStart,
                                                                   cdsEnd,
                                                                   exonStarts,
                                                                   exonEnds)))

    def absolute_score(self, seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds):
        return self._score_split_alignment(self._split_into_codons(self._fetch_splice(seqname,
                                                                   strand,
                                                                   cdsStart,
                                                                   cdsEnd,
                                                                   exonStarts,
                                                                   exonEnds)))

    def train_coding(self, seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds):
        return self._train_matrix(self._split_into_codons(self._fetch_splice(seqname,
                                                                   strand,
                                                                   cdsStart,
                                                                   cdsEnd,
                                                                   exonStarts,
                                                                   exonEnds)),
                                  self.coding_matrix)

    def train_noncoding(self, seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds):
        return self._train_matrix(self._split_into_codons(self._fetch_splice(seqname,
                                                                   strand,
                                                                   cdsStart,
                                                                   cdsEnd,
                                                                   exonStarts,
                                                                   exonEnds)),
                                  self.noncoding_matrix)

if __name__ == "__main__":
    species = ["mm9", "rn4", "cavPor2", "oryCun1", "hg18", "panTro2", "ponAbe2",
               "rheMac2", "calJac1", "otoGar1", "tupBel1", "sorAra1", "eriEur1",
               "canFam2", "felCat3", "equCab1", "bosTau3", "dasNov1", "loxAfr1",
               "echTel1", "monDom4", "ornAna1", "galGal3", "anoCar1", "xenTro2",
               "tetNig1", "fr2", "gasAcu1", "oryLat1", "danRer5"]

    csf = CSFMultizMaf("/comp_sync/data/foreign/ucsc/20100921_multiz30way", "mm9")

    csf.blank_matrix(species)

    import sys, MySQLdb

    mysql_conn = MySQLdb.connect(host = "genome-mysql.cse.ucsc.edu",
                                 user = "genome",
                                 passwd = "",
                                 db = "mm9")

    db_conn = mysql_conn.cursor(MySQLdb.cursors.DictCursor)

    db_conn.execute ("SELECT * FROM refGene WHERE chrom = 'chr10' AND cdsEnd - cdsStart > 0 LIMIT 10")

    for record in db_conn.fetchall():
        print "training coding sequence %s" % (record["name"],)

        # seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds
        csf.train_coding(record["chrom"],
                         record["strand"],
                         int(record["cdsStart"]),
                         int(record["cdsEnd"]),
                         map(int, record["exonStarts"].split(",")[:-1]),
                         map(int, record["exonEnds"].split(",")[:-1]))

    db_conn.execute ("SELECT * FROM refGene WHERE chrom = 'chr10' AND cdsEnd - cdsStart = 0 LIMIT 10")

    for record in db_conn.fetchall():
        print "training non-coding sequence %s" % (record["name"],)

        # seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds
        csf.train_noncoding(record["chrom"],
                         record["strand"],
                         int(record["cdsStart"]),
                         int(record["cdsEnd"]),
                         map(int, record["exonStarts"].split(",")[:-1]),
                         map(int, record["exonEnds"].split(",")[:-1]))

    csf.crunch_matrix()
    
    csf.save_matrix("/tmp/matrix.pickle")

    db_conn.execute ("SELECT * FROM refGene WHERE chrom = 'chr10'")

    for record in db_conn.fetchall():
        print "scoring %s" % (record["name"],)

        # seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds
        print csf.absolute_score(record["chrom"],
                         record["strand"],
                         int(record["cdsStart"]),
                         int(record["cdsEnd"]),
                         map(int, record["exonStarts"].split(",")[:-1]),
                         map(int, record["exonEnds"].split(",")[:-1]))
