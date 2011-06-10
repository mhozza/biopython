# Copyright 2011 by Andrew Sczesnak.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.AlignIO import MafIO

class _RFC_CSF_Shared():
    def _multiseq_to_seq_by_id(self, multiseq):
        if not isinstance(multiseq, MultipleSeqAlignment):
            raise ValueError("Expected a MultipleSeqAlignment object")

        seq_by_id = {}

        for seqrec in multiseq:
            seq_by_id[seqrec.id] = str(seqrec.seq)

        return seq_by_id

class _MultizMafAccess():
    def _get_mafindexer(self, seqname):
        """Returns a MafIndex object for seqname, loading one when necessary.

        Returns an existing one when possible.
        """

        if not seqname in self._mafindexers:
            maf_file = "%s/%s.maf" % (self._maf_dir, seqname)
            maf_index_file = "%s/%s.mafindex" % (self._maf_dir, seqname)

            self._mafindexers[seqname] = MafIO.MafIndex(maf_index_file,
                                                        maf_file,
                                                        self._target_seqname + "." + seqname)

        return self._mafindexers[seqname]

    def _fetch_splice(self, seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds):
        # verify the provided exon coordinates
        if len(exonStarts) <> len(exonEnds):
            raise ValueError("Every position in starts must have a match in ends")

        for exonstart, exonend in zip(exonStarts, exonEnds):
            if exonstart >= exonend:
                raise ValueError("Exon coordinates invalid (%s >= %s)" % (exonstart, exonend))

        # apply this logic only to transcripts containing a coding sequence
        if cdsEnd - cdsStart > 0:
            # accept only the exons containing the CDS, and trim the ends
            kept_starts = []
            kept_ends = []

            in_cds = False

            for start, end in zip(exonStarts, exonEnds):
                # CDS is fully contained within the exon
                if cdsStart >= start and cdsEnd <= end:
                    kept_starts.append(cdsStart)
                    kept_ends.append(cdsEnd)
                # exon contains the cdsStart: trim left end
                elif cdsStart >= start and cdsStart <= end:
                    kept_starts.append(cdsStart)
                    kept_ends.append(end)

                    in_cds = True
                # exon contains the cdsEnd: trim right end
                elif cdsEnd >= start and cdsEnd <= end:
                    kept_starts.append(start)
                    kept_ends.append(cdsEnd)

                    in_cds = False
                # exon is fully contained within the CDS
                elif start > cdsStart and end < cdsEnd:
                    kept_starts.append(start)
                    kept_ends.append(end)
                # exon contains no coding sequence
                else:
                    if in_cds == True:
                        raise ValueError("A non-coding exon may not occur within a CDS")
        elif cdsEnd - cdsStart == 0:
            kept_starts = exonStarts
            kept_ends = exonEnds
        elif cdsEnd - cdsStart < 0:
            raise ValueError("CDS length may not be less than zero!")

        # get a MafIndex object
        idx = self._get_mafindexer(seqname)

        # we don't specify the strand here because we have to assume that the
        # reference strand is positive, and it is if we're using multiz and have
        # specified the correct reference sequence when initializing this class
        multiseq = idx.get_spliced(kept_starts, kept_ends)

        # convert this MultipleSeqAlignment object into a dict
        seq_by_id = {}

        for seqrec in multiseq:
            # we're going to have to assume that the ID is in the format of
            # species.chromosome, otherwise we won't get very far
            species = seqrec.id.split(".", 1)[0]

            # we take into account strandedness at this point
            sequence = str(seqrec.seq.reverse_complement()) if strand == "-" else str(seqrec.seq)

            seq_by_id[species] = sequence

        return seq_by_id