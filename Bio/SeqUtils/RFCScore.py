# Copyright 2011 by Andrew Sczesnak.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Implements the reading frame conservation (RFC) score test described by
Manolis Kellis in:
    
    http://web.mit.edu/manoli/www/thesis/Chapter2.html

"""

from __rfc_csf_shared import _MultizMafAccess, _RFC_CSF_Shared

class RFCScore(_RFC_CSF_Shared):
    def __init__(self, target_seqname):
        self._target_seqname = target_seqname
        
    def _score_alignment(self, seq_by_id):
        def _convert_to_123(sequence, frame):
            cur_num = frame
            new_seq = []
            append = new_seq.append
            
            for i in sequence:
                if i == "-":
                    append(i)
                else:
                    append(cur_num)
                    
                    if cur_num == 3:
                        cur_num = 1
                    else:
                        cur_num += 1
                        
            return "".join(map(str, new_seq))

        if self._target_seqname not in seq_by_id:
            raise ValueError("Specified target_seqname (%s) not in this alignment." % \
                             self._target_seqname)
   
        # convert target_seqname to 123s
        ref_seq_123 = _convert_to_123(seq_by_id[self._target_seqname], 1)

        # for every sequence
        rfc_by_id = {}
        
        for seqid, seq in seq_by_id.iteritems():
            # convert to 123s in all frames
            for_all_frames = [_convert_to_123(seq, x) for x in xrange(1, 4)]
            
            # for every frame, find % in common between reference/foreign
            perc_result = []
            
            for for_seq_123 in for_all_frames:
                gapless_pos = 0
                in_common = 0
                
                for x, y in zip(ref_seq_123, for_seq_123):
                    if x != "-" and y != "-":
                        gapless_pos += 1
                        
                        if x == y:
                            in_common += 1
                
                try:
                    perc_result.append(float(in_common) / float(gapless_pos))
                except ZeroDivisionError:
                    perc_result.append(0)
                
            rfc_by_id[seqid] = max(perc_result)
            
        return rfc_by_id
        
    def score_alignment(self, multiseq):
        return self._score_alignment(self._multiseq_to_seq_by_id(multiseq))
        
    def _sliding_window(self, seq_by_id, win_len, win_overlap):
        def _slice_seq_by_id(seq_by_id, from_pos, to_pos):
            new_seq_by_id = {}
            
            for seqid, seq in seq_by_id.iteritems():
                new_seq_by_id[seqid] = seq[from_pos:to_pos]
                
            return new_seq_by_id
                
        rfc_by_id = {}
        
        for start_pos in xrange(0, len(seq_by_id[self._target_seqname]), win_overlap):
            this_seq_by_id = _slice_seq_by_id(seq_by_id, start_pos, start_pos + win_len)
            
            this_score = self._score_alignment(this_seq_by_id)
            
            for seqid, score in this_score.iteritems():
                try:
                    rfc_by_id[seqid].append(score)
                except KeyError:
                    rfc_by_id[seqid] = [score]
                    
        # take the average
        for seqid, scorelist in rfc_by_id.iteritems():
            try:
                rfc_by_id[seqid] = float(sum(scorelist)) / float(len(scorelist))
            except ZeroDivisionError:
                rfc_by_id[seqid] = 0

        return rfc_by_id
                
    def sliding_window(self, multiseq, win_len, win_overlap):
        return self._sliding_window(self._multiseq_to_seq_by_id(multiseq), win_len, win_overlap)
        
class RFCMultizMaf(RFCScore, _MultizMafAccess):
    def __init__(self, maf_dir, target_seqname):
        self._maf_dir = maf_dir
        self._mafindexers = {}

        RFCScore.__init__(self, target_seqname)
        
    def score_alignment(self, seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds):
        return self._score_alignment(self._fetch_splice(seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds))
        
    def sliding_window(self, seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds, win_len, win_overlap):
        return self._sliding_window(self._fetch_splice(seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds), win_len, win_overlap)
                
if __name__ == "__main__":
    rfc = RFCMultizMaf("/comp_sync/data/foreign/ucsc/20100921_multiz30way", "mm9")

    import MySQLdb

    mysql_conn = MySQLdb.connect(host = "genome-mysql.cse.ucsc.edu",
                                 user = "genome",
                                 passwd = "",
                                 db = "mm9")

    db_conn = mysql_conn.cursor(MySQLdb.cursors.DictCursor)

    db_conn.execute ("SELECT * FROM refGene WHERE chrom = 'chr10' LIMIT 50")

    for record in db_conn.fetchall():
        print "scoring %s" % (record["name"],)

        # seqname, strand, cdsStart, cdsEnd, exonStarts, exonEnds
        print rfc.sliding_window(record["chrom"],
                                 record["strand"],
                                 int(record["cdsStart"]),
                                 int(record["cdsEnd"]),
                                 map(int, record["exonStarts"].split(",")[:-1]),
                                 map(int, record["exonEnds"].split(",")[:-1]),
                                 100,
                                 50)
        