# Copyright 2011 by Andrew Sczesnak.  All rights reserved.
# Revisions Copyright 2011 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.AlignIO support for the "maf" multiple alignment format.

You are expected to use this module via the Bio.AlignIO functions(or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

"""
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Generic import Alignment
from Bio.Align import MultipleSeqAlignment
from Interfaces import SequentialAlignmentWriter
import shlex

class MafWriter(SequentialAlignmentWriter):
    """Accepts a MultipleSeqAlignment object, writes a MAF file.
    
    More later.
    """
    
    def write_header(self):
        self.handle.write ("##maf version=1 scoring=none\n")
        self.handle.write ("# generated by Biopython\n\n")

    def _write_record (self, record):
        fields = ["s",
                  #In the MAF file format, spaces are not allowed in the id
                  "%-40s" % record.id.replace(" ","_"),
                  "%15s" % record.annotations.get("start", 0),
                  "%5s" % record.annotations.get("size", len(str(record.seq).replace("-",""))),
                  record.annotations.get("strand", "+"),
                  "%15s" % record.annotations.get("srcSize", 0),
                  str(record.seq)]
        self.handle.write(" ".join(fields) + "\n")

    def write_alignment(self, alignment):
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an alignment object")
        
        if len(set([len(x) for x in alignment])) > 1:
            raise ValueError("Sequences must all be the same length")

        all_ids = [x.id for x in alignment]

        if len(all_ids) != len(set(all_ids)):
            raise ValueError("Identifiers in each MultipleSeqAlignment must be unique")

        # for now, use ._annotations private property, but restrict keys to those
        # specifically supported by the MAF format, according to spec
        try:
            anno = " ".join(["%s=%s" % (x, y) for x, y in alignment._annotations.items() if x in ("score", "pass")])
        except AttributeError:
            anno = "score=0.00"

        self.handle.write("a %s\n" % (anno,))

        recs_out = 0

        for record in alignment:
            self._write_record(record)

            recs_out += 1

        self.handle.write("\n")

        return recs_out
        
def MafIterator(handle, seq_count = None, alphabet = single_letter_alphabet):
    in_a_bundle = False
    
    annotations = []
    records = []
    
    while True:
        # allows parsing of the last bundle without duplicating code
        try:
            line = handle.next()
        except StopIteration:
            line = ""
        
        if in_a_bundle:
            if line.startswith("s"):
                # add a SeqRecord to the bundle
                line_split = line.strip().split()

                if len(line_split) <> 7:
                    raise ValueError("Error parsing alignment - 's' line must have 7 fields")

                # s (literal), src (ID), start, size, strand, srcSize, text (sequence)
                anno = {"start": int(line_split[2]),
                        "size": int(line_split[3]),
                        "strand": line_split[4],
                        "srcSize": int(line_split[5])}
                        
                sequence = line_split[6]
                
                # interpret a dot/period to mean same the first sequence
                if "." in sequence:
                    if not records:
                        raise ValueError("Found dot/period in first sequence of alignment")
                        
                    ref = str(records[0].seq)
                    new = []
                    
                    for (s, r) in zip(sequence, ref):
                        new.append(r if s == "." else s)
                             
                    sequence = "".join(new)
                    
                records.append(SeqRecord(Seq(sequence, alphabet),
                               id = line_split[1],
                               name = line_split[1],
                               description = "",
                               annotations = anno))
            elif line.startswith("e") or \
                 line.startswith("i") or \
                 line.startswith("q"):
                # not implemented
                pass
            elif not line.strip():
                # end a bundle of records
                if seq_count is not None:
                    assert len(records) == seq_count
                    
                alignment = MultipleSeqAlignment(records, alphabet)
                #TODO - Introduce an annotated alignment class?
                #See also Bio/AlignIO/FastaIO.py for same requirement.        
                #For now, store the annotation a new private property:
                alignment._annotations = annotations
                
                yield alignment
                
                in_a_bundle = False
                
                annotations = []
                records = []
            else:
                raise ValueError("Error parsing alignment - unexpected line:\n%s" % (line,))
        elif line.startswith("a"):
            # start a bundle of records
            in_a_bundle = True
            
            annotations = dict([x.split("=") for x in line.strip().split()[1:]])
                
            if len([x for x in annotations.keys() if x not in ("score", "pass")]) > 0:
                raise ValueError("Error parsing alignment - invalid key in 'a' line")
        elif line.startswith("#"):
            # ignore comments
            pass
        elif not line:
            break

class MafIndex():
    def __init__(self, sqlite_file, maf_file, target_seqname):
        """Indexes or loads the index of a MAF file."""
        
        self.target_seqname = target_seqname
        
        # imports
        import os, itertools
        
        try:
            from sqlite3 import dbapi2 as _sqlite
            from sqlite3 import IntegrityError as _IntegrityError
            from sqlite3 import OperationalError as _OperationalError
        except ImportError:
            from Bio import MissingPythonDependencyError
            raise MissingPythonDependencyError("Requires sqlite3, which is "
                                               "included Python 2.5+")
        
        # make sure maf_file exists, then open it up
        if os.path.isfile(maf_file):
            self._maf_fp = open(maf_file, "r")
        else:
            raise ValueError("Error opening %s -- file not found" % (maf_file,))
        
        # if sqlite_file exists, use the existing db, otherwise index the file
        if os.path.isfile(sqlite_file):            
            con = _sqlite.connect(sqlite_file)
            self._con = con
            
            try:
                db_target = con.execute("SELECT value FROM meta_data WHERE key = 'target_seqname'").fetchone()[0]      
                
                if db_target != target_seqname:
                    raise ValueError("Provided database indexed for %s, expected %s" % (db_target, target_seqname))

                record_count = int(con.execute("SELECT value FROM meta_data WHERE key = 'record_count'").fetchone()[0])          
                
                if record_count == -1:
                    raise ValueError("Unfinished/partial database provided")
                    
                records_found = int(con.execute("SELECT COUNT(*) FROM offset_data").fetchone()[0])
                
                if records_found <> record_count:
                    raise ValueError("Expected %s records, found %s.  Corrupt index?" % (record_count, records_found))
            except _OperationalError, err:
                raise ValueError("Not a Biopython index database? %s" % err)
        else:
            con = _sqlite.connect(sqlite_file)
            self._con = con

            # generator function, returns index information for each bundle
            def _maf_index():
                line = self._maf_fp.readline()

                while line:
                    if line.startswith("a"):
                        # note the offset
                        offset = self._maf_fp.tell() - len(line)
                        
                        # search the following lines for a match to target_seqname
                        while True:
                            line = self._maf_fp.readline()
                            
                            if not line.strip() or line.startswith("a"):
                                raise ValueError("Target for indexing (%s) not found in this bundle" % (target_seqname,))
                            elif line.startswith("s"):
                                # s (literal), src (ID), start, size, strand, srcSize, text (sequence)
                                line_split = line.strip().split()
                                
                                if line_split[1] == target_seqname:
                                    start = int(line_split[2])
                                    end = int(line_split[2]) + int(line_split[3])
    
                                    yield (self._ucscbin(start, end), start, end, offset)
                                    
                                    break                    

                    line = self._maf_fp.readline()

            # make the tables
            con.execute("CREATE TABLE meta_data (key TEXT, value TEXT);")
            con.execute("INSERT INTO meta_data (key, value) VALUES ('record_count', -1);")
            con.execute("INSERT INTO meta_data (key, value) VALUES ('target_seqname', '%s');" % (target_seqname,))
            con.execute("CREATE TABLE offset_data (bin INTEGER, start INTEGER, end INTEGER, offset INTEGER);")
            
            insert_count = 0
            
            mafindex_func = _maf_index()
                        
            while True:
                batch = list(itertools.islice(mafindex_func, 100))
                if not batch: break
            
                con.executemany("INSERT INTO offset_data (bin, start, end, offset) VALUES (?,?,?,?);", batch)
                con.commit()
                
                insert_count += len(batch)
                
            con.execute("CREATE INDEX IF NOT EXISTS bin_index ON offset_data(bin);")
            con.execute("CREATE INDEX IF NOT EXISTS start_index ON offset_data(start);")
            con.execute("CREATE INDEX IF NOT EXISTS end_index ON offset_data(end);")
            
            con.execute("UPDATE meta_data SET value = '%s' WHERE key = 'record_count'" % (insert_count,))
            
            con.commit()
            
        # lastly, setup a MafIterator pointing at the open maf_file
        self._mafiter = MafIterator(self._maf_fp)

    @staticmethod
    def _region2bin(start, end):
        """Finds bins that a region may belong to.
        
        Converts a region to a list of bins that it may belong to, including largest
        and smallest bins.
        """

        bins = [0, 1]

        bins.extend(range (1 + (start >> 26), 2 + ((end - 1) >> 26)))
        bins.extend(range (9 + (start >> 23), 10 + ((end - 1) >> 23)))
        bins.extend(range (73 + (start >> 20), 74 + ((end - 1) >> 20)))
        bins.extend(range (585 + (start >> 17), 586 + ((end - 1) >> 17)))

        return set(bins)

    @staticmethod
    def _ucscbin(start, end):
        """Returns the smallest bin a given region will fit into.
        
        Taken from http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
        """

        binOffsets = [512+64+8+1, 64+8+1, 8+1, 1, 0]

        _binFirstShift = 17
        _binNextShift = 3

        startBin = start
        endBin = end - 1

        startBin >>= _binFirstShift
        endBin >>= _binFirstShift

        for i in xrange(0, len(binOffsets)):
            if startBin == endBin:
                return binOffsets[i] + startBin

            startBin >>= _binNextShift
            endBin >>= _binNextShift

        return 0
            
    def get_record(self, offset):
        """Retrieves a single MAF record located at the offset provided."""
        
        self._maf_fp.seek(offset)
        return self._mafiter.next()
        
    def search(self, starts, ends):
        """Searches index database for MAF records overlapping ranges provided."""
        
        # verify the provided exon coordinates
        if len(starts) <> len(ends):
            raise ValueError("Every position in starts must have a match in ends")
        
        for exonstart, exonend in zip(starts, ends):
            if exonstart >= exonend:
                raise ValueError("Exon coordinates invalid (%s >= %s)" % (exonstart, exonend))
        
        con = self._con
        
        # search for every exon, put hits in a dict where offset is the key
        # this will automatically prevent retrieval of duplicate records
        hits = {}
        
        for exonstart, exonend in zip(starts, ends):
            possible_bins = ", ".join(map(str, self._region2bin (exonstart, exonend)))
            
            result = con.execute("SELECT * FROM offset_data WHERE bin IN (%s) AND "
            "(end BETWEEN %s AND %s OR %s BETWEEN start AND end) ORDER BY start ASC;" \
            % (possible_bins, exonstart, exonend, exonend))
    
            rows = result.fetchall()
            
            for i in rows:
                hits[i[3]] = (i[1], i[2])
        
        # iterate through hits, fetching alignments from the MAF file and checking
        # to be sure we've retrieved the expected record
        for offset, (rec_start, rec_end) in hits.items():
            fetched = self.get_record (int(offset))
            
            for record in fetched:
                if record.id == self.target_seqname:
                    start = record.annotations["start"]
                    end = record.annotations["start"] + record.annotations["size"]
                    
                    if not (start == rec_start and end == rec_end):
                        raise ValueError("Expected %s-%s @ offset %s, found %s-%s" % \
                        (rec_start, rec_end, offset, start, end))
                        
            yield fetched
            
    def get_spliced(self, starts, ends, strand = "+"):
        """Returns a multiple alignment of the exact sequence range provided.
        
        Accepts two lists of start and end positions on target_seqname, representing
        exons to be spliced in silico.  Returns a MultipleSeqAlignment of the
        desired sequenced spliced together.
        """

        # validate strand
        if strand not in ("+", "-"):
            raise ValueError("Strand must be + or -")
                
        # pull all alignments that span the desired intervals
        fetched = [x for x in self.search(starts, ends)]
        
        # keep track of the expected letter count
        expected_letters = sum([y - x for x, y in zip(starts,ends)])
        
        # if there's no alignment, return filler for the assembly of the length given
        if len(fetched) == 0:
            return MultipleSeqAlignment([SeqRecord(Seq("N" * expected_letters),
                                                   id = self.target_seqname)])
        
        # find the intersection of all IDs in these alignments
        all_seqnames = set([x.id for y in fetched for x in y])
                
        # split every record by base position
        split_by_position = dict([(x, {}) for x in all_seqnames])
        
        # keep track of what the total number of (unspliced) letters should be
        total_rec_length = 0
        
        # track first strand encountered on the target seqname
        ref_first_strand = None
        
        for multiseq in fetched:
            # find the target_seqname in this MultipleSeqAlignment and use it to
            # set the parameters for the rest of this iteration
            for seqrec in multiseq:
                if seqrec.id == self.target_seqname:
                    try:
                        if ref_first_strand == None:
                            ref_first_strand = seqrec.annotations["strand"]
                            
                            if ref_first_strand not in ("+", "-"):
                                raise ValueError("Strand must be + or -")
                        elif ref_first_strand != seqrec.annotations["strand"]:
                            raise ValueError("Encountered strand='%s' on target seqname, expected '%s'" % \
                            (seqrec.annotations["strand"], ref_first_strand))
                    except KeyError:
                        raise ValueError("No strand information for target seqname (%s)" % \
                        (self.target_seqname,))
                        
                    rec_length = len(seqrec)
                    rec_start = seqrec.annotations["start"]
                    rec_end = seqrec.annotations["start"] + seqrec.annotations["size"]

                    total_rec_length += rec_end - rec_start
                    
                    # blank out these positions for every seqname
                    for seqrec in multiseq:
                        for pos in xrange(rec_start, rec_end):
                            split_by_position[seqrec.id][pos] = ""
                    
                    break
            else:
                raise ValueError("Did not find %s in alignment bundle" % (self.target_seqname,))
                
            # the true, chromosome/contig/etc position in the target seqname
            real_pos = rec_start
            
            for gapped_pos in xrange(0, rec_length):
                for seqrec in multiseq:
                    # keep track of this position's value for the target seqname
                    if seqrec.id == self.target_seqname: track_val = seqrec.seq[gapped_pos]
                    
                    split_by_position[seqrec.id][real_pos] += seqrec.seq[gapped_pos]
                    
                # increment the real_pos counter only when non-gaps are found in
                # the target_seqname, and we haven't reached the end of the record
                if track_val != "-" and real_pos < rec_end - 1: real_pos += 1
                
        # make sure the number of bp entries equals the sum of the record lengths
        if len(split_by_position[self.target_seqname]) <> total_rec_length:
            raise ValueError("Target seqname (%s) has %s records, expected %s" % \
            (self.target_seqname, len(split_by_position[self.target_seqname]), total_rec_length))

        # translates a position in the target_seqname sequence to its gapped length        
        realpos_to_len = dict([(x, len(y)) for x, y in split_by_position[self.target_seqname].items() if len(y) > 1])

        # splice together the exons            
        subseq = {}
        
        for seqid in all_seqnames:
            seq_split = split_by_position[seqid]
            seq_splice = []
            
            filler_char = "N" if seqid == self.target_seqname else "-"

            # iterate from start to end, taking bases from split_by_position when
            # they exist, using N or - for gaps when there is no alignment.            
            append = seq_splice.append
            
            for exonstart, exonend in zip(starts, ends):
                for real_pos in xrange(exonstart, exonend):
                    # if this seqname has this position, add it
                    if real_pos in seq_split:
                        append(seq_split[real_pos])
                    # if not, but it's in the target_seqname, add length-matched filler
                    elif real_pos in realpos_to_len:
                        append(filler_char * realpos_to_len[real_pos])
                    # it's not in either, so add a single filler character
                    else:
                        append(filler_char)
                        
            subseq[seqid] = "".join(seq_splice)

        # make sure we're returning the right number of letters
        if len(subseq[self.target_seqname].replace("-", "")) != expected_letters:
            raise ValueError("Returning %s letters for target seqname (%s), expected %s" % \
            (len(subseq[self.target_seqname].replace("-", "")), self.target_seqname, expected_letters))
                
        # check to make sure all sequences are the same length as the target seqname
        ref_subseq_len = len(subseq[self.target_seqname])
        
        for seqid, seq in subseq.items():
            if len(seq) <> ref_subseq_len:
                raise ValueError("Returning length %s for %s, expected %s" % \
                (len(seq), seqid, ref_subseq_len))
                
        # finally, build a MultipleSeqAlignment object for our final sequences
        result_multiseq = []
        
        for seqid, seq in subseq.items():
            seq = Seq(seq)
            
            seq = seq if strand == ref_first_strand else seq.reverse_complement()
            
            result_multiseq.append(SeqRecord(seq,
                                             id = seqid,
                                             name = seqid,
                                             description = ""))
                                             
        return MultipleSeqAlignment(result_multiseq)
        
