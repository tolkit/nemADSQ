#!/usr/bin/env python3
"""
Identify long reads with telomeric repeats at their start or end.
This script relies on the strandedness of the input motif and long reads by assuming that
the forward motif will only occur on the right end of the reads, 
while the reverse complement of the motif will occur on the left side.

Usage:
        filter_telomeric_reads.py [--reads FILE] [--string STR] 
                                  [--out FILE] [--times INT]
                                  [--buffered INT] [--min_len INT] 

options:
    -r FILE, --reads FILE   filename for uncompressed reads to check.
    -s STR, --string STR    telomeric repeat
                            [Default: TTAGGC]
    --times INT             minimum number of contiguous occurrences.
                            [Default: 3]
    -b INT, --buffered INT  Length of sequence to be kept in memory
                            [Default: 10000000]
    -o FILE, --out FILE     filename for telomeric reads.
                            [Default: telomericReads.fasta]
"""

import sys
import os
from docopt import docopt

__author__ = "Pablo Manuel Gonzalez de la Rosa"
__version__ = '1.2.0'

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def reverse_complement_sequence(seq):
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

search_space_inside_read = 15
def trim_find(seq, min_occur, motif_size, motif, rev_motif, longer_motif, longer_rev_motif):
    matches_start = seq.find(longer_rev_motif, 0, motif_size * search_space_inside_read)
    matches_end = seq.find(longer_motif, -motif_size * search_space_inside_read)
    is_telomeric = 0
    if matches_start >= 0:
        is_telomeric = 1
    elif matches_end >= 0:
        is_telomeric = 1
    return is_telomeric



if __name__ == "__main__":
    args                = docopt(__doc__)
    motif               = args['--string']
    fastafile           = args['--reads']
    min_occur           = int(args['--times'])
    buffer_size         = int(args['--buffered'])
    outfile             = args['--out']
    rev_motif           = reverse_complement_sequence(motif)
    longer_motif        = motif * min_occur
    longer_rev_motif    = rev_motif * min_occur
    motif_size          = len(motif)
    non_telomeric_reads     = ''
    seq_in_mem_len      = 0

    # parse to read stdin or read to file
    if fastafile == "-":
        infile = sys.stdin
        print("Reading from STDIN")
    else:
        infile = open(fastafile, 'r')



    if os.path.exists(outfile):
        print("[WARNING] Overwriting existing file: {}".format(outfile))
        os.remove(outfile)

    for name, seq, qual in readfq(infile):
        sequence_is_telomeric = trim_find(seq, min_occur, motif_size,
        motif, rev_motif, longer_motif, longer_rev_motif)
        if sequence_is_telomeric == 1:
            non_telomeric_reads += ">%s\n" % name
            non_telomeric_reads += "%s\n" % seq
            seq_in_mem_len += len(seq)
            if seq_in_mem_len > buffer_size:
                with open(outfile, 'a') as ofh:
                    ofh.writelines(non_telomeric_reads)
                non_telomeric_reads = ''
                seq_in_mem_len = 0

    # Write sequence left in buffer after loop ends
    if seq_in_mem_len > 0:
        with open(outfile, 'a') as ofh:
                    ofh.writelines(non_telomeric_reads)

    # If no read had telomeres, create an empty file
    if not os.path.exists(outfile):
        open(outfile, 'a').close()
        print("[WARNING] No reads were identified without telomeric repeat")
    

    # parse to close opened file
    if fastafile != "-":
        infile.close()
