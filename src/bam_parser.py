import os
import re
import pysam
import numpy as np

from data_structures import Alignment


_CIGAR_RE = re.compile(r'\d+[A-Za-z=]')

def _parse_cigar(cigar_string: str, is_reverse: bool) -> tuple[int, int, int, int, list[tuple[int, int]], list[tuple[int, int]]]:
    op_pairs = _CIGAR_RE.findall(cigar_string)

    ref_pos = 0
    query_pos = 0
    INS_list = []
    DEL_list = []
    min_del_len = 45
    min_ins_len = 50

    for op_pair in op_pairs:
        op = op_pair[-1]
        length = int(op_pair[:-1])

        if op == "S":
            query_pos += length
        elif op in "M=X":
            ref_pos += length
            query_pos += length
        elif op == "I":
            if length > min_ins_len:
                INS_list.append((ref_pos, length))
            query_pos += length
        elif op == "N":
            ref_pos += length
        elif op == "D":
            if length > min_del_len:
                if is_reverse:
                    #DEL_list.append((ref_pos - length, length))
                    DEL_list.append((ref_pos, length))
                else:
                    DEL_list.append((ref_pos, length))
            ref_pos += length

    if op_pairs[0][-1] == "S":
        query_start = int(op_pairs[0][:-1]) + 1
    else:
        query_start = 0
    if op_pairs[-1][-1] == "S":
        query_end = query_pos - int(op_pairs[-1][:-1])
    else:
        query_end = query_pos
    
    return (ref_pos, query_pos, query_start, query_end, INS_list, DEL_list)


def get_alignments_from_samfile(samfile: pysam.AlignmentFile, max_nm: float
    ) -> tuple[dict[str, Alignment], dict[str, list[int, int, str, int]], dict[str, list[int, int, str, int]], dict[str, np.ndarray]]:
    '''Extract the reads, cigar insertions/deletions and coverage from the contents of the samfile.'''
    reads = {}
    all_INS = {}
    all_DEL = {}
    coverage = {}

    for query in samfile.fetch():
        chrom = query.reference_name
        if not chrom or not query.cigarstring or not query.query_name or query.mapping_quality < 30:
            continue

        seq = query.query_sequence
        ref_len, query_len, query_start, query_end, INS_list, DEL_list = _parse_cigar(cigar_string=query.cigarstring, is_reverse=query.is_reverse)

        ref_start = query.reference_start
        ref_end = query.reference_start + ref_len
        phase = query.get_tag("HP") if query.has_tag("HP") else 0

        INS_list = [(ins[0] + ref_start, ins[1], seq[ins[0]:ins[0]+ins[1]], phase) for ins in INS_list]
        DEL_list = [(d[0] + ref_start, d[1], "", phase) for d in DEL_list]
        all_INS.setdefault(chrom, []).extend(INS_list)
        all_DEL.setdefault(chrom, []).extend(DEL_list)

        # Increment the bases between the reference start and end
        coverage.setdefault(chrom, np.zeros((250000000)))[ref_start:ref_end+1] += 1
        # Remove one from the places of deletions (that are larger than the threshold)
        for del_start, del_end, _, _ in DEL_list:
            coverage[chrom][del_start:del_end+1] -= 1
        # Filter out alignments without supplimentary alignments
        if not query.has_tag("SA"):
            continue
        # Filter out alignments with an edit distance ratio above a threshold
        if query.has_tag("NM") and query.get_tag("NM") / query.reference_length > max_nm:
            continue
        # Mirror the start and endpoints of the alignment if on the reverse strand
        if query.is_reverse:
            query_start, query_end = query_len - query_end, query_len - query_start
            ref_start, ref_end = ref_end, ref_start
        
        alignment = Alignment(
            query_name= query.query_name,
            is_primary= not (query.is_secondary or query.is_supplementary),
            chrom= chrom,
            query_start= query_start,
            query_end= query_end,
            ref_start= ref_start+1,
            ref_end= ref_end,
            strand= "-" if query.is_reverse else "+",
            cigar= query.cigarstring,
            mapq= query.mapping_quality,
            phase= phase,
            seq= query.query_sequence
        )
        reads.setdefault(query.query_name, []).append(alignment)

    # Sort alignments by query start position
    for alignments in reads.values():
        alignments.sort(key=lambda a: a.query_start)
        
    return reads, all_INS, all_DEL, coverage