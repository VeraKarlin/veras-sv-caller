import numpy as np
import statistics
from sklearn.cluster import DBSCAN
from data_structures import Alignment, BPInfo, SVInfo, SVType


def _add_group(sv_calls: dict[str, dict[str, list]], chrom: str, group: list, sv_type: SVType, coverage_dict: dict):
    window = 50
    starts = [int(e[0]) for e in group]
    if sv_type == SVType.INS:
        ends = [int(e[0]) for e in group] 
    elif sv_type == SVType.DEL:
        ends = [int(e[0]) + int(e[1]) for e in group] 

    seqs = [e[2] for e in group]
    seqs.sort(key=lambda s: len(s))
    phase_list = [0, 0, 0]
    for e in group:
        phase_list[int(e[3])] += 1

    start_pos = int(statistics.median(starts))
    first_cov_before = statistics.median(coverage_dict[chrom][start_pos-window:start_pos])
    first_cov_after = statistics.median(coverage_dict[chrom][start_pos:start_pos+window])
    first_bp = BPInfo(chrom=chrom, pos=start_pos, strand=0, direction=0, cov_before=first_cov_before, cov_after=first_cov_after)
    
    end_pos = int(statistics.median(ends))
    second_cov_before = statistics.median(coverage_dict[chrom][end_pos-window:end_pos])
    second_cov_after = statistics.median(coverage_dict[chrom][end_pos:end_pos+window])
    second_bp = BPInfo(chrom=chrom, pos=end_pos, strand=0, direction=0, cov_before=second_cov_before, cov_after=second_cov_after)
    
    phase = phase_list.index(max(phase_list))
    phase_ratio = max(phase_list) / sum(phase_list)
    sv_info = SVInfo(
        first=first_bp,
        second=second_bp,
        sv_type=sv_type, 
        sequence=seqs[len(seqs)//2], 
        support=len(group),
        phase=phase,
        phase_ratio=phase_ratio,
        sv_pipeline="cigar"
    )
    sv_calls.setdefault(chrom, {}).setdefault(chrom, []).append(sv_info)


def cluster_cigar_variants(cigar_var_dict: dict[str, list[tuple[int, int]]], split_var_dict: dict[str, list[tuple[int, int]]], 
                           sv_calls: dict[str, dict[str, SVInfo]], sv_type: SVType, coverage_dict: dict, eps: float, min_samples: int) -> dict[str, dict[str, SVInfo]]:
    '''Cluster the SVs which have been determined from the cigar strings with the SVs of that type that have been called from split reads.'''
    # Merge dictionaries:
    combined_var_dict = cigar_var_dict.copy()
    for key, val in split_var_dict.items():
        combined_var_dict.setdefault(key, []).extend(val)

    # Cluster the variants of each chromosome:
    for chrom in combined_var_dict.keys():
        variants = [variant[:2] for variant in combined_var_dict[chrom]]
        if not variants:
            continue
        X = np.array(variants)
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        labels = dbscan.fit(X).labels_
        unique_labels = set(labels)
        var_array = np.array(combined_var_dict[chrom])
        clusters = {int(label): var_array[labels == label].tolist() for label in unique_labels}
        for cluster in clusters.values():
            _add_group(sv_calls=sv_calls, chrom=chrom, group=cluster, sv_type=sv_type, coverage_dict=coverage_dict)
            
    return sv_calls