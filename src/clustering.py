import numpy as np
from sklearn.cluster import DBSCAN

from data_structures import Alignment
from dataclasses import replace


def _pair_alignments(reads: dict[str, list[Alignment]]) -> list[tuple[Alignment, Alignment]]:
    '''Pair adjacent aligmments within a read.'''
    pairs = []
    for name, alignments in reads.items():
        for i in range(len(alignments) - 1):
            pairs.append((alignments[i], alignments[i+1]))
    return pairs


def _flip_alignments(pairs: list[tuple[Alignment, Alignment]]) -> list[tuple[Alignment, Alignment]]:
    '''Reorder alignments of each pair based on chromosome number and position in read.'''
    new_list = []
    for orig_align_1, orig_align_2 in pairs:
        align_1 = replace(orig_align_1)
        align_2 = replace(orig_align_2)
        if (align_1.chrom == align_2.chrom and align_1.ref_end > align_2.ref_start) or align_1.chrom > align_2.chrom:
            align_1.ref_start, align_1.ref_end = align_1.ref_end, align_1.ref_start
            align_2.ref_start, align_2.ref_end = align_2.ref_end, align_2.ref_start
            new_list.append((align_2, align_1))
        else:
            new_list.append((align_1, align_2))
    return new_list


def _group_pairs_by_chrom(pairs: list[tuple[Alignment, Alignment]]) -> dict[str, dict[str, tuple[Alignment, Alignment]]]:
    '''Group pairs by chromosomes and store them in a dictionary.'''
    grouped_pairs = {}
    for pair in pairs:
        grouped_pairs.setdefault(pair[0].chrom, {}).setdefault(pair[1].chrom, []).append(pair)
    return grouped_pairs


def cluster_positions(reads: dict[str, list[Alignment]], eps:float, min_samples:int) -> dict[str, dict[str, dict[int, tuple[Alignment, Alignment]]]]:
    '''Cluster the alignments with DBSCAN into clusters of alignment pairs.'''
    all_pairs = _pair_alignments(reads)
    flipped_pairs = _flip_alignments(all_pairs)
    grouped_pairs = _group_pairs_by_chrom(flipped_pairs)

    grouped_clusters = {}
    for chrom_1 in grouped_pairs:
        for chrom_2 in grouped_pairs[chrom_1]:
            pairs = grouped_pairs[chrom_1][chrom_2]
            X = np.array([(pair[0].ref_end, pair[1].ref_start) for pair in pairs])
            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            labels = dbscan.fit(X).labels_
            unique_labels = set(labels)
            if not unique_labels:
                continue

            # Create a dict where the labels are the cluster indexes and the values the Alignment pairs
            clusters = {int(label): np.array(pairs)[labels == label].tolist() for label in unique_labels}
            grouped_clusters.setdefault(chrom_1, {})[chrom_2] = clusters

    return grouped_clusters