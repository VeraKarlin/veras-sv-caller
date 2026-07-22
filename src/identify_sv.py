from statistics import median, mean
import pysam
from data_structures import Alignment, SVInfo, BPInfo, SVType


def _check_coverage(chrom: str, pos: int, cov_dict: dict) -> tuple[float, float]:
    window = 50
    cov_before = median(cov_dict[chrom][pos-window:pos]) if cov_dict[chrom][pos-window:pos].any() else 0
    cov_after = median(cov_dict[chrom][pos:pos+window]) if cov_dict[chrom][pos:pos+window].any() else 0
    return (cov_before, cov_after)


def get_cluster_info(cluster_dict: dict[str, dict[str, dict[int, list[tuple[Alignment, Alignment]]]]], coverage_dict: dict,
                     ) -> dict[str, dict[str, list[SVInfo]]]:
    '''Determine SV type and properties of each cluster. Insertions and deletions are stored separately from the other SV calls.'''
    sv_calls = {}
    del_dict = {}
    ins_dict = {}
    for first_chrom in cluster_dict.keys():
        for second_chrom in cluster_dict[first_chrom].keys():
            clusters = list(cluster_dict[first_chrom][second_chrom].items())
            clusters.sort(key=lambda p: mean((p[1][0][0].ref_end, p[1][0][1].ref_start)))
            for label, pairs in clusters:
                if label == -1:
                    continue
                same_strand = 0
                all_pos_1 = []
                all_pos_2 = []
                phase_list = [0, 0, 0]
                # The following vairables range from -1 to 1, with 0 representing a total mix
                strand_1 = 0
                strand_2 = 0
                tot_direction_1 = 0
                tot_direction_2 = 0
                
                for align_1, align_2 in pairs:
                    all_pos_1.append(align_1.ref_end)
                    all_pos_2.append(align_2.ref_start)
                    # Measure the average strand
                    strand_1 += 1 if align_1.strand == "+" else - 1
                    strand_2 += 1 if align_2.strand == "+" else - 1
                    # Measure if the strands are the same
                    if align_1.strand == align_2.strand:
                        same_strand += 1
                    # Measure the directions of the aligned segment
                    tot_direction_1 += 1 if align_1.ref_start > align_1.ref_end else -1
                    tot_direction_2 += 1 if align_2.ref_start < align_2.ref_end else -1
                    # Check phase
                    phase_list[align_1.phase] += 1
                    phase_list[align_2.phase] += 1

                n = len(pairs)
                direction_1 = tot_direction_1 / n
                direction_2 = tot_direction_2 / n
                same_strand_ratio = same_strand / n

                first_bp_pos = int(median(all_pos_1))
                first_bp_cov = _check_coverage(first_chrom, first_bp_pos, coverage_dict)
                first_bp = BPInfo(
                    chrom=first_chrom,
                    pos=first_bp_pos,
                    strand=strand_1 / n,
                    direction=direction_1,
                    cov_before=first_bp_cov[0],
                    cov_after=first_bp_cov[1]
                )

                second_bp_pos = int(median(all_pos_2))
                second_bp_cov = _check_coverage(second_chrom, second_bp_pos, coverage_dict)
                second_bp = BPInfo(
                    chrom=second_chrom,
                    pos=second_bp_pos,
                    strand=strand_2 / n,
                    direction=direction_2,
                    cov_before=second_bp_cov[0],
                    cov_after=second_bp_cov[1]
                )

                phase = phase_list.index(max(phase_list))
                phase_ratio = max(phase_list) / sum(phase_list)
                is_inverted = same_strand_ratio <= 0.25
                sv_type = None

                if first_chrom != second_chrom:
                    sv_type = SVType.BND
                    sv_info = SVInfo(
                        first=second_bp, 
                        second=first_bp, 
                        sv_type=sv_type, 
                        sequence="<INS>", 
                        support=len(pairs), 
                        phase=phase, 
                        phase_ratio=phase_ratio,
                        sv_pipeline="soft_comp_bnd_",
                        is_inverted=is_inverted
                    )
                    sv_calls.setdefault(second_chrom, {}).setdefault(first_chrom, []).append(sv_info)
                elif is_inverted:
                    sv_type = SVType.INV
                elif direction_1 >= 0.8 and direction_2 <= -0.8:
                    if abs(second_bp.pos - first_bp.pos) >= 50:
                        sv_type = SVType.DUP
                    else:
                        sv_type = SVType.INS
                elif direction_1 <= -0.8 and direction_2 >= 0.8:
                    sv_type = SVType.DEL
                elif same_strand_ratio > 0.25:
                    sv_type = SVType.INS
                else:
                    print(f"No SV type identifed for SV at {first_chrom}:{first_bp_pos}")
                    continue

                if sv_type == SVType.INS:
                    for align_1, align_2 in pairs:
                        ins_seq = align_1.seq[align_1.query_end:align_2.query_start]
                        # Deal with "insertions" with a length of 0
                        ins_dict.setdefault(first_chrom, []).append((align_1.ref_end, len(ins_seq), ins_seq, phase))
                    continue
                if sv_type == SVType.DEL:
                    for align_1, align_2 in pairs:
                        del_dict.setdefault(first_chrom, []).append((align_1.ref_end, align_2.ref_start - align_1.ref_end + 1, "", phase))
                    continue

                sv_info = SVInfo(
                    first=first_bp,
                    second=second_bp, 
                    sv_type=sv_type, 
                    sequence="<INS>", 
                    support=len(pairs), 
                    phase=phase, 
                    phase_ratio=phase_ratio,
                    sv_pipeline="soft_",
                    is_inverted=is_inverted
                    )
                sv_calls.setdefault(first_chrom, {}).setdefault(second_chrom, []).append(sv_info)

    return sv_calls, ins_dict, del_dict
