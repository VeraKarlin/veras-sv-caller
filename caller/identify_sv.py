from data_structures import Alignment, SVInfo, BPInfo, SVType
from statistics import median, mean
import pysam


def check_coverage(chrom: str, pos: int, cov_dict: dict) -> tuple[float, float]:
    window = 50
    cov_before = median(cov_dict[chrom][pos-window:pos]) if cov_dict[chrom][pos-window:pos].any() else 0
    cov_after = median(cov_dict[chrom][pos:pos+window]) if cov_dict[chrom][pos:pos+window].any() else 0
    return (cov_before, cov_after)


def get_cluster_info(samfile: pysam.Samfile, cluster_dict: dict[str, dict[str, dict[int, list[tuple[Alignment, Alignment]]]]], coverage_dict: dict
                     ) -> dict[str, dict[str, list[SVInfo]]]:
    sv_calls = {}
    for first_chrom in cluster_dict.keys():
        sv_calls[first_chrom] = {}
        for second_chrom in cluster_dict[first_chrom].keys():
            clusters = list(cluster_dict[first_chrom][second_chrom].items())
            clusters.sort(key=lambda p: mean((p[1][0][0].ref_end, p[1][0][1].ref_start)))

            sv_calls[first_chrom][second_chrom] = []
            for label, pairs in clusters:
                if label == -1:
                    continue
                same_strand = 0
                all_pos_1 = []
                all_pos_2 = []
                # The following vairables ranges from -1 to 1, with 0 representing a total mix
                strand_1 = 0
                strand_2 = 0
                tot_direction_1 = 0
                tot_direction_2 = 0
                phase_list = [0, 0, 0]
                
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
                first_bp_cov = check_coverage(first_chrom, first_bp_pos, coverage_dict)
                first_bp = BPInfo(
                    chrom=first_chrom,
                    pos=first_bp_pos,
                    strand=strand_1 / n,
                    direction=direction_1,
                    cov_before=first_bp_cov[0],
                    cov_after=first_bp_cov[1]
                )

                second_bp_pos = int(median(all_pos_2))
                second_bp_cov = check_coverage(second_chrom, second_bp_pos, coverage_dict)
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

                sv_set = set()
                if same_strand_ratio <= 0.25:
                    sv_set.add(SVType.INV)
                    # Check if inversion with deletions:
                    if len(sv_calls[first_chrom][second_chrom]) > 0:
                        last_call = sv_calls[first_chrom][second_chrom][-1] 
                        if SVType.INV in last_call.structural_variants:
                            bp_1 = last_call.first
                            bp_2 = first_bp
                            bp_3 = last_call.second
                            bp_4 = second_bp
                            # If the inversions are overlapping:
                            if bp_2.pos < bp_3.pos:
                                # Remove the previous inversion call
                                sv_calls[first_chrom][second_chrom].pop()
                                if abs(bp_2.pos - bp_1.pos) > 50:
                                    sv_calls[first_chrom][second_chrom].append(SVInfo(
                                        first=bp_1, 
                                        second=bp_2, 
                                        structural_variants={SVType.DEL}, 
                                        sequence="", 
                                        support=len(pairs), 
                                        phase=phase, 
                                        phase_ratio=phase_ratio, 
                                        sv_pipeline="soft_")
                                    )
                                sv_calls[first_chrom][second_chrom].append(SVInfo(
                                    first=bp_2, 
                                    second=bp_3, 
                                    structural_variants={SVType.INV}, 
                                    sequence="", 
                                    support=len(pairs), 
                                    phase=phase, 
                                    phase_ratio=phase_ratio, 
                                    sv_pipeline="soft_")
                                )
                                if abs(bp_4.pos - bp_3.pos) > 50:
                                    sv_calls[first_chrom][second_chrom].append(SVInfo(
                                        first=bp_3, 
                                        second=bp_4, 
                                        structural_variants={SVType.DEL}, 
                                        sequence="", 
                                        support=len(pairs), 
                                        phase=phase, 
                                        phase_ratio=phase_ratio, 
                                        sv_pipeline="soft_")
                                    )
                                continue
                
                if first_chrom != second_chrom:
                    sv_set.add(SVType.BND)

                    sv_info = SVInfo(
                        first=second_bp, 
                        second=first_bp, 
                        structural_variants=sv_set, 
                        sequence="<INS>", 
                        support=len(pairs), 
                        phase=phase, 
                        phase_ratio=phase_ratio,
                        sv_pipeline="soft_comp_bnd_"
                    )
                    sv_calls.setdefault(second_chrom, {}).setdefault(first_chrom, []).append(sv_info)

                elif direction_1 >= 0.8 and direction_2 <= -0.8:
                    if abs(second_bp.pos - first_bp.pos) >= 50:
                        sv_set.add(SVType.DUP)
                    else:
                        sv_set.add(SVType.INS)
                elif direction_1 <= -0.8 and direction_2 >= 0.8:
                    sv_set.add(SVType.DEL)
                elif same_strand_ratio > 0.25:
                    sv_set.add(SVType.INS)
                # TODO: Add sequence if insertion
                sv_info = SVInfo(
                    first=first_bp, 
                    second=second_bp, 
                    structural_variants=sv_set, 
                    sequence="<INS>", 
                    support=len(pairs), 
                    phase=phase, 
                    phase_ratio=phase_ratio,
                    sv_pipeline="soft_"
                    )
                sv_calls[first_chrom][second_chrom].append(sv_info)
    return sv_calls
