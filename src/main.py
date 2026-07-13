import os
import pysam
import typer
from typing import Annotated
from pathlib import Path

import bam_parser
import clustering
import identify_sv
import cigar_variants
import write_vcf
from src.data_structures import SVType


app = typer.Typer()
@app.command()
def main(
    bam: Annotated[Path, typer.Option(help="Path to input bam file.")],
    output: Annotated[str, typer.Option(help="Path to output vcf file.")] = None,
    max_nm: Annotated[float, typer.Option(help="Maximum threshold for ratio between edit distance and alignment length.")] = 0.05,
    eps: Annotated[float, typer.Option(help="Epsilon for DBSCAN clustering.")] = 50,
    min_samples: Annotated[int, typer.Option(help="Minumum number of samples for DBSCAN clustering.")] = 3,
):  
    try:
        samfile = pysam.AlignmentFile(str(bam), 'rb')
    except OSError as e:
        print(f"ERROR: {e}")
        raise typer.Exit(code=1)
    print(f"\nSamfile loaded!")

    ## Parse bam file ##
    reads, INS_dict, DEL_dict, coverage_dict = bam_parser.get_alignments_from_samfile(samfile=samfile, max_nm=max_nm)
    print(f"Samfile parsed!")

    ## Cluster split reads ##
    grouped_clusters = clustering.cluster_positions(reads, eps=eps, min_samples=min_samples)
    print(f"Clusters clustered!")
    
    ## Identify SV from split reads ##
    sv_calls, split_INS_dict, split_DEL_dict = identify_sv.get_cluster_info(cluster_dict=grouped_clusters, coverage_dict=coverage_dict, samfile=samfile)
    print(f"SVs called!")

    ## Add SV calls from cigar evidence ##
    sv_calls = cigar_variants.cluster_cigar_variants(cigar_var_dict=INS_dict, split_var_dict=split_INS_dict, sv_calls=sv_calls, 
                                                     sv_type=SVType.INS, coverage_dict=coverage_dict, eps=eps, min_samples=min_samples)
    print(f"Insertions done!")
    sv_calls = cigar_variants.cluster_cigar_variants(cigar_var_dict=DEL_dict, split_var_dict=split_DEL_dict, sv_calls=sv_calls, 
                                                     sv_type=SVType.DEL, coverage_dict=coverage_dict, eps=eps, min_samples=min_samples)
    print(f"Deletions done!")

    ## Write VCF file ##
    sample_id = ".".join(os.path.basename(bam).split(".")[:-1])
    if output == None:
        output = f"{sample_id}.veras_sv_caller.vcf"
    write_vcf.write_vcf_file(samfile=samfile, sv_calls=sv_calls, output_path=output)
    samfile.close() 
    print(f"VCF file written!\n")


if __name__ == "__main__":
    app()
