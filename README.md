veras-sv-caller is a structural variant caller for long-read sequencing.

## Usage
    python src/main.py --bam input.bam
Additional options:
    --output       name of output VCF file (default: path of input file)
    --max-nm       maximum allowed ratio between alignment copy distance and length (default: 0.05)
    --eps          eps parameter for DBSCAN (default: 50)
    --min-samples  min_samples parameter for DBSCAN (default: 3)

## Requirements
* Python ==3.14.6
* pysam >=0.24.0
* numpy >=2.5.1
* scikit-learn >=1.9.0
* typer >=0.26.8 
