import argparse

def make_copy_without_sv(vcf_path: str, excluded_SV_string: str, output_path: str) -> None:
    excluded_SV_set = {sv.upper() for sv in excluded_SV_string.split(",")}
    with open(vcf_path, "r") as file:
        with open(output_path, "w") as out:
            for line in file.readlines():
                if line.startswith("#"):
                    out.write(line)
                    continue
                var = line.split("SVTYPE=")[-1].split(";")[0].split(":")[0]
                if var in excluded_SV_set:
                    continue
                out.write(line)
        out.close
    file.close()


def main():
    parser = argparse.ArgumentParser(
        prog='MakeCopyWithoutSV',
        description='Reads a VCF file and removes the desired structural variants.'
    )
    parser.add_argument('--in', help='path to input vcf file')
    parser.add_argument('--svs', help='names of SVs to be excluded, separated by commas (e.g. INS,DEL,BND,DUP)')
    parser.add_argument('--out', help='path to filtered output vcf file')
    args = vars(parser.parse_args())
    make_copy_without_sv(args["in"], args["svs"], args["out"])
    print(f"Filtered VCF file {args["out"]} created!")


if __name__ == "__main__":
    main()