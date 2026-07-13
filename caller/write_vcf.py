import pysam
from data_structures import SVType, SVInfo


def generate_VCF_header(file, contig_info):
	# General header
	file.write("##fileformat=VCFv4.2\n")
	file.write("##source=test_SV_caller")
	import time
	file.write("##fileDate=%s\n"%(time.strftime('%Y-%m-%d %H:%M:%S %w-%Z',time.localtime())))
	for contig in contig_info:
		file.write("##contig=<ID=%s,length=%d>\n"%(contig[0], contig[1]))

	# Specific header
	# ALT
	file.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
	file.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
	file.write("##ALT=<ID=DUP:TANDEM,Description=\"Tandem duplication\">\n")
	file.write("##ALT=<ID=DUP:INV,Description=\"Inverted tandem duplication\">\n")
	file.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
	file.write("##ALT=<ID=INS,Description=\"Insertion\">\n")
	file.write("##ALT=<ID=BND,Description=\"Break end of translocation\">\n")

	# INFO
	file.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variant\">\n")
	file.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n")
	file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
	file.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
	file.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n")
	file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
	file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n")
	file.write("##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"Confidence interval around inserted/deleted material between breakends\">\n")
	# file.write("##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n")
	file.write("##INFO=<ID=RE,Number=1,Type=Integer,Description=\"Number of read support this record\">\n")
	file.write("##INFO=<ID=STRAND,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n")
	file.write("##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Supporting read names of SVs (comma separated)\">\n")
	file.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency.\">\n")
	file.write("##FILTER=<ID=q5,Description=\"Quality below 5\">\n")
	# FORMAT
	file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	file.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# High-quality reference reads\">\n")
	file.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# High-quality variant reads\">\n")
	file.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">\n")
	file.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"# Genotype quality\">\n")

	#file.write("##CommandLine=\"cuteSV %s\"\n"%(" ".join(argv)))

def write_sv_lines(vcf_file, sv_calls: dict[str, dict[str, list[SVInfo]]]):
	lines = []
	i = 0
	for first_chrom in sv_calls.keys():
		for second_chrom in sv_calls[first_chrom].keys():
			calls = list(sv_calls[first_chrom][second_chrom])
			calls.sort(key=lambda x: int(x.first.pos))

			for call in calls:
				if call.first.cov_before + call.second.cov_after < 2:
					continue
				cov_ratio = 2 * call.support / (call.first.cov_before + call.second.cov_after)
				if cov_ratio < 0.25:
					continue
				line = []

				# CHROM  
				#line.append(call.first.chrom)
				line.append(first_chrom)

				# POS  
				line.append(str(call.first.pos))

				# ID  
				line.append(f"sv_caller.{call.sv_type.name}.{i}.{call.sv_pipeline}")
				i += 1

				# REF  
				line.append("N")

				# ALT
				if call.sv_type == SVType.INS:
					line.append(call.sequence)
				elif call.sv_type == SVType.BND:
					if call.sv_type == SVType.INV:
						line.append(f"[{second_chrom}:{call.second.pos}[N")
					else:
						line.append(f"N]{second_chrom}:{call.second.pos}]")
				else:
					line.append(f"<{call.sv_type.name}>")

				# QUAL  
				line.append("60")

				# FILTER
				line.append("PASS")

				# INFO  
				info = []
				info.append(f"SVTYPE={call.sv_type.name}")
				if call.sv_type == SVType.DEL:
					# The length should be negative if the SV is a deletion
					info.append(f"SVLEN={call.first.pos - call.second.pos - 1}")
				elif call.sv_type == SVType.INS:
					info.append(f"SVLEN={len(call.sequence)}")
				elif call.sv_type == SVType.BND:
					pass
				elif call.sv_type == SVType.INV or call.sv_type == SVType.DUP:
					info.append(f"SVLEN={call.second.pos - call.first.pos}")
				if call.sv_type == SVType.BND:
					info.append(f"END={call.first.pos}")
				else:
					info.append(f"END={call.second.pos}")
				info.append(f"SUPPORT={call.support}")
				info.append(f"COVERAGE={call.first.cov_before},{call.first.cov_after},{call.second.cov_before},{call.second.cov_after}")
				if call.sv_type == SVType.BND:
					info.append(f"CHR2='{second_chrom}'")
				info.append(f"PHASE={call.phase}")
				
				info.append(f"VAF={cov_ratio:.3f}")
				line.append(";".join(info))

				# FORMAT
				line.append("GT:GQ:PL:AD")

				# SAMPLE
				sample = []
				#if call.phase_ratio < 0.75:
				if cov_ratio > 0.75:
					sample.append("1|1")
				else:
					sample.append("0|1")
				sample.append("60")
				sample.append(str(round(call.phase_ratio, 2)))
				line.append(":".join(sample))

				# Combine all lines
				lines.append("\t".join(line))

	vcf_file.write("\n".join(lines))


def write_vcf_file(samfile: pysam.Samfile, sv_calls: dict, output_path: str):
	chr_name_list = list()
	contig_info = list()
	ref_ = samfile.get_index_statistics()
	for i in ref_:
		chr_name_list.append(i[0])
		local_ref_len = samfile.get_reference_length(i[0])
		contig_info.append([i[0], local_ref_len])

	with open(output_path, 'w') as vcf_file:
		generate_VCF_header(file=vcf_file, contig_info=contig_info)
		vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
		write_sv_lines(vcf_file=vcf_file, sv_calls=sv_calls)
		vcf_file.close()