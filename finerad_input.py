#!/usr/bin/env python2
# -*- coding: utf-8 -*-



'''
The program takes as input the .haplotypes.tsv file from Stacks, the.alleles.loci
from ipyrad, or the .alleles from pyrad to parse the haplotypes (phased SNPs) of
each sample in the dataset and convert it into a haplotype matrix ready for 
analysis with fineRADstructure (http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html)

fineRADstructure is a more sensistive method than Structure or Admixture because uses
the entire haplotype information instead of a single randomly sampled SNP

Type python finerad_input.py -h to show the help.
'''


__author__      = "Edgardo M. Ortiz"
__version__     = "1.0"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2017-10-23"


import sys
import argparse


def stacks_to_finerad(filename, outfile, minsample):
	output = open(outfile, "w")
	locus_num = 0
	with open(filename) as stacks_haps:
		for line in stacks_haps:
			# Write header of output file with only sample names
			if line.startswith("Catalog"):
				header = line.strip("\n").split("\t")[2:]
				output.write("\t".join(header)+"\n")

			# Skip locus if it is invariable
			elif "consensus" in line:
				continue

			# Process if locus is variable
			else:
				genotypes = line.strip("\n").replace("-","").split("\t")[2:]

				# Skip locus if any sample has more than two alleles
				if any(g.count("/") > 1 for g in genotypes):
					continue
				else:

					# Write row of haplotypes to output if locus contains enough samples
					if len(genotypes) - genotypes.count("") >= minsample:
						output.write("\t".join(genotypes)+"\n")
						locus_num += 1
	output.close()
	print str(locus_num)+" loci for "+str(len(genotypes))+" samples written to "+outfile


def pyrad_to_finerad(filename, outfile, minsample):
	haplotypes = {}

	# First pass to the file to get all sample names
	with open(filename) as pyrad_alleles:
		for line in pyrad_alleles:
			if not line.startswith("/"):
				sample = line.split()[0].strip(">")[:-2] # Sample name does not include last two characters
				if sample not in haplotypes:
					haplotypes[sample] = []

	# Second pass to the file to parse haplotypes
	with open(filename) as pyrad_alleles:
		locus = []
		locus_num = 0
		for line in pyrad_alleles:
			if not line.startswith("/"):
				locus.append(line.strip("\n"))
			else:

				# Skip if locus doesn't contain enough samples
				if len(locus)/2 >= minsample: # Every sample has two lines
					snp_raw_coords = []
					for pos in range(len(line)):
						if line[pos] in ["*","-"]:
							snp_raw_coords.append(pos)

					# Skip locus if it is invariable
					if snp_raw_coords == []:
						locus = []

					# Process if locus is variable
					else:

						# But first check that SNPs (columns) don't contain deletions, or more than 50% Ns
						snp_coords = []
						for pos in snp_raw_coords:
							snp = ""
							for allele in locus:
								snp += allele[pos]
							if snp.count("-") == 0:
								if snp.count("N") <= len(locus)/4:
									snp_coords.append(pos)

						# Add an empty cell to each sample in haplotypes
						for sample in haplotypes:
							haplotypes[sample].append("")

						# Fill the last empty cells with the haplotypes
						for allele in locus:
							sample = allele.split()[0].strip(">")[:-2]
							hap = ""
							hap = "".join([hap+allele[c] for c in snp_coords])

							# Don't add haplotype if it is more than 50% Ns
							if hap.count("N") >= len(hap)/2:
								hap = ""

							# Add haplotypes to samples
							if haplotypes[sample][-1] == "":
								haplotypes[sample][-1] = hap
							else:
								if haplotypes[sample][-1] != hap:
									haplotypes[sample][-1] = haplotypes[sample][-1]+"/"+hap
						locus_num += 1
						locus = []
				else:
					locus = []

	# Remove samples without data
	for sample in haplotypes.keys():
		if haplotypes[sample].count("") == len(haplotypes[sample]):
			del haplotypes[sample]

	# Write header
	output = open(outfile, "w")
	output.write("\t".join(sorted(haplotypes.keys()))+"\n")

	# Write genotypes
	final_locus_num = 0
	for l in range(locus_num):
		genotypes = []
		for sample in sorted(haplotypes.keys()):
			genotypes.append(haplotypes[sample][l])

		# Final check for minsample
		if len(genotypes) - genotypes.count("") >= minsample:
			output.write("\t".join(genotypes)+"\n")
			final_locus_num += 1
	output.close()
	print str(final_locus_num)+" loci for "+str(len(genotypes))+" samples written to "+outfile


def main():
	parser = argparse.ArgumentParser(description="Converts haplotype data from Stacks, pyrad or ipyrad into a haplotype matrix for analysis with fineRADstructure (http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html)")
	parser.add_argument("-i", "--input", action="store", dest="filename", required=True,
		help="Name of the haplotypes.tsv file from Stacks, the .alleles file from pyrad, or the .alleles.loci file from ipyrad to be converted.")
	parser.add_argument("-t", "--type", action="store", dest="data_type", default="",
		help="Source of the input file, options: stacks, pyrad, ipyrad, default=guess")
	parser.add_argument("-n", "--minsample", action="store", type=int, default=2,
		help="Minimum number of samples in a locus, default=2")
	parser.add_argument("-o", "--output", action="store", dest="outfile", default="",
		help="Name for the output file, default=input name + .minsample + .finerad")
	args = parser.parse_args()

	filename = args.filename
	data_type = args.data_type
	minsample = args.minsample
	if args.outfile == "":
		outfile = filename+".min"+str(minsample)+".finerad"
	else:
		outfile = args.outfile

	if "pyrad" in data_type.lower() or "alleles" in filename:
		pyrad_to_finerad(filename, outfile, minsample)
	elif "stacks" in data_type.lower() or "haplotypes.tsv" in filename:
		stacks_to_finerad(filename, outfile, minsample)
	else:
		print "Unknown input format, use -h for help"

if __name__ == "__main__":
    main()
