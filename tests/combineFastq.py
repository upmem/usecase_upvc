#!/usr/bin/env python

import sys, os, random

def printUsage():
	sys.stderr.write("Usage: {} set1prefix set2prefix outprefix [seed]\n".format(sys.argv[0]))
	sys.exit(1)

if len(sys.argv) < 4 or len(sys.argv) > 5:
	printUsage()

matprefix = sys.argv[1]
patprefix = sys.argv[2]
outprefix = sys.argv[3]
if len(sys.argv) == 5:
	seed = sys.argv[4]
else:
	seed = None

random.seed(seed)

mat_1 = "{}1.fq".format(matprefix)
if not os.path.exists(mat_1):
	mat_1 = "{}1.fastq".format(matprefix)
mat_2 = "{}2.fq".format(matprefix)
if not os.path.exists(mat_2):
	mat_2 = "{}2.fastq".format(matprefix)

pat_1 = "{}1.fq".format(patprefix)
if not os.path.exists(pat_1):
	pat_1 = "{}1.fastq".format(patprefix)

pat_2 = "{}2.fq".format(patprefix)
if not os.path.exists(pat_2):
	pat_2 = "{}2.fastq".format(patprefix)

out1 = "{}1.fastq".format(outprefix)
out2 = "{}2.fastq".format(outprefix)
counter = 0

if not os.path.exists(mat_1) or not os.path.exists(mat_2) or not os.path.exists(pat_1) or not os.path.exists(pat_2):
	printUsage()

with open(mat_1) as matfp1, open(mat_2) as matfp2, open(pat_1) as patfp1, open(pat_2) as patfp2, open(out1, 'w') as outfp1, open(out2, 'w') as outfp2:
	# get the header lines of the maternal and paternal pairs (starting with @)
	mat1buf = matfp1.readline()
	mat2buf = matfp2.readline()
	pat1buf = patfp1.readline()
	pat2buf = patfp2.readline()

	# while we still have sequences to add, keep iterating
	while (mat1buf != "" or pat1buf != ""):
		counter += 1
		# If we have gotten all entries of the paternal pair, or if we randomly choose the maternal pair, then use the maternal pair
		if (pat1buf == "") or (mat1buf != "" and random.random() < 0.5):
			mat1buf = matfp1.readline() + matfp1.readline() + matfp1.readline()
			mat2buf = matfp2.readline() + matfp2.readline() + matfp2.readline()
			outfp1.write("@{}\n{}".format(str(counter), mat1buf))
			outfp2.write("@{}\n{}".format(str(counter), mat2buf))
			mat1buf = matfp1.readline()
			matfp2.readline()
		else:
			pat1buf = patfp1.readline() + patfp1.readline() + patfp1.readline()
			pat2buf = patfp2.readline() + patfp2.readline() + patfp2.readline()
			outfp1.write("@{}\n{}".format(str(counter), pat1buf))
			outfp2.write("@{}\n{}".format(str(counter), pat2buf))
			pat1buf = patfp1.readline()
			patfp2.readline()


sys.stdout.write("Finished processing {} reads\n".format(str(counter)))
