#!/usr/bin/env python

import sys, random
from GenomeParser import GenomeParser

SIZE_INDEL = 5

if len(sys.argv) < 4 or len(sys.argv) > 5:
	sys.stderr.write("Usage: {} genome.fasta variants.vcf output [seed]\n".format(sys.argv[0]))
	sys.exit(1)

def isweakspot(v):
	return (v[:2]*6) == v

rgfilename = sys.argv[1]
vcfname = sys.argv[2]
outputname = sys.argv[3]

if len(sys.argv) == 5:
	seed = sys.argv[4]
else:
	seed = None

random.seed(seed)

variantsWritten = {}

lpos = 0
with open(rgfilename) as rg, open(vcfname) as srcvcf, open(outputname, 'w') as trgvcf:
	gp = GenomeParser(rgfilename)
	rgLine = rg.readline()
	while(rgLine.startswith('>')):
		rgLine = rg.readline()
	rgLinelen = len(rgLine)

	for line in srcvcf:
		if line.startswith('#'):
			trgvcf.write(line)
		else:
			if int(random.random()*10) == 1:
				pos = int(line.split('\t')[1])
				if abs(pos-lpos) > 50:
					zone = gp.getNuc('{}-{}'.format(str(pos), str(pos+11)))
					if not isweakspot(zone) and pos not in variantsWritten:
						columns = line.split('\t')
						if len(columns[3]) <= SIZE_INDEL:
							alt = columns[4].split(',')
							good_alt = False
							for a in alt:
								if len(a) <= SIZE_INDEL:
									columns[4] = a
									good_alt = True
									break;
							if good_alt == True:
								newline = '\t'.join(columns)
								variantsWritten[pos] = True
								lpos = pos
								trgvcf.write(newline)

sys.stdout.write("{} variants written.\n".format(len(variantsWritten)))
