"""
Mapping utils used for mapping chromosome coordinates to genes or other gene elements.

Created on 5/4/2014
"""

from classes import *
from retrieve_annotations import *


def read_indel(fn):
	'''read indels from .indel file into a list of Indels'''
	indels = []
	with open (fn) as f:
		for line in f:
			if not line.startswith('#'):
				sl = line.strip().split('\t')
				chrom = 'chr'+sl[2]
				refStart = int(float(sl[6]))
				refEnd = int(float(sl[7]))
				confidence = float(sl[9])
				typeStr = sl[10]
				coordinates = [chrom, refStart, refEnd]
				indel = Indel(coordinates, confidence, typeStr)
				indels.append(indel)
	return indels

def read_indel_iter(fn):
	"""yield generator for Indels if the file is large"""
	indels = []
	with open (fn) as f:
		for line in f:
			if not line.startswith('#'):
				sl = line.strip().split('\t')
				chrom = 'chr'+sl[2]
				refStart = int(float(sl[6]))
				refEnd = int(float(sl[7]))
				confidence = float(sl[9])
				typeStr = sl[10]
				coordinates = [chrom, refStart, refEnd]
				indel = Indel(coordinates, confidence, typeStr)
				yield indel

def mapping(indel_fn, annotation_fn='knownGenes.txt'):
	'''a draft function for crudely mapping indels to genes'''
	genes = load_annotations(annotation_fn)
	print 'Number of genes in the referrence genome: ',len(genes)
	for indel in read_indel_iter(indel_fn):
		for gene in genes:
			if indel.position in gene:
				print 'indel position:%s, mapped to gene %s'%(indel.position, gene.refSeqId)
				yield (indel, gene)

## test:
c = 0
for indel, gene in mapping('NA12878_to_hg19.indel'):
	c += 1
print 'total # of mapped indel-gene pairs:', c
