"""
Mapping utils used for mapping chromosome coordinates to genes or other gene elements.

Created on 5/4/2014
"""

from classes import *
from retrieve_annotations import *
import pygal
import os

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

def _mapping(indel_fn, annotation_fn='knownGenes.txt'):
	'''a draft function for crudely mapping indels to genes'''
	genes = load_annotations(annotation_fn)
	print 'Number of genes in the referrence genome: ', len(genes)
	## counters:
	a = 0 # indel count
	b = 0 # indels mapped to genes
	c = 0 # indels mapped to CDS
	d = 0 # indels mapped to exons
	genes_indel = {} # genes with insertion/deletion
	for indel in read_indel_iter(indel_fn):
		a += 1
		position = indel.position
		mapped_g = False # whether mapped to genes
		for gene in genes:
			if position in gene:
				mapped_g = True
				gene_name = gene.name
				genes_indel[gene_name] = [indel.typeStr, 'gene']
				if position in gene.CDS:
					c += 1
					genes_indel[gene_name][1] = 'CDS'
					for exon in gene.exons:
						if position in exon:
							genes_indel[gene_name][1] = 'exon'
							d += 1
							break
		if mapped_g:
			b += 1
	return genes_indel, a, b, c, d

def output_wraper(genes_indel, a, b, c, d):
	try:
		os.mkdir('output')
	except WindowsError:
		pass
	os.chdir('./output')
	plots(10,5,4,3)
	return


def plots(a, b, c, d):
	a = float(a)
	pie_chart = pygal.Pie()
	pie_chart.title = 'Indels mapped to genomic elements (in %)'
	pie_chart.add('Intergenic', 100*(a-b)/a)
	pie_chart.add('Gene (non-coding region)', 100*(b-c)/a)
	pie_chart.add('Exon', 100*d/a)
	pie_chart.add('Intron', 100*(c-d)/a)
	pie_chart.render_to_file('pie_chart.svg')
	return

def write_output(genes_indel, a, b, c, d): 
	return

## test:
output_wraper({},10,5,4,3)
