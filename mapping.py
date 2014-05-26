"""
Mapping utils used for mapping chromosome coordinates to genes or other gene elements.

Created on 5/4/2014
"""
from annotations import *
from classes import *

import os
import argparse as ap
import datetime as dt
from collections import Counter
import cookielib, urllib2
## dependency:
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import poster
import pygal as pg


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
	'''a  function for mapping indels to genes'''
	genes = load_annotations(annotation_fn)
	indels = []
	## counters:
	a = 0 # indel count
	b = 0 # indels mapped to genes
	c = 0 # indels mapped to CDS
	d = 0 # indels mapped to exons
	genes_indel = {} # genes with insertion/deletion
	insertion_lengths = []
	deletion_lengths = []
	print 'Mapping INDELs in %s to genome' % indel_fn
	for indel in read_indel_iter(indel_fn):
		a += 1
		if indel.chrom == 'chr23':
			indel.chrom = 'chrX'
		position = indel.position
		chrom = indel.chrom
		typeStr = indel.typeStr
		if typeStr == 'insertion':
			insertion_lengths.append(len(indel))
		elif typeStr == 'deletion':
			deletion_lengths.append(len(indel))

		indel.mapped_genes = []
		indel.mapped_cds = []
		indel.mapped_exons = []
		for gene in genes[chrom]:
			if position in gene:
				gene_name = gene.name
				indel.mapped_genes.append(str(gene))
				genes_indel[gene_name] = [typeStr, 'gene']
				if position in gene.CDS:
					indel.mapped_cds.append(str(gene))
					genes_indel[gene_name][1] = 'CDS'
					for exon in gene.exons:
						if position in exon:
							indel.mapped_exons.append(str(gene))
							genes_indel[gene_name][1] = 'exon'
							d += 1
							break
		if len(indel.mapped_genes) != 0:
			b += 1
		if len(indel.mapped_cds) != 0:
			c += 1
		indels.append(indel)
	return indels, genes_indel, a, b, c, d, insertion_lengths, deletion_lengths


def list_pprint(l):
	"""pretty print a list """
	if len(l) == 0:
		return 'NA'
	else:
		return ','.join(l)


def indels2bed(indels, indel_fn):
	"""convert a list of indels to bed format. """
	outfn = indel_fn[0:-6] + '.bed'
	header = ['chrom', 'chromStart', 'chromEnd', 'type', 'confScore',
		'mappedGenes', 'mappedCDS', 'mappedExons']
	with open (outfn, 'w') as out:
		out.write('#fields:\t')
		out.write('\t'.join(header) + '\n')
		for indel in indels:
			items = indel.position
			items.extend([indel.typeStr, indel.confidence, 
				list_pprint(indel.mapped_genes), list_pprint(indel.mapped_cds), list_pprint(indel.mapped_exons)])
			items = [str(i) for i in items]
			out.write('\t'.join(items) + '\n')
	return


def format_gene_sets(genes_indel):
	"""convert the genes_indel dictionary to a gene-set library format"""
	d_term_genes = {}
	for gene in genes_indel:
		term = '_'.join(genes_indel[gene])
		if term not in d_term_genes:
			d_term_genes[term] = [gene]
		else:
			d_term_genes[term].append(gene)
	return d_term_genes
		

def d2gmt(d, gmtfn):
	"""write a dictionary into gmt file format"""
	with open (gmtfn, 'w') as out:
		for term in d:
			out.write(term + '\tna\t')
			out.write('\t'.join(d[term]) + '\n')
	return


def enrichr_link(genesStr, meta=''):
	"""post a gene list to enrichr server and get the link."""
	cj = cookielib.CookieJar()
	opener = poster.streaminghttp.register_openers()
	opener.add_handler(urllib2.HTTPCookieProcessor(cookielib.CookieJar()))

	params = {'list':genesStr,'description':meta}
	datagen, headers = poster.encode.multipart_encode(params)
	url = "http://amp.pharm.mssm.edu/Enrichr/enrich"
	request = urllib2.Request(url, datagen, headers)
	urllib2.urlopen(request)

	x = urllib2.urlopen("http://amp.pharm.mssm.edu/Enrichr/share")
	response = x.read()
	split_strings = response.split('"')
	linkID = split_strings[3]
	share_url_head = "http://amp.pharm.mssm.edu/Enrichr/enrich?dataset="
	link = share_url_head + linkID
	return link


def genomic_distribution(indel_fn):
	"""get chromosome distribution of insertion and deletion"""
	print 'Calculating genomic distribution statistics'
	chroms_i = []
	chroms_d = []
	for indel in read_indel_iter(indel_fn):
		chrom = int(indel.chrom[3:]) # convert chrom to int
		if indel.typeStr == 'insertion':
			chroms_i.append(chrom)
		else:
			chroms_d.append(chrom)
		c_chroms_i = Counter(chroms_i)
		c_chroms_d = Counter(chroms_d)
	return c_chroms_i, c_chroms_d


def plot_pie(a, b, c, d):
	a = float(a)
	pie_chart = pg.Pie()
	pie_chart.title = 'Indels mapped to genomic elements (in %)'
	pie_chart.add('Intergenic', 100*(a-b)/a)
	pie_chart.add('Gene (non-coding region)', 100*(b-c)/a)
	pie_chart.add('Exon', 100*d/a)
	pie_chart.add('Intron', 100*(c-d)/a)
	pie_chart.render_to_file('indel_gene_element_mapping.svg')
	return


def plot_bars(c_chroms_i, c_chroms_d):
	"""bar plots for chromosome distribution of INDELs"""
	bar_plot = pg.Bar(spacing=0, x_title='Chromosomes', y_title='Number of INDELs')
	bar_plot.title = 'Chromosome distribution of INDELs'
	chroms = c_chroms_i.keys()
	chroms.sort()
	bar_plot.x_labels = map(str, chroms)
	bar_plot.add('Insertion', [c_chroms_i[x] for x in chroms])
	bar_plot.add('Deletion', [c_chroms_d[x] for x in chroms])
	bar_plot.add('INDELs', [c_chroms_i[x] + c_chroms_d[x] for x in chroms])
	bar_plot.render_to_file('chromosome_distribution.svg')
	return


def plot_hist(insertion_lengths, deletion_lengths):
	"""plot length distribution of INDELs with matplotlib"""
	fig = plt.figure(figsize=(12, 8))
	ax1 = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)
	
	insertion_lengths = np.array(insertion_lengths)
	deletion_lengths = np.array(deletion_lengths)
	indel_lengths = np.concatenate((insertion_lengths, deletion_lengths))
	min_len, max_len = indel_lengths.min(), indel_lengths.max() # get the range for lengths of INDELs
	bins = np.linspace(min_len, max_len, 20) # for histograms
	ax1.hist(insertion_lengths, bins=bins, stacked=True, label='insertion', color='b')
	ax1.hist(deletion_lengths, bins=bins, stacked=True, label='deletion', color='r')
	ax1.legend(loc='best')
	ax1.set_xlabel('Lengths', fontsize=18)
	ax1.set_ylabel('Counts', fontsize=18)

	x = np.linspace(min_len, max_len) # for PDF
	arrays = [insertion_lengths, deletion_lengths, indel_lengths]
	labels = ['insertion', 'deletion', 'INDEL']
	colors = ['b', 'r', 'm']
	for array, l, c in zip(arrays, labels, colors):
		kde = gaussian_kde(array)
		ax2.plot(x, kde(x), label=l, color=c, lw=2)
	ax2.legend(loc='best')
	ax2.yaxis.tick_right() 
	ax2.yaxis.set_label_position("right") # put y-label on the right
	ax2.set_xlabel('Lengths', fontsize=18)
	ax2.set_ylabel('Probability', fontsize=18)
	fig.tight_layout()
	plt.savefig('lengths_distributions.pdf')
	return


def write_output(indel_fn, d_term_genes, a, b, c, d):
	"""write the running infomation and results into a txt file"""
	outfn = indel_fn[0:-6] + '.out'
	with open (outfn, 'w') as out:
		## overall summary statistics:
		out.write('Summary statistics: \n')
		out.write('-'*8 + '\n')
		out.write('Total Number of INDELs: %s\n'%a)
		out.write('INDELs mapped to genes: %s\n'%b)
		out.write('INDELs mapped to CDS: %s\n'%c)
		out.write('INDELs mapped to exons: %s\n'%d)
		out.write('\n')
		## write Enrichr links:
		out.write('Enrichment analysis for INDELs affected genes: \n')
		out.write('-'*8 + '\n')
		print 'Performing enrichment analysis for the gene lists'
		for term in d_term_genes:
			link = enrichr_link('\n'.join(d_term_genes[term]), meta=term)
			out.write(term + ':\t' + link + '\n')
	return


def output_wraper(indel_fn, d_term_genes, a, b, c, d, insertion_lengths, deletion_lengths, indels):
	c_chroms_i, c_chroms_d = genomic_distribution(indel_fn)
	try:
		os.mkdir('output')
	except WindowsError:
		pass
	os.chdir('./output')
	write_output(indel_fn, d_term_genes, a, b, c, d)
	indels2bed(indels, indel_fn)
	d2gmt(d_term_genes, 'genes_affected_by_indels.gmt')
	plot_pie(a, b, c, d)
	plot_bars(c_chroms_i, c_chroms_d)
	plot_hist(insertion_lengths, deletion_lengths)
	pwd = os.getcwd()
	print "All results are stored in %s" % pwd
	return 


def main():
	## parsing argvs
	parser = ap.ArgumentParser()
	parser.add_argument('indel filename', help='the name of the indel file', type=str)
	parser.add_argument('-a','--annotation', help='the name of gene annotation file name', default='knownGenes.txt', type=str)
	args = vars(parser.parse_args())
	indel_fn = args['indel filename']
	annotation_fn = args['annotation']

	start_time = dt.datetime.now()
	print 'Thank you for using pyHgMapper'
	print start_time
	print '-'*10
	## map indels:
	indels, genes_indel, a, b, c, d, insertion_lengths, deletion_lengths = mapping(indel_fn, annotation_fn=annotation_fn)
	d_term_genes = format_gene_sets(genes_indel)
	## output results:
	output_wraper(indel_fn, d_term_genes, a, b, c, d, insertion_lengths, deletion_lengths, indels)
	end_time = dt.datetime.now()
	print '-'*10
	print end_time
	print 'Finished, time elapsed: ', end_time - start_time
	return


if __name__ == '__main__':
	main()
