"""
Scripts to retrieve human genome annotation files and process them into desirable formats

Created on 4/27/2014
"""

__author__ = """\n""".join(['Zichen Wang <zichen.wang@mssm.edu>'])

import os
import urllib2, gzip
import cPickle as pickle
from classes import Gene

def retrieve_annotation(url=None):
	'''function for retrieving UCSC gene table, not used'''
	os.chdir('./annotations')
	if url is None:
		## url for UCSC knownGene table
		url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz'
	try:
		u = urllib2.urlopen(url)
	except Exception as e:
		print e, "\ncouldn't download the genome annotation file, please check the url!"
	else: ## if no network errors occur
		meta = u.info()
		fn = url.split('/')[-1]
		if fn in os.listdir(os.getcwd()):
			print 'The file %s is already downloaded!' % fn
			return
		else: ## file is not downloaded yet
			file_size = int(meta["Content-Length"])
			file_size_dl = 0
			block_size = 8192
			print "Downloading: %s Bytes: %s from UCSC" % (fn, file_size)
			with open (fn, 'wb') as out:
				while True:
					buffer = u.read(block_size)
					if not buffer:
					    break
					file_size_dl += len(buffer)
					out.write(buffer)
			print 'Finished downloading file: ', fn
	return


def parse_gene_table(fn):
	os.chdir('./annotations')
	genes = []
	if fn in os.listdir(os.getcwd()):
		print 'Reading information of genes from %s'%fn
		if fn.endswith('.gz'):
			with gzip.open(fn, 'rb') as f:
				for line in f:
					if not line.startswith('#'):
						G = Gene(line)
						genes.append(G)
		else:
			with open (fn,'r') as f:
				for line in f:
					if not line.startswith('#'):
						G = Gene(line)
						genes.append(G)
		print 'Finished reading gene table'
		print 'Pickling genes data'
		pfn = fn.split('.')[0] + '.p'
		pickle.dump(genes, open(pfn, 'wb'))
		os.chdir('../') # cd back to the working dir
		print 'Finished '
		return genes
	else:
		raise IOError('File %s not found'%fn)


def load_annotations(fn):
	os.chdir('./annotations')
	pfn = fn.split('.')[0] + '.p'
	if pfn in os.listdir(os.getcwd()):
		print 'Loading annotations from %s'%pfn
		genes = pickle.load(open(pfn, 'rb'))
		os.chdir('../')
	else:
		print 'Loading annotations from %s'%fn
		genes = parse_gene_table(fn)
	return genes

