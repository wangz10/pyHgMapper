"""
Scripts to retrieve human genome annotation files and process them into desirable formats

Created on 4/27/2014
"""

__author__ = """\n""".join(['Zichen Wang <zichen.wang@mssm.edu>'])

import os
import urllib2, gzip

def retrieve_annotation(url=None):
	os.chdir('./annotations')
	if url is None:
		## url for UCSC knownGene tabel
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

os.chdir('./annotations')

with gzip.open('knownGenes.txt.gz', 'rb') as f:
	for line in f:
		sl = line.strip().split('\t')
		if line.startswith('#'):
			header = sl
			for i in header:
				print i
		else:
			break


