"""
All classes used.

Created on 4/27/2014
"""


class GeneElement(object):
	"""docstring for GeneElement"""
	def __init__(self, chrom, start, end):
		self.chrom = chrom
		self.position = [chrom, int(start), int(end)]		

	def __len__(self):
		size = abs(self.position[1] - self.position[2])
		return size

	def __contains__(self, coordinates):
		
		if coordinates[0] == self.chrom and coordinates[1] < max(self.position[1], self.position[2]) \
			and coordinates[1] > min(self.position[1], self.position[2]):
			return True
		else:
			return False

class Gene(GeneElement):
	"""docstring for Gene"""
	def __init__(self, sl): # sl is the splited line in gene table
		self.refSeqId = sl[1]
		self.chrom = sl[2]
		self.position = [self.chrom, int(sl[4]), int(sl[5])] ## chrom, strand, txStart, txEnd 
		self.CDS = GeneElement(sl[2], sl[6], sl[7])
		self.exons = []
		for s, e in zip(sl[9].split(','), sl[10].split(',')):
			Exon = GeneElement(self.chrom, s, e)
			self.exons.append(Exon)
		self.name = sl[12]

	def __str__(self):
		return self.refSeqId


