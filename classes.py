"""
All classes used.

Created on 4/27/2014
"""

def _check_coordinates(coordinates):
	if len(coordinates) == 3 and type(coordinates[0])== str \
		and type(coordinates[1]) == int and type(coordinates[2]) == int:
		return coordinates
	else:
		raise ValueError ('The coordinates are not in a correct format.')

class GeneElement(object):
	"""docstring for GeneElement"""
	def __init__(self, chrom, start, end):
		self.chrom = chrom
		self.position = [chrom, int(start), int(end)]		

	def __len__(self):
		size = abs(self.position[1] - self.position[2])
		return size

	def __contains__(self, coordinates):
		'''
		The format of coordinates should be:
		coordinates = [chrom, start(int), end(int)]
		'''
		coordinates = _check_coordinates(coordinates)
		start = min(self.position[1], self.position[2])
		end = max(self.position[1], self.position[2])
		if coordinates[0] == self.chrom:
			if start < coordinates[1] < end or start < coordinates[2] < end:
				return True
			else: 
				return False
		else:
			return False

	def intersect(self, coordinates):
		"""return the number of bases in the intersection between
		GeneElement and a given coordinates"""
		if coordinates not in self:
			return 0
		else:
			start = min(self.position[1], self.position[2])
			end = max(self.position[1], self.position[2])
			coor_bases = set(range(min(coordinates[1],coordinates[2]),max(coordinates[1],coordinates[2])))
			elem_bases = set(range(start,end))
			overlap_size = len(coor_bases & elem_bases)
			return overlap_size

	def percentage_overlapped(self, coordinates):
		"""calculate the percentage of the overlapping bases over 
		the bases/length of GeneElement"""
		if coordinates not in self:
			return 0.
		else:
			overlap_size = intersect(self, coordinates)
			return float(overlap_size)/len(self)

class Gene(GeneElement):
	"""docstring for Gene"""
	def __init__(self, line): # line is the line in UCSC gene table
		sl = line.strip().split('\t') # sl is the splited line in gene table
		self.chrom = sl[2]
		GeneElement.__init__(self, sl[2], int(float(sl[4])), int(float(sl[5]))) ## chrom, strand, txStart, txEnd 
		self.refSeqId = sl[1]
		self.CDS = GeneElement(sl[2], sl[6], sl[7])
		self.exonCount = int(sl[8])
		self.exons = []
		for s, e in zip(sl[9].strip(',').split(','), sl[10].strip(',').split(',')):
			Exon = GeneElement(self.chrom, s, e)
			self.exons.append(Exon)
		self.name = sl[12]

	def __str__(self):
		return self.refSeqId

class Indel(GeneElement):
	"""docstring for Indel"""
	def __init__(self, coordinates, confidence, typeStr):
		GeneElement.__init__(self, *coordinates)
		self.confidence = confidence
		if typeStr.startswith('insertion'):
			self.typeStr = 'insertion'
		elif typeStr.startswith('deletion'):
			self.typeStr = 'deletion'
		else:
			raise ValueError('Unrecognizable typeStr %s'%typeStr)

		