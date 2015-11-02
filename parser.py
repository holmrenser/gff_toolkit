__author__ = 'rensholmer'

from .gff import Gff
from .gffsubpart import GffSubPart

class Parser(object):
	"""
	Parser object gff3 formatted annotation files.
	The parse() method return the processed Gff object.
	Has several methods for known incorrectly formatted files: _ratt() and _manual()
	"""
	_filetypes = ['standard','ratt','manual']
	def __init__(self,gff_file,filetype='standard',fasta_file=None,remove_noncoding=False):
		if filetype not in self._filetypes:
			raise TypeError()
		self.filename = gff_file
		self.filetype = filetype
		self.fasta_file = fasta_file
		self.remove_noncoding = remove_noncoding
		self.gff = Gff(filename=gff_file)
		self.parsed = False
	def _readlines(self):
		with open(self.filename,'rU') as fh:
			for line in fh:
				if line[0] == '#' or not line.strip():
					continue
				parts = line.strip().split('\t')
				yield GffSubPart(*parts)
	def _ratt(self):
		remove = []
		featuremap = {}
		childmap = {}
		for feature in self.gff.getitems(featuretype='mRNA'):
			oldID = feature.ID
			newID = feature.attributes['locus_tag'][0]
			if oldID == newID:
				continue
			if newID not in childmap:
				childmap[newID] = []

			if oldID not in featuremap:
				featuremap[(feature.seqid,oldID)] = newID
			else:
				raise ValueError('old ID already exists!')

			for subfeature in self.gff.get_children(feature,featuretype='CDS'):
				subfeature.parents = [newID]
				oldsubID = subfeature.ID
				newsubID = '{0}.{1}'.format(newID,subfeature.featuretype)
				if newsubID not in childmap[newID]:
					childmap[newID].append(newsubID)
				if oldsubID not in featuremap:
					featuremap[(subfeature.seqid,oldsubID)] = newsubID
				elif newsubID != featuremap[oldsubID]:
					raise ValueError('Old sub ID has multiple new IDs!')
				subfeature.ID = newsubID
				subfeature.attributes['ID'] = [newsubID]
				subfeature.attributes['Parent'] = [newID]
			feature.ID = newID
			feature.attributes['ID'] = [newID]
		for x in featuremap:
			self.gff.features[x[0]][featuremap[x]] = self.gff.features[x[0]].pop(x[1])
		for c in childmap:
			for t in self.gff[c]:
				t.children = childmap[c]
	def _manual(self):
		pass

	def parse(self):
		for subpart in self._readlines():
			self.gff.update(subpart)
		self.gff.set_children()
		if self.filetype == 'ratt':
			self._ratt()
		if self.fasta_file:
			self.gff.add_fasta(self.fasta_file)
			self.gff.getseq(topfeaturetype='mRNA',subfeaturetype='CDS')
		if self.remove_noncoding:
			if not self.fasta_file:
				e = 'Can not remove non coding annotations without sequence data'
				raise NotImplementedError(e)
			for feature in self.gff.getitems(featuretype='mRNA'):
				try:
					pep = feature.pep
				except TranslateError:
					remove.append(feature.ID)
					continue
				if pep[0] != 'M':
					remove.append(feature.ID)
				elif pep[-1] != '*':
					remove.append(feature.ID)
				elif '*' in pep[1:-1]:
					remove.append(feature.ID)
			self.gff.remove(remove)
		self.gff.make_index()
		self.parsed = True
		return self.gff

def linereader(gff_file,strict=True):
	"""
	Reads a single gff formatted line
	"""
	with open(gff_file,'rU') as fh:
		for line in fh:
			if line[0] == '#' or not line.strip():
				continue
			parts = line.strip().split('\t')
			yield GffSubPart(*parts,strict=strict)

def parser(gff_file,filetype='standard',fasta_file=None,remove_noncoding=False):
	"""
	Wrapper function for Parser object.
	"""
	p = Parser(gff_file,filetype=filetype,fasta_file=fasta_file,remove_noncoding=remove_noncoding)
	return p.parse()

