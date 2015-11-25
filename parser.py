__author__ = 'rensholmer'

from .gff import Gff
from .gffsubpart import GffSubPart,TranslateError
from itertools import groupby

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
	#@profile
	def _ratt(self):
		print '_ratt()'
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
				featuremap.setdefault(newID,set()).add(feature.ID)
			else:
				raise ValueError('old ID already exists!')

			for subfeature in self.gff.get_children(feature,featuretype='CDS'):
				subfeature.parents = [newID]
				oldsubID = subfeature.ID
				newsubID = '{0}.{1}'.format(newID,subfeature.featuretype)
				if newsubID not in childmap[newID]:
					childmap[newID].append(newsubID)
				featuremap.setdefault(newsubID,set()).add(oldsubID)

				subfeature.ID = newsubID
				subfeature.attributes['ID'] = [newsubID]
				subfeature.attributes['Parent'] = [newID]
				subfeature.source = 'ratt'
			feature.ID = newID
			feature.source = 'ratt'
			feature.attributes['ID'] = [newID]
			try:
				del feature.attributes['featflags']
			except:
				pass
			del feature.attributes['locus_tag']
			del feature.attributes['gene']
			del feature.attributes['transl_table']
				
		for x in featuremap:
			for y in featuremap[x]:
				self.gff.name_index[x] = self.gff.name_index[y]
			for y in featuremap[x]:
				self.gff.name_index.pop(y,set())
		for c in childmap:
			for t in self.gff[c]:
				t.children = childmap[c]
	def _manual(self):
		pass
	def _remove_noncoding(self):
		print '_remove_noncoding()'
		if not self.fasta_file:
			e = 'Can not remove non coding annotations without sequence data'
			raise NotImplementedError(e)
		remove = []
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
		print 'Remove {0} genes because they do not encode a protein'.format(len(remove))
		self.gff.remove(remove)
	def _remove_partial(self):
		print '_remove_partial()'
		remove = []
		for feature in self.gff.getitems(featuretype='mRNA'):
			if feature.attributes['partial'][0] == 'True':
				remove.append(feature.ID)
		print 'Remove {0} genes because they are a partial transfer'.format(len(remove))
		self.gff.remove(remove)
	def parse(self):
		for subpart in self._readlines():
			self.gff.update(subpart)
		self.gff.set_children()
		if self.filetype == 'ratt':
			self._ratt()
			#self._remove_partial()
		if self.fasta_file:
			self.gff.add_fasta(self.fasta_file)
			self.gff.getseq(topfeaturetype='mRNA',subfeaturetype='CDS')
		if self.remove_noncoding:
			self._remove_noncoding()			
		self.gff.make_index()
		self.parsed = True
		return self.gff

def fasta_iter(fasta_file):
	with open(fasta_file,'rU') as filehandle:
		faiter = (x[1] for x in groupby(filehandle, lambda line: line[0] == '>'))
		for header in faiter:
			header = header.next()[1:].strip()
			seq = ''.join(s.strip() for s in faiter.next())
			yield header,seq

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

