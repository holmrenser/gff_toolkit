__author__ = 'rensholmer'

from .gff import Gff
from .gffsubpart import GffSubPart

class Parser(object):
	"""
	"""
	_filetypes = ['standard','ratt','manual']
	def __init__(self,gff_file,filetype='standard'):
		if filetype not in self._filetypes:
			raise TypeError()
		self.filename = gff_file
		self.filetype = filetype
		self.gff = Gff()
		self.parse()
	def _readlines(self):
		with open(self.filename,'rU') as fh:
			for line in fh:
				if line[0] == '#' or not line.strip():
					continue
				parts = line.strip().split('\t')
				yield GffSubPart(*parts)
	def parse(self):
		for subpart in self._readlines():
			self.gff.update(subpart)
		return self.gff

def linereader(gff_file,strict=True):
	with open(gff_file,'rU') as fh:
		for line in fh:
			if line[0] == '#' or not line.strip():
				continue
			parts = line.strip().split('\t')
			yield GffSubPart(*parts,strict=strict)

#@profile
def parser(gff_file,strict=True,ratt=False,fasta_file=None):
	gff = Gff(filename=gff_file)
	for subfeature in linereader(gff_file,strict=strict):
		gff.update(subfeature)
	gff.set_children()
	if ratt:
		remove = []
		featuremap = {}
		childmap = {}
		for feature in gff.getitems(level='mRNA'):
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

			for subfeature in gff.get_children(feature,featuretype='CDS'):
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
			gff.features[x[0]][featuremap[x]] = gff.features[x[0]].pop(x[1])
		for c in childmap:
			for t in gff[c]:
				t.children = childmap[c]

	if fasta_file != None:
		gff.add_fasta(fasta_file)
		gff.getseq(toplevel='mRNA',sublevel='CDS')
		if ratt:
			for feature in gff.getitems(level='mRNA'):
				try:
					pep = feature.pep
				except:
					remove.append(feature.ID)
					continue
				if pep[0] != 'M':
					remove.append(feature.ID)
				elif pep[-1] != '*':
					remove.append(feature.ID)
				elif '*' in pep[1:-1]:
					remove.append(feature.ID)
			gff.remove(remove)
	gff.make_index()
	return gff
