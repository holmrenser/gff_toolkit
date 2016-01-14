__author__ = 'rensholmer'

from gff import Gff
from gffsubpart import GffSubPart,TranslateError
from itertools import groupby
import uuid

class Parser(object):
	"""
	Parser object gff3 formatted annotation files.
	The parse() method return the processed Gff object.
	Has several methods for known incorrectly formatted files: _ratt() and _manual(), which are called by specifying filetype
	"""
	_filetypes = ['standard','ratt','manual']
	def __init__(self,gff_file,filetype='standard',fasta_file=None,remove_noncoding=False,limit=None,author=None):
		if filetype not in self._filetypes:
			raise TypeError()
		self.filename = gff_file
		if filetype not in self._filetypes:
			e = '{0} is not a valid filetype'.format(filetype)
			raise TypeError(e)
		self.filetype = filetype
		self.fasta_file = fasta_file
		self.remove_noncoding = remove_noncoding
		if limit:
			if not isinstance(limit,dict):
				e = '{0} is not a valid type for limit'.format(type(limit))
				raise TypeError(e)
			else:
				#reformat limit values into lists if they are not already lists
				self.limit = {}
				for key,value in limit.iteritems():
					if isinstance(value,(list,tuple)):
						self.limit[key] = value
					elif isinstance(value,basestring):
						self.limit[key] = [value]
				#self.limit = {key:[value] for key,value in limit.iteritems() if isinstance(value,basestring) or if isinstance(value,(list,tuple))}
		else:
			self.limit = limit
 		self.gff = Gff(filename=gff_file)
 		self.author = author
		#self.parsed = False
	def _readlines(self):
		"""
		Generator that yields GffSubPart objects formatted according to filetype. 
		Limit is a dictionary that can limit the lines parsed based on some attributes like 
		"""
		with open(self.filename,'rU') as fh:
			for line in fh:
				if line[0] == '#' or not line.strip():
					continue
				parts = line.strip().split('\t')
				#print parts
				if self.limit != None:
					#print self.limit
					#Series of list comprehensions to check if any of the relevant attributes of sub are in the limit dictionary
					#An empty list is falsy, so if the relevant value is not in the limit dictionary, continue to the next line
					if not [x for x in self.limit.get('seqid',[parts[0]]) if x == parts[0]]:
						continue
					elif not [x for x in self.limit.get('source',[parts[1]]) if x == parts[1]]:
						continue
					elif not [x for x in self.limit.get('featuretype',[parts[2]]) if x == parts[2]]:
						continue
					elif not [x for x in self.limit.get('strand',[parts[6]]) if x == parts[6]]:
						continue
				sub = GffSubPart(*parts,filetype=self.filetype)
				#print sub.stringify().strip()
				yield sub
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
		ann = {}
		for sub in self.gff.getitems(featuretype='CDS'):
			if sub.attributes.get('Name',False):
				ann.setdefault(sub.ID,[]).append(sub)
		for gene,cds_list in ann.iteritems():
			self._assert_equal(cds_list)
			template = cds_list[0]
			exon_counter = 0
			#GffSubPart elements
			seqid = template.seqid
			source = template.source
			start = min(s.start for s in cds_list)
			end = max(s.end for s in cds_list)
			score = template.score
			strand = template.strand
			phase = template.phase
			#GffSubPart attributes
			gene_ID = uuid.uuid4()
			mrna_ID = '{0}.1'.format(gene_ID)
			name = template.attributes.get('Name',False)[0]
			author = self.author#template.attributes.get('Created by',False)[0]
			gene_attributes = 'ID={0};Name={1};Created by={2}'.format(gene_ID,name,author)
			mrna_attributes = 'ID={0};Parent={1};Name={2};Created by={3}'.format(mrna_ID,gene_ID,name,author)
			#make the gene GffSubPart
			gene = GffSubPart(seqid,source,'gene',start,end,score,strand,phase,gene_attributes)
			#make the mRNA GffSubPart
			mrna = GffSubPart(seqid,source,'mRNA',start,end,score,strand,phase,mrna_attributes)
			self.gff.update(gene)
			self.gff.update(mrna)
			if gene.strand == '+':
					reverse = True
			else:
					reverse = False
			for cds in sorted(cds_list,key=lambda x:x.get_start(),reverse=reverse):
				exon_counter += 1
				old_ID = cds.ID
				cds.ID = '{0}.CDS{1}'.format(mrna_ID,exon_counter)
				cds.attributes['ID'] = ['{0}.CDS{1}'.format(mrna_ID,exon_counter)]
				#print self.gff.name_index[old_ID]
				self.gff.name_index.setdefault(cds.ID,set()).add(cds._key)
				self.gff.name_index[old_ID] = {x for x in self.gff.name_index[old_ID] if x != cds._key}
				if len(self.gff.name_index[old_ID]) == 0:
					self.gff.name_index.pop(old_ID)
				cds.parents = [mrna_ID]
				cds.attributes['Parent'] = [mrna_ID]
				cds.attributes.pop('modified by',[])
				cds.attributes.pop('created by',[])
				cds.attributes['Created by'] = [self.author]
				#gff.update(sub)
		#self.gff.set_children()

	@staticmethod
	def _assert_equal(subslist):
		first = subslist[0]
		assert all(x.seqid == first.seqid for x in subslist)
		assert all(x.source == first.source for x in subslist)
		assert all(x.strand == first.strand for x in subslist)
		try:
			assert all(x.attributes['Name'] == first.attributes['Name'] for x in subslist),[x.attributes['Name'] for x in subslist]
		except KeyError as e:
			print first
			raise e
		except AssertionError as e:
			raise e
		try:
			pass
			#assert all(x.attributes['Created by'] == first.attributes['Created by'] for x in subslist)
		except KeyError as e:
			print first
			raise e
		
	def _remove_noncoding(self):
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
				if feature.seq[0:3] != 'CTG':
					remove.append(feature.ID)
			elif pep[-1] != '*':
				remove.append(feature.ID)
			elif '*' in pep[1:-1]:
				remove.append(feature.ID)
		print 'Remove {0} genes because they do not encode a protein'.format(len(remove))
		self.gff.remove(remove)

	def _remove_partial(self):
		#print '_remove_partial()'
		remove = []
		for feature in self.gff.getitems(featuretype='mRNA'):
			if feature.attributes['partial'][0] == 'True':
				remove.append(feature.ID)
		#print 'Remove {0} genes because they are a partial transfer'.format(len(remove))
		self.gff.remove(remove)

	def parse(self):
		for subpart in self._readlines():
			self.gff.update(subpart)
		
		if self.filetype == 'ratt':
			self.gff.set_children()
			self._ratt()
		elif self.filetype == 'manual':
			self._manual()
			self.gff.set_children()
		else:
			self.gff.set_children()
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

def parser(gff_file,filetype='standard',fasta_file=None,remove_noncoding=False,limit=None,author=None):
	"""
	Wrapper function for Parser object.
	"""
	p = Parser(gff_file,filetype=filetype,fasta_file=fasta_file,remove_noncoding=remove_noncoding,limit=limit,author=author)
	return p.parse()

