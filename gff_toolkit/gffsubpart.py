__author__ = 'rensholmer'

import re
import pprint

class TranslateError(Exception):
	pass
class SeqError(Exception):
	pass
class CoordinateError(Exception):
	pass
class IDError(Exception):
	pass

class GffSubPart(object):
	"""
	Work in progress: contained by Gff object
	basically this is one line of a gff file, with parents/children indicated by ID
	some methods to determine overlap between other GffSubPart instances
	TODO:
	work on getter/setter for seq
	work on getter/setter for ID, ratt uses locus_tag for ID of transferred annotations, currently this is handled by parser.py
	Make this a sort of 'baseclass', for instance only mRNA should have a pep option. Do this by subclassing
	"""
	_featuretypes = ('gene','mRNA','CDS','exon','three_prime_UTR','five_prime_UTR','match','match_part','protein_match','transcript_match','biological_region','polypeptide')
	_filetypes = ('gff','tbl')
	def __init__(self,*args,**kwargs):
		"""
		Letsa go!
		"""
		self.seqid = None
		self.source = None
		self.featuretype = None
		self.start = None
		self.end = None
		self.score = None
		self.strand = None
		self.phase = None
		#self.parents = [] # list of GffSubPart IDs
		self.children = [] #list of GffSubPart IDs
		self.comments = None #Not implemented
		self._pep = ''
		self._seq = ''
		self._key = None #Unique ID
		self.container = None # This will be a Gff object so that the GffSubFeature can find its child/parent objects
		self.attributes = {}
		if len(args) >= 8:
			self.seqid = args[0]
			self.source = args[1]
			self.featuretype = args[2]
			self.start = int(args[3])
			self.end = int(args[4])
			if not self.end >= self.start:
				e = 'End must be greater than or equal to Start: {0}\t{1}\t{2}'.format(self.seqid,self.start,self.end)
				raise CoordinateError(e)
			if args[5] == '.':
				self.score = args[5]
			else:
				self.score = float(args[5])
			self.strand = args[6]
			if args[7] == '.':
				self.phase = args[7]
			else:
				self.phase = int(args[7])
		elif len(args) > 0:
			e = 'Line is not proper gff: ' + '\t'.join(args)
			raise AttributeError(e)
		if len(args) == 9:
			#split column 9 on ';' to get attributes
			#split attributes on '=' to get key:value
			self.attributes = dict(a.split('=',1) for a in args[8].split(';'))
			#split value on ',' to get values
			for attribute,value in self.attributes.iteritems():
				self.attributes[attribute] = value.split(',')
		if not self.ID:
			if kwargs.get('strict',False):
				e = 'No ID found in attributes: ' + '\t'.join(args)
				raise AttributeError(e)
			elif kwargs.get('filetype','standard') == 'manual':
				newID = self.attributes.get('Name',False)
				if newID:
					self.ID = newID[0]
				else:
					e = 'Cannot process manual annotation without Name and ID: {0}'.format(self.stringify().strip())
					raise IDError(e)
			else:
				self.ID = '{0}.{1}'.format(self.parents[0],self.featuretype)
	def __str__(self):
		"""
		Lets print this
		"""
		while True:
			try:
				seq = self.seq
				pep = self.pep
			except:
				seq = ''
				pep = ''
			s = {	'ID':self.ID,'seqid':self.seqid,'source':self.source,'type':self.featuretype,'start':self.start,'end':self.end,
					'score':self.score,'strand':self.strand,'phase':self.phase,'attributes':self.attributes,'children':self.children,'parents':self.parents,'seq':seq,'pep':pep }

			return pprint.pformat(s,indent=1)
		#code below can be used to print gff style, currently done by stringify() method
		s = [self.seqid,self.source,self.featuretype,self.start,self.end,self.score,self.strand,self.phase]
		s = [str(x) for x in s]
		s = '\t'.join(s)
		a = []
		for key in self.attributes:
			a.append('{0}={1}'.format(key,','.join(self.attributes[key])))
		a = ';'.join(a)
		return '{0}\t{1}\n'.format(s,a)
		#return '{0}\t{1}\t{2}'.format(s,a,self.seq)
	def __repr__(self):
		"""
		"""
		#return self.__str__()
		return '{0}({1})'.format(self.__class__,self.__dict__)
	def __eq__(self, other):
		"""
		Compare two GffSubFeatures on content or ID
		"""
		if isinstance(other,GffSubPart):
			return self.__dict__ == other.__dict__
		elif isinstance(other,basestring):
			return self.ID == other
		else:
			e = 'Cannot compare GffSubPart to object of type {0}'.format(type(other))
			raise NotImplementedError(e)
	def todict(self):
		"""
		Return a dict
		"""
		try:
			seq = self.seq
			pep = self.pep
		except:
			seq = ''
			pep = ''
		dic = {	'ID':self.ID,
				'seqid':self.seqid,
				'source':self.source,
				'type':self.featuretype,
				'start':self.start,
				'end':self.end,
				'score':self.score,
				'strand':self.strand,
				'phase':self.phase,
				'attributes':self.attributes,
				'children':self.children,
				'parents':self.parents,
				'seq':seq,
				'pep':pep,
				'file':self.container.filename}
		return dic
	def fromdict(self,dic):
		pass
	@property
	def _start(self):
		"""
		This is a hack for the intervaltree: it prevents NULL intervals when start == end
		"""
		return self.start - 1

	@property
	def target(self):
		_target = self.attributes.get('Target',[False])[0]
		if not _target:
			return
		_target = _target.split()
		target_dict = {'ID':_target[0],'start':_target[1],'end':_target[2]}
		return target_dict
	@target.setter
	def target(self,value):
		if not isinstance(value,dict):
			raise ValueError()
		if not 'ID' in value.keys() or 'start' in value.keys() or 'end' in value.keys():
			raise ValueError()
		_target = '{0} {1} {2}'.format(value['ID'],value['start'],value['end'])
		self.attributes['Target'] = [_target]
	@target.deleter
	def target(self):
		self.attributes.pop('Target',False)
	
	@property
	def ID(self):
		return self.attributes.get('ID',[None])[0]
	@ID.setter
	def ID(self,value):
		self.attributes['ID'] = [value]

	@property
	def forward(self):
		return self.strand == '+'

	@property
	def reverse(self):
		return self.strand == '-'

	@property
	def parents(self):
		return self.attributes.get('Parent',[])
	@parents.setter
	def parents(self,value):
		if not isinstance(value,(tuple,list)):
			e = '{0} is not a valid type for parents property'.format(type(value))
			raise TypeError(e)
		self.attributes['Parent'] = value

	@property
	def siblings(self):
		"""
		return sorted list of siblings
		"""
		sibs = []
		for parent_ID in self.parents:
			for parent in self.container[parent_ID]:
				for sib in self.container.get_children(parent,featuretype=self.featuretype):
					sibs.append(sib)
		return sorted(sibs,key = lambda x: x.get_start(),reverse = self.reverse)
	
	@property
	def seq(self):
		'''
		if cds.start + 1 == cds.end and index == len(children) - 1:
			cds.seq = self.seq[cds.seqid][cds.start-1]
		else:
			cds.seq = self.seq[cds.seqid][cds.start-1:cds.end]
		'''
		if self.featuretype == 'CDS':
			#if self.start + 1 == self.end and self.siblings.index(self) == len(self.siblings) - 1:
			#	seq = self.container.seq[self.seqid][self.start-1]
			#else:
			seq = self.container.seq[self.seqid][self.start-1:self.end]
			if self.reverse:
				seq = self._revcomp(seq)
			#print '+',seq
			return seq
		elif self.featuretype == 'mRNA':
			seq = ''
			cds_list = sorted(self.container.get_children(self,featuretype='CDS'),key = lambda x: x.get_start(),reverse=self.reverse)
			for i,cds in enumerate(cds_list):
				if i == 0 and cds.phase != '.':
					#if cds.phase != '.':
					#	pass
					seq += cds.seq#[cds.phase:]
				else:
					seq += cds.seq
			return seq
		else:
			raise NotImplementedError()

	@property
	def pep(self):
		if self.featuretype == 'mRNA':
			cdss = sorted(self.get_children(self),key = lambda x: x.get_start(),reverse=self.reverse)
			cds = cdss[0]
			phase = cds.phase
			return self._translate(self.seq[phase:])
		elif self.featuretype == 'CDS':
			phase = self.phase
			if phase =='.':
				phase = 0
			return self._translate(self.seq[phase:])
		else:
			raise NotImplementedError()

	@property
	def gff_fields(self):
		fields = (self.seqid,self.source,self.featuretype,self.start,self.end,self.score,self.strand,self.phase)
		fields = [str(f) for f in fields]
		attributes = ';'.join(['{0}={1}'.format(key,','.join(value)) for key,value in self.attributes.iteritems() if value])
		#for key,value in self.attributes.iteritems():
		#	attributes += '{0}={1}'.format(key,','.join(value)) + ';'
		return fields + [attributes]

	@property
	def gtf_fields(self):
		fields = (self.seqid,self.source,self.featuretype,self.start,self.end,self.score,self.strand,self.phase)
		fields = [str(f) for f in fields]
		attributes = '; '.join(['{0} "{1}"'.format(key,','.join(value)) for key,value in self.attributes.iteritems() if value])
		return fields + [attributes]
	

	def get_children(self,*args,**kwargs):
		for sub in self.container.get_children(self,*args,**kwargs):
			yield sub
	
	def stringify(self,filetype='gff'):
		"""
		Return gff style tab separated string 
		"""
		if filetype not in self._filetypes:
			e = '{0} is not a valid filetype'.format(filetype)
			raise TypeError(e)
		elif filetype == 'gff':
			s = (self.seqid,self.source,self.featuretype,self.start,self.end,self.score,self.strand,self.phase)
			s = (str(x) for x in s)
			s = '\t'.join(s)
			a = []
			for key,value in self.attributes.iteritems():
				a.append('{0}={1}'.format(key,','.join(value)))
			a = ';'.join(a)
			return '{0}\t{1}\n'.format(s,a)
		elif filetype == 'tbl':
			lines = ['{0}\t{1}\t{2}'.format(self.start,self.end,self.featuretype)]
			if self.featuretype == 'gene':
				name = self.attributes.get('Name',False)[0]
				if name:
					lines.append('\t\t\tgene\t{0}'.format(name))
				lines.append('\t\t\tlocus_tag\t{0}'.format(self.ID))
			elif self.featuretype in ['mRNA','CDS']:
				lines.append('\t\t\tproduct\tNone')
			return '\n'.join(lines)
	def set_attribute(self,key,value):
		if isinstance(value,basestring):
			self.attributes.setdefault(key,[]).append(value)
		elif hasattr(value,'__iter__'):
			for v in value:
				self.attributes.setdefault(key,[]).append(value)
		else:
			e = '{0} is not a valid type for GffSubPart attribute values'.format(type(value))
			raise TypeError(e)

	def getattribute(self,key):
		for value in self.attributes.get(key,[]):
			yield value

	def set_start(self,value):
		"""
		Used in Gff._change_cds to change start/stop based on strand
		:param value: int
		:return:
		"""
		if self.strand == '+':
			self.start = value
		else:
			self.end = value
	def get_start(self):
		"""
		:return: start if strand == + else end
		"""
		return self.start if self.strand == '+' else self.end

	def get_end(self):
		"""
		:return: end if strand == + else start
		"""
		return self.end if self.strand == '+' else self.start

	def set_end(self,value):
		"""
		Used in Gff._change_cds to change start/stop based on strand
		:param value:  int
		:return:
		"""
		if self.strand == '+':
			self.end = value
		else:
			self.start = value


	def _find_orf(self,fasta_dic=None):
		orfs = []
		if not self.seq:
			raise Exception('Can not find ORF if self.seq not set')
		frame1 = (0,self.seq[0:])
		frame2 = (1,self.seq[1:-2])
		frame3 = (2,self.seq[2:-1])
		for offset,frame in frame1,frame2,frame3:
			codons = re.findall('...',frame)
			starts = [index * 3 for index,codon in enumerate(codons) if codon == 'ATG']
			if not starts:
				starts = [index * 3 for index,codon in enumerate(codons) if codon == 'CTG']
				if not starts:
					continue
			stops = [(index+1) * 3 for index,codon in enumerate(codons) if codon in ['TAA','TGA','TAG']]
			if not stops:
				continue
			
			start = min(starts) + offset
			stop = min(stops) + offset
			if stop - start > 0: # test for positive mRNA length
				orfs.append((start,stop))
		if not orfs:
			return False
		longest_orf = max(orfs,key = lambda x: x[1]-x[0])
		print 'SEQ',self.seq
		print 'ORF',self.seq[longest_orf[0]:longest_orf[1]]
		return longest_orf

	def getnested(self,reverse=False,featuretype=None):
		"""
		generator to get all nested subfeatures.
		Reverse == True starts at the tips of the tree: CDS/exon --> mRNA --> gene
		Reverse == False starts at the top of the tree: gene --> mRNA --> CDS/exon
		"""
		if not reverse and (featuretype == None or featuretype == self.featuretype):
			yield self
		if self.children:
			for child in self.children:
				for c in child.getnested(reverse=reverse,featuretype=featuretype):
					#if featuretype == None or featuretype == self.featuretype:
					yield c
		if reverse and (featuretype == None or featuretype == self.featuretype):
			yield self
	@staticmethod
	def _revcomp(seq):
		comp = {'A':'T','T':'A','G':'C','C':'G','N':'N',
				'a':'t','t':'a','g':'c','c':'g','n':'n'}
		return  ''.join([comp[n] for n in seq[::-1]])
	@staticmethod
	def _translate(seq):
		p = {'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
		     'AGG': 'R', 'AGC': 'S', 'GTA': 'V',
		     'AGA': 'R', 'ACT': 'T', 'GTG': 'V',
		     'AGT': 'S', 'CCA': 'P', 'CCC': 'P',
		     'GGT': 'G', 'CGA': 'R', 'CGC': 'R',
		     'TAT': 'Y', 'CGG': 'R', 'CCT': 'P',
		     'GGG': 'G', 'GGA': 'G', 'GGC': 'G',
		     'TAA': '*', 'TAC': 'Y', 'CGT': 'R',
		     'TAG': '*', 'ATA': 'I', 'CTT': 'L',
		     'ATG': 'M', 'CTG': 'L', 'ATT': 'I',
		     'CTA': 'L', 'TTT': 'F', 'GAA': 'E',
		     'TTG': 'L', 'TTA': 'L', 'TTC': 'F',
		     'GTC': 'V', 'AAG': 'K', 'AAA': 'K',
		     'AAC': 'N', 'ATC': 'I', 'CAT': 'H',
		     'AAT': 'N', 'GTT': 'V', 'CAC': 'H',
		     'CAA': 'Q', 'CAG': 'Q', 'CCG': 'P',
		     'TCT': 'S', 'TGC': 'C', 'TGA': '*',
		     'TGG': 'W', 'TCG': 'S', 'TCC': 'S',
		     'TCA': 'S', 'GAG': 'E', 'GAC': 'D',
		     'TGT': 'C', 'GCA': 'A', 'GCC': 'A',
		     'GCG': 'A', 'GCT': 'A', 'CTC': 'L',
		     'GAT': 'D'}
		codons = re.findall('...',seq)
		pep = [p[codon.upper()] if 'N' not in codon.upper() else 'X' for codon in codons]
		return ''.join(pep)

	def match(self,other):
		if not isinstance(other,GffSubPart):
			e = '{0} is not of type GffSubPart, can not match'
			raise TypeError(e)
		if not self.strand == other.strand:
			return False
		if self.start >= other.start and self.start <= other.end:
			return True
		elif other.start >= self.start and other.start <= self.end:
			return True
		else:
			return False
	def exact_match(self,other):
		if not isinstance(other,GffSubPart):
			e = '{0} is not of type GffSubPart, can not match'
			raise TypeError(e)
		if self.strand == other.strand and self.start == other.start and self.end == other.end:
			return True
		else:
			return False

	def _match(self,other):
		if not isinstance(other,GffSubPart):
			e = '{0} is not of type GffSubPart, can not match'
			raise TypeError(e)
		if self.strand == other.strand:
			if self.start == other.start and self.end == other.end:
				return 2
			elif self.start >= other.start and self.start <= other.end:
				return 1
			elif other.start >= self.start and other.start <= self.end:
				return 1
			else:
				return 0
		else:
			return 0


