__author__ = 'rensholmer'

import re
import pprint

class TranslateError(Exception):
	pass
class SeqError(Exception):
	pass
class CoordinateError(Exception):
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
	_featuretypes = ['gene','mRNA','CDS','exon','three_prime_UTR','five_prime_UTR','match','match_part','protein_match','transcript_match','biological_region']
	_filetypes = ['gff','tbl']
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
		self.parents = [] # list of GffSubPart IDs
		self.children = [] #list of GffSubPart IDs
		self.comments = None #Not implemented
		self._pep = ''
		self._seq = ''
		self._key = None #Unique ID
		self.container = None # This will be a Gff object so that the GffSubFeature can find its child/parent objects
		if len(args) == 9:
			self.seqid = args[0]
			self.source = args[1]
			self.featuretype = args[2]
			self.start = min(int(args[3]),int(args[4]))#int(args[3])
			self.end = max(int(args[3]),int(args[4]))#int(args[4])
			if not self.end > self.start:
				e = 'End must be larger than Start: {0}\t{1}\t{2}'.format(self.seqid,self.start,self.end)
				#raise CoordinateError(e)#,[self.seqid,self.start,self.end]
				self.end += 1
			if args[5] == '.':
				self.score = args[5]
			else:
				self.score = int(args[5])
			self.strand = args[6]
			if args[7] == '.':
				self.phase = args[7]
			else:
				self.phase = int(args[7])
			att = args[8].split(';')
			self.attributes = { a.split('=')[0]: a.split('=')[1].split(',') for a in att}
		elif len(args) > 0:
			e = 'Line is not proper gff: ' + '\t'.join(args)
			raise AttributeError(e)
		if 'ID' in self.attributes:
			self.ID = self.attributes['ID'][0]
		self.ID = self.attributes.get('ID',[None])[0]
		if not self.ID:
			if kwargs.get('strict',False):
				e = 'No ID found in attributes: ' + '\t'.join(args)
				raise AttributeError(e)
			else:
				self.ID = '{0}.{1}'.format(self.attributes['Parent'][0],self.featuretype)
		for parent in self.attributes.get('Parent',[]):
			self.parents.append(parent)
		#self._key = (self.ID,self.seqid,self.start,self.end,self.strand)
	def __str__(self):
		"""
		Lets print this
		"""
		while True:
			s = {	'ID':self.ID,'seqid':self.seqid,'source':self.source,'type':self.featuretype,'start':self.start,'end':self.end,
					'score':self.score,'strand':self.strand,'phase':self.phase,'attributes':self.attributes,'children':self.children,'parents':self.parents,'seq':self.seq }

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
	def stringify(self,filetype='gff'):
		"""
		Return gff style tab separated string 
		"""
		if filetype not in self._filetypes:
			e = '{0} is not a valid filetype'.format(filetype)
			raise TypeError(e)
		if filetype == 'gff':
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
	@property
	def seq(self):
		return self._seq
	@seq.setter
	def seq(self,value):
		self._seq = value
	@seq.deleter
	def seq(self):
		self._seq = ''

	@property
	def pep(self):
		if self._pep:
			return self._pep
		if self._seq:
			self._translate()
			return self.pep
		e = 'Cannot determine protein sequence because DNA is unknown for feature {0}'.format(self.ID)
		raise TranslateError(e)
	@pep.setter
	def pep(self,value):
		self._pep = value
	@pep.deleter
	def pep(self):
		self._pep = ''

	def set_attribute(self,key,value):
		if isisintance(value,basestring):
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
		frame2 = (1,self.seq[1:])
		frame3 = (2,self.seq[2:])
		for offset,frame in frame1,frame2,frame3:
			codons = re.findall('...',frame)
			starts = [index * 3 for index,codon in enumerate(codons) if codon == 'ATG']
			if not starts:
				continue
			stops = [index * 3 for index,codon in enumerate(codons) if codon in ['TAA','TGA','TAG']]
			if not stops:
				continue
			start = min(starts) + offset
			stop = min(stops) + offset + 2

			if stop - start > 0: # test for positive mRNA length
				orfs.append((start,stop))
		if not orfs:
			return False
		longest_orf = max(orfs,key = lambda x: x[1]-x[0])
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
	def _revcomp(self):
		comp = {'A':'T','T':'A','G':'C','C':'G','N':'N',
				'a':'t','t':'a','g':'c','c':'g','n':'n'}
		self.seq =  ''.join([comp[n] for n in self.seq[::-1]])

	def _translate(self):
		#print 'translate',self.ID,self._pep
		if not self.seq:
			raise NotImplementedError('Cannot determine self.pep without self.seq')
		if not len(self.seq) % 3 == 0:
			raise TranslateError('Length of sequence not dividable by 3')
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
		codons = re.findall('...',self.seq)
		pep = [p[codon.upper()] if 'N' not in codon.upper() else 'X' for codon in codons]
		self.pep = ''.join(pep)


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

class Gene(GffSubPart):
	"""
	Should inherit from GffSubPart
	"""
	pass

class Mrna(GffSubPart):
	"""
	Should inherit from GffSubPart
	"""
	def __init__(self):
		super(mRNA,self).__init__
		self._pep = ''
	@property
	def pep(self):
		return self._pep
	@pep.setter
	def pep(value):
		self._pep = value
	@pep.deleter
	def pep():
		self._pep = ''
	
class Cds(GffSubPart):
	"""
	Should inherit from GffSubPart
	"""
	pass

class Exon(GffSubPart):
	"""
	Should inherit from GffSubPart
	"""
	pass

class Utr(GffSubPart):
	"""
	Should inherit from GffSubPart
	"""
	pass

class Match(GffSubPart):
	"""
	Should inherit from GffSubPart
	"""
	pass

class MatchPart(GffSubPart):
	"""
	Should inherit from GffSubPart
	"""
	pass

