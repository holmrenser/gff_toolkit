__author__ = 'rensholmer'

from gffsubpart import GffSubPart
from itertools import groupby #groupby for fasta parsing
#import blist #use sorted list object for indexing
import pprint
from intervaltree import IntervalTree, Interval

class Gff(object):
	"""
	Work in progess: holds GffSubParts object
	"""
	_combos = [{'gene':{'mRNA':['CDS','exon','five_prime_UTR','three_prime_UTR']}},
				{'match':'match_part'},
				{'protein_match':'match_part'}]
	_featuretypes = ['gene','mRNA','CDS','exon','five_prime_UTR','three_prime_UTR',
	          'match','protein_match','transcript_match','match_part',
	           'biological_region']
	def __init__(self,*args,**kwargs):
		"""
		Fire it up!
		"""
		self.features = {} #dict with {uniqueID1:feature1,uniqueID2:feature2,...}   OLD:{seqid:[GffSubPart1,GffSubPart2,etc]}
		self.seq = {} #sequencedict with {header:seq}
		#self._typecounts = {l:0 for l in self._featuretypes}
		self._uniqueID = 0
		self.filename = kwargs.get('filename','')
		self.name_index = {} #dict with {ID:set(uniqueID1,uniqueID2,),..} to access features based on non unique ID
		self.position_index = {}#bx.intervaltree
		self.type_index = {l:set() for l in self._featuretypes} #dict with {featuretype:set(uniqueID1,uniqueID2,),...} to access features based on featuretype
	def __iter__(self):
		"""
		Lets loop over all subfeatures
		"""
		for f in self.features.values():
			yield f
	def __getitem__(self,key):
		"""
		Allow square bracket access based on ID, like this: gff[ID]
		"""
		if not isinstance(key,basestring):
			e = 'Object of type {0} is not a valid key'.format(type(key))
			raise TypeError(e)
		for uniqueID in self.name_index[key]:
			yield self.features[uniqueID]
	def __str__(self):
		return self.filename #temporary
	def __repr__(self):
		return self.__str__()

	@property
	def uniqueID(self):
		self._uniqueID += 1
		return self._uniqueID
	
	def typecounts(self):
		return pprint.pprint({k:len(v) for k,v in self.type_index.iteritems() if len(v) > 0})
	
	def stringify(self):
		"""
		Returns entire object as gff formatted string
		"""
		string = []
		for sub in self:
			string.append(sub.stringify())
		return ''.join(string)

	def getitems(self,seqid=None,start=None,end=None,strand=None,featuretype=None):
		"""
		Select based on seqid (scf0002,chr1,etc.) and featuretype (gene,mRNA,etc.)
		"""
		if featuretype != None:
			if isinstance(featuretype,basestring) and featuretype in self._featuretypes:
				featuretype = [featuretype]
			elif isinstance(featuretype,(list,tuple)):
				pass
			else:
				e = '{0} is not a valid type for featuretype'.format(type(featuretype))
				raise TypeError(e)
		if seqid == None:
			if start != None:
				e = 'Can not provide start when no seqid is provided'
				raise NotImplementedError(e)
			elif end != None:
				e = 'Can not provide end when no seqid is provided'
				raise NotImplementedError(e)
			elif strand != None:
				e = 'Can not provide strand when no seqid is provided'
				raise NotImplementedError(e)
			else:
				seqids = self.position_index.keys()
				#for f in featuretype:
				#	for uniqueID in self.type_index[f]:
				#		sub = self.features[uniqueID]
				#		if f == None or sub.featuretype in f:
				#			yield sub
		elif not isinstance(seqid,basestring):
			e = '{0} is not a valid type for seqid'.format(type(seqid))
			raise TypeError(e)
		else:
			seqids = [seqid]

		for seqid in seqids:
			if seqid not in self.position_index:
				return #False
			if (start == None or end == None) and start != end:
				raise Exception()
			if start == None:
				start = 0
			if end == None:
				end = 10e10
			for f in self.position_index[seqid].search(start,end):
				sub = self.features[f.data['ID']]
				if featuretype == None or sub.featuretype in featuretype:
					yield sub

	def getseq(self,feature=None,subfeaturetype=None,topfeaturetype=None):
		if isinstance(feature,basestring):
			features = self[feature]
		elif isinstance(feature,GffSubPart):
			features = [feature]
		elif feature == None:
			features = self.getitems(featuretype=topfeaturetype)
		else:
			raise TypeError('feature is not of type GffSubPart, or String')

		for f in features:
			if f.seq:
				print f.ID,'allready has seq!'
				f.seq = ''
			if f.strand == '+':
				reverse = False
			else:
				reverse = True
			children = self.get_children(f,featuretype=subfeaturetype)
			children = sorted([c for c in children],key = lambda x: x.get_start(),reverse=reverse)
			for c in children:
				c.seq = self.seq[c.seqid][c.start-1:c.end]
				if reverse:
					c._revcomp()
				f.seq += c.seq

	def remove(self,key,nested=True):
		"""
		Not stable, sometimes doesn't remove things?
		"""
		if isinstance(key,basestring):
			keys = self[key]
		elif isinstance(key,GffSubPart):
			keys = [key]
		elif isinstance(key,(list,tuple)):
			keys = []
			for k in key:
				if isinstance(k,basestring):
					keys += [x for x in self[k]]
				elif isinstance(k,GffSubPart):
					keys.append(k)
				else:
					e = '{0} is not a valid type'.format(k)
					raise TypeError(e)
		else:
			e = '{0} is not a valid type'.format(key)
			raise TypeError(e)
		if nested:
			nestedkeys = []
			for key in keys:
				nk = [x for x in self.get_children(key,reverse=True)]
				nestedkeys += nk
		for k in set(nestedkeys):
			self.type_index[k.featuretype].remove(k._key)
			self.position_index[k.seqid].discard(Interval(k.start,k.end,{'ID':k._key}))
			self.name_index[k.ID].remove(k._key)
			if len(self.name_index[k.ID]) == 0:
				del self.name_index[k.ID]
			del self.features[k._key]
		return True

	def names(self,seqid=None, featuretype=None):
		"""
		:return:
		"""
		for f in self.getitems(seqid=seqid,featuretype=featuretype):
			yield f.ID

	def update(self,subfeature):
		"""
		:param subfeature: GffSubFeature object
		:return:
		"""
		if not isinstance(subfeature,GffSubPart):
			raise NotImplementedError()
		ID = self.uniqueID

		subfeature._key = ID

		self.features[ID] = subfeature
		self.name_index.setdefault(subfeature.ID,set()).add(ID)
		self.type_index[subfeature.featuretype].add(ID)

		interval = Interval(subfeature.start,subfeature.end,{'ID':ID})
		self.position_index.setdefault(subfeature.seqid,IntervalTree()).add(interval)
		subfeature.container = self

	def make_index(self):
		pass

		#for subfeature in self:
		#	if subfeature.seqid not in self.index[subfeature.strand]:
		#		self.index[subfeature.strand][subfeature.seqid] = blist.sortedlist()
		#	self.index[subfeature.strand][subfeature.seqid].update((((subfeature.start,subfeature.end),subfeature.ID),))

	def set_children(self):
		"""
		Sets the children attribute of all subfeatures
		:return:
		"""
		for sub in self:
			for p_name in sub.parents:
				for p_obj in self[p_name]:
					if sub.ID not in p_obj.children:
						p_obj.children.append(sub.ID)

	def get_children(self,key,reverse=False,featuretype=None):
		"""
		:param key: subfeature ID or subfeature object
		:param reverse: reverses return order. I.e.: reverse=True return CDS->mRNA->gene. reverse=False returns gene->mRNA->CDS
		:param featuretype: string or list with featuretypes to be returned
		:return: nested generator of subfeature objects
		"""
		if isinstance(key,GffSubPart):
			keys = [key]
		elif isinstance(key,basestring):
			keys = [k for k in self[key]]
		elif isinstance(key,(list,tuple)):
			for k in key:
				keys = []
				if isinstance(k,GffSubPart):
					keys.append(k)
				elif isinstance(k,basestring):
					keys += [x for x in self[k]]
				else:
					e = '{0} is not a valid key for Gff.get_children()'.format(k)
					raise TypeError(e)
			keys = [key]
		else:
			print type(key),key
			e = '{0} is not a valid key for Gff.get_children()'.format(key)
			raise TypeError(e)
		if featuretype != None:
			if isinstance(featuretype,basestring):
				if featuretype in self._featuretypes:
					featuretype = [featuretype]
				else:
					e = '{0} is not a valid featuretype'.format(featuretype)
					raise TypeError(e)
			elif isinstance(featuretype,(list,tuple)):
				pass
			else:
				e = '{0} is not a valid type for featuretype'.format(type(featuretype))
				raise TypeError(e)

		for k in keys:
			if not reverse and (featuretype == None or k.featuretype in featuretype):
				yield k
			for child in k.children:
				for nested_child in self.get_children(child):
					if featuretype == None or nested_child.featuretype in featuretype:
						yield nested_child
			if reverse and (featuretype == None or k.featuretype in featuretype):
				yield k
	#def get_overlap(self,subfeature,stranded=True):
	def get_overlap(self,seqid,start,end,strand=None,featuretype=None):
		"""
		Returns all subfeatures of featuretype that 
		"""
		if seqid not in self.features:
			#trick to end generator if searching for unknown seqid
			return
		if start > end:
			e = 'start cannot be bigger than end, start: {0} -- end: {1}'.format(start,end)
			raise ValueError(e)
		if strand:
			if strand not in ['+','-']:
				e = '{0} is not a valid strand'.format(strand)
			strands = [strand]
		else:
			strands = ['+','-']
		
		for strand in strands:
			found = 0
			if seqid not in self.index[strand]:
				continue
			for i in self.index[strand][seqid]:
				if found > 0:
					break
				for ii in self[i[1]]:
					if featuretype == None:
						pass
					elif ii.featuretype != featuretype:
						continue
					if start <= ii.start <= end or start <= ii.end <= end:
						found += 1
						yield ii
					elif ii.start <= start <= ii.end or ii.start <= end <= ii.end:
						found += 1
						yield ii
					elif found > 0:
						break
	def add_fasta(self,filename):
		"""
		:param filehandle: fasta formatted DNA sequence file
		:return:
		"""
		with open(filename,'rU') as filehandle:
			faiter = (x[1] for x in groupby(filehandle, lambda line: line[0] == '>'))
			for header in faiter:
				header = header.next()[1:].strip()
				seq = ''.join(s.strip() for s in faiter.next())
				self.seq[header] = seq
	def write_tbl(self):
		dic = {}
		for x in self.getitems(featuretype='gene'):
			string = ''
			for y in self.get_children(x,featuretype=['gene','mRNA','CDS']):
				string += y.stringify(filetype='tbl') + '\n'
			dic.setdefault(x.seqid,[]).append(string)
			#print string
		for s in dic:
			print '>Feature {0}'.format(s)
			for t in dic[s]:
				print t
		#for uniqueID in self.type_index['gene']:
		#	gene = self.features[uniqueID]
		#	print gene.stringify(filetype='tbl')

	def _range_map(self,subfeature):
		"""
		:param subfeature: GffSubFeature object with children
		:return: dictionary with {subfeature.seq.coordinate : scaffold.seq.coordinate}
		"""
		pos = 0
		cds = []
		range_map = {}
		children = self.get_children(subfeature,featuretype='CDS')
		if subfeature.strand == '+':
			reverse = False
			step = 1
		else:
			reverse = True
			step = -1
		children = sorted([c for c in children],key = lambda x:x.start,reverse=reverse)
		for c in children:
			mrna_range = range(pos,pos+c.end-c.start+1)
			genome_range = range(c.get_start(),c.get_end()+step,step)
			x = {a:b for a,b in zip(mrna_range,genome_range)}
			range_map.update(x)
			pos += c.end - c.start  + 1
		return range_map

	def _change_cds(self,subfeature,genome_orf):
		"""
		:param self:
		:param subfeature: subfeature with CDS children
		:param orf: tuple with genome start and genome stop of ORF in subfeature mRNA
		:return: True if succesful change else False
		"""
		found_start = False
		found_stop = False
		remove = []
		children = self.get_children(subfeature,featuretype='CDS')
		forward = subfeature.strand == '+'
		if forward:
			s1 = genome_orf[0]
			s2 = genome_orf[1]
			children = sorted(children,key=lambda x: x.start)
		else:
			s1 = genome_orf[1]
			s2 = genome_orf[0]
			children = sorted(children,key=lambda x: x.end, reverse=True)
		for c in children:
			c1 = c.get_start()
			c2 = c.get_end()

			forward_start =  c1 <= s1 <= c2
			forward_stop = c1 <= s2 <= c2
			reverse_start = c1 >= s1 >= c2
			reverse_stop = c1 >= s2 >= c2

			#change start
			if (forward and forward_start) or (not forward and reverse_start) and not found_stop:
				print 'FOUND START',(c1,s1,c2)
				c.set_start(s1)
				print c.start,c.end
				print c.seq[c.start-1:c.end]
				found_start = True
			#change stop and continue loop from the top
			if (forward and forward_stop) or (not forward and reverse_stop) and found_start and not found_stop:
				c.set_end(s2)
				found_stop = True
				#remove zero length exons
				if c.start == c.end:
					remove.append(c)
				continue
			#remove zero length exons
			if c.start == c.end:
				remove.append(c)
			#remove exon before start
			if not found_start:
				remove.append(c)
			#remove exon after stop
			if found_stop and found_start:
				remove.append(c)
		self.remove(remove)
		self.getseq(subfeature,topfeaturetype=subfeature.featuretype,subfeaturetype='CDS')
		if len(subfeature.seq) < 6 or len(subfeature.seq) % 3 != 0:
			return False
		else:
			return True
	def fix_orf(self,subfeature,min_len=0):
		orf = subfeature._find_orf()
		if not orf:
			return False
		if max(orf) - min(orf) < min_len:
			return False
		range_map = self._range_map(subfeature)
		if not range_map:
			return False
		genome_orf = range_map[orf[0]],range_map[orf[1]] #self._map_start_stop(orf,range_map)
		if subfeature.strand == '-':
			genome_orf = (genome_orf[1],genome_orf[0])
		if not genome_orf:
			return False
		change = self._change_cds(subfeature,genome_orf)
		print 'changed',change
		return change

