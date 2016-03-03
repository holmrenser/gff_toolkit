__author__ = 'rensholmer'

from gffsubpart import GffSubPart
from itertools import groupby #groupby for fasta parsing
import pprint
from intervaltree import IntervalTree, Interval

class Gff(object):
	"""
	Work in progess: holds GffSubParts object
	"""
	_combos = [{'gene':{'mRNA':['CDS','exon','five_prime_UTR','three_prime_UTR']}},
				{'match':'match_part'},
				{'protein_match':'match_part'}]
	_featuretypes = ('gene','mRNA','CDS','exon','five_prime_UTR','three_prime_UTR',
	          'match','protein_match','transcript_match','match_part',
	           'biological_region','polypeptide')
	def __init__(self,*args,**kwargs):
		"""
		Fire it up!
		"""
		self.features = {} #dict with {_key1:feature1,_key2:feature2,...}   OLD:{seqid:[GffSubPart1,GffSubPart2,etc]}
		self.seq = {} #sequencedict with {header:seq}
		self._removedkeys = set()
		self._uniqueID = 0 #unique IDs for subfeature._key
		self.filename = kwargs.get('filename','')
		self.name_index = {} #dict with {ID:set(_key1,_key2,),..} to access features based on non unique ID
		self.position_index = {}#Intervaltree
		self.type_index = {l:set() for l in self._featuretypes} #dict with {featuretype:set(uniqueID1,uniqueID2,),...} to access features based on featuretype
	def __iter__(self):
		"""
		Iterate over all values
		"""
		for item in self.features.values():
			yield item
	def __getitem__(self,key):
		"""
		Allow square bracket access based on ID, like this: gff[ID]
		"""
		if not isinstance(key,basestring):
			e = 'Object of type {0} is not a valid key'.format(type(key))
			raise TypeError(e)
		for uniqueID in self.name_index[key]:
			yield self.features[uniqueID]
	def __len__(self):
		return len(self.features.values())
	def __str__(self):
		return self.filename #temporary
	def __repr__(self):
		return self.__str__()
	def __add__(self,other):
		"""
		Add two Gff objects together
		"""
		new_gff = Gff()
		for self_obj in self:
			new_gff.update(self_obj)
		for other_obj in other:
			new_gff.update(other_obj)
		#new_gff.set_children()
		return new_gff

	def split(self):
		for seqid in self.position_index:
			new_gff = Gff()
			new_gff.seq[seqid] = self.seq[seqid]
			for sub in self.getitems(seqid=seqid):
				new_gff.update(sub)
			yield new_gff

	@property
	def uniqueID(self):
		self._uniqueID += 1
		return self._uniqueID
	
	def typecounts(self):
		return {k:len(v) for k,v in self.type_index.iteritems() if len(v) > 0}
	
	def stringify(self):
		"""
		Args:
			None
		Returns:
			Entire Gff object as gff formatted string
		"""
		string = []
		for sub in self:
			string.append(sub.stringify())
		return ''.join(string)

	def getitems(self,seqid=None,start=None,end=None,strand=None,featuretype=None):
		"""
		Args:
			seqid: blabla
		Returns:
			GffSubPart
		"""
		coords = (start,end)
		if coords.count(None) == 1:
			e = '({0},{1}) are not valid coordinates'.format(*coords)
			raise Exception(e)
		if featuretype != None:
			if isinstance(featuretype,basestring) and featuretype in self._featuretypes:
				featuretype = [featuretype]
			elif isinstance(featuretype,(list,tuple)):
				pass
			else:
				e = '{0} is not a valid type for featuretype'.format(type(featuretype))
				raise TypeError(e)
		if seqid == None:
			if coords.count(None) != 2:
				e = 'Can not provide end when no seqid is provided'
				raise NotImplementedError(e)
			elif strand != None:
				e = 'Can not provide strand when no seqid is provided'
				raise NotImplementedError(e)
			else:
				seqids = self.position_index.keys()

		elif not isinstance(seqid,basestring):
			e = '{0} is not a valid type for seqid'.format(type(seqid))
			raise TypeError(e)
		else:
			seqids = [seqid]

		for seqid in seqids:
			if seqid not in self.position_index:
				continue #False
			if start == None and end == None:
				subs = (self.features[_key] for _key in (interval.data['_key'] for interval in self.position_index[seqid].items()) if _key not in self._removedkeys)
			else:
				subs = (self.features[_key] for _key in (interval.data['_key'] for interval in self.position_index[seqid].search(start,end)) if _key not in self._removedkeys)
			for sub in subs:
				if (featuretype == None or sub.featuretype in featuretype) and (strand == None or sub.strand == strand) and not (sub.start == end or sub.end == start):
					yield sub
		return

	def getseq(self,feature=None,subfeaturetype=None,topfeaturetype=None):
		"""
		This is replaced by a combination of properties on GffSubPart: seq,pep and siblings
		"""
		if isinstance(feature,basestring):
			features = self[feature]
		elif isinstance(feature,GffSubPart):
			features = [feature]
		elif feature == None:
			features = self.getitems(featuretype=topfeaturetype)
		else:
			raise TypeError('feature is not of type GffSubPart, or String')

		for feature in features:
			if feature.seq:
				feature.seq = ''
			if feature.strand == '+':
				reverse = False
			else:
				reverse = True
			children = self.get_children(feature,featuretype=subfeaturetype)
			children = sorted([c for c in children],key = lambda x: x.get_start(),reverse=reverse)
			for index,cds in enumerate(children):
				#if cds.end - cds.start == 1 and this is the last cds, only select one nucleotide, otherwise two
				if cds.start + 1 == cds.end and index == len(children) - 1:
					cds.seq = self.seq[cds.seqid][cds.start-1]
				else:
					cds.seq = self.seq[cds.seqid][cds.start-1:cds.end]
				if reverse:
					cds._revcomp()
				#print cds.seq
				feature.seq += cds.seq

	def remove(self,key,nested=True):
		"""
		Args:
			key: string or GffSubPart or list
			nested: bool
		Returns:
			None
		"""
		if isinstance(key,basestring):
			keys = self[key]
		elif isinstance(key,GffSubPart):
			keys = [key]
		elif hasattr(key,'__iter__'):#isinstance(key,(list,tuple)):
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

		for key in keys:
			if nested:
				nestedkeys = list(self.get_children(key,reverse=True))
			else:
				nestedkeys = [key]
			remove_parents = set()
			for nestedkey in nestedkeys:
				self.type_index[nestedkey.featuretype].remove(nestedkey._key)
				self.position_index[nestedkey.seqid].discardi(nestedkey.start,nestedkey.end,{'_key':nestedkey._key})
				#self.position_index[k.seqid].remove(Interval(k.start,k.end,{'_key':k._key}))
				self.name_index[nestedkey.ID].remove(nestedkey._key)
				if len(self.name_index[nestedkey.ID]) == 0:
					del self.name_index[nestedkey.ID]
				parents = {pp for p in nestedkey.parents for pp in self[p]}
				for p in parents:
					if nestedkey.ID not in self.name_index:
						p.children.remove(nestedkey.ID)
					if not p.children and p not in nestedkeys:
						remove_parents.add(p)
				del self.features[nestedkey._key]
				self._removedkeys.add(nestedkey._key)
			if remove_parents:
				self.remove(list(remove_parents),nested=False)

		return True

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

		interval = Interval(subfeature._start,subfeature.end,{'_key':ID})
		self.position_index.setdefault(subfeature.seqid,IntervalTree()).add(interval)
		subfeature.container = self

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

	def get_parents(self,key,reverse=True,featuretype=None):
		"""
		"""
		pass

	def get_children(self,key,reverse=False,featuretype=None,seen=None):
		"""
		:param key: subfeature ID or subfeature object
		:param reverse: reverses return order. I.e.: reverse=True return CDS->mRNA->gene. reverse=False returns gene->mRNA->CDS
		:param featuretype: string or list with featuretypes to be returned
		:return: nested generator of subfeature objects

		TODO: add something that prevents double yields
		"""
		if seen == None:
			seen = set()
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
		#print [s.ID for s in seen]
		for k in keys:
			if k._key in self._removedkeys:
				continue
			if not reverse and (featuretype == None or k.featuretype in featuretype) and k not in seen:
				seen.add(k)
				yield k
			for child in k.children:
				for nested_child in self.get_children(child,seen=seen):
					if featuretype == None or nested_child.featuretype in featuretype and k not in seen:
						seen.add(nested_child)
						yield nested_child
			if reverse and (featuretype == None or k.featuretype in featuretype) and k not in seen:
				seen.add(k)
				yield k

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
		"""
		Args: None
		Return:
			.tbl formatted string
		"""
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
		children = sorted((c for c in children),key = lambda x:x.start,reverse=reverse)
		for c in children:
			mrna_range = range(pos,pos+c.end-c.start+1)
			genome_range = range(c.get_start(),c.get_end()+step,step)
			assert len(mrna_range) == len(genome_range)
			x = {a:b for a,b in zip(mrna_range,genome_range)}
			range_map.update(x)
			pos += c.end - c.start + 1
		range_map[max(range_map) + 1] = range_map[max(range_map)] + step
		return range_map

	def _change_cds(self,subfeature,genome_orf):
		"""
		:param self:
		:param subfeature: subfeature with CDS children
		:param orf: tuple with genome start and genome stop of ORF in subfeature mRNA
		:return: True if succesful change else False
		"""
		new_start = genome_orf[0]
		new_stop = genome_orf[1]
		forward = subfeature.strand == '+'
		reverse = subfeature.strand == '-'
		remove = set()
		for cds in self.get_children(subfeature,featuretype='CDS'):
			found_start = False
			found_stop = False
			if cds.end < new_start or cds.start > new_stop:
				#print 'remove'
				remove.add(cds)
				continue
			found_start = cds.start <= new_start <= cds.end
			found_stop = cds.start <= new_stop <= cds.end
			if found_start:
				print 'found start, stop ==',found_stop
				if found_stop:
					if reverse:
						new_start += 1
				cds.start = new_start
				cds.phase = 0
			if found_stop:
				print 'found stop, start ==',found_start
				if found_start:
					if forward:
						new_stop -= 1
				cds.end = new_stop

		'''
		children = self.get_children(subfeature,featuretype='CDS')
		forward = subfeature.strand == '+'
		if forward:
			new_start = genome_orf[0]
			new_end = genome_orf[1] #- 1
			children = sorted(children,key=lambda x: x.get_start())
		else:
			new_start = genome_orf[1]
			new_end = genome_orf[0] #+ 1
			children = sorted(children,key=lambda x: x.get_start(), reverse=True)

		print 'GENOME ORF',genome_orf
		for c in children:
			old_start = c.get_start()
			old_end = c.get_end()


			forward_start = old_start <= new_start <= old_end
			reverse_start = old_start >= new_start >= old_end
			forward_stop = old_start <= new_end <= old_end
			reverse_stop = old_start >= new_end >= old_end
			#print '---'
			#change start
			if (forward and forward_start) or (not forward and reverse_start) and not found_stop:
				print 'start',old_start,(new_start,old_end)
				found_start = True


				#fix problem with exons of length 1
				if abs(new_start - old_end) == 1:
				#	if forward:
				#		pass
					#new_start -= 1
					#c.set_end(old_end + 1)
					print 'this just happened',old_start,(new_start,old_end)

				#if old_end == new_start:
				#	if forward:
				#		pass
				#		#new_start -= 1
				#	else:
				#		pass
				#		#new_start += 1
				#	print 'this just happened',old_start,(new_start,old_end)
				c.set_start(new_start)
			#change stop and continue loop from the top
			if (forward and forward_stop) or (not forward and reverse_stop) and found_start and not found_stop:
				print 'end',(old_start,new_end),old_end
				found_stop = True
				#fix problem if start and stop are in the same exon
				if forward_start:
					new_end -= 1
				if reverse_start:
					new_end += 1
				#if abs(old_start - old_end) == 1:
				#	if forward:
				#		pass
						#new_end += 1
				#	else:
				#		pass
						#new_end -= 1
				#remove exons of length 0
				if old_start == new_end:
					remove.append(c)
				else:
					c.set_end(new_end)
				continue
			#remove exon before start
			if not found_start:
				remove.append(c)
			#remove exon after stop
			if found_stop and found_start:
				remove.append(c)
		'''
		#remove non coding exons
		print [x.ID for x in remove]
		if remove:
			self.remove(remove,nested=False)
		#recalculate transcript coding sequence
		self.getseq(subfeature,topfeaturetype=subfeature.featuretype,subfeaturetype='CDS')
		#assertions to make sure ORF is fixed
		assert len(subfeature.seq) >= 6,(subfeature.ID,len(subfeature.seq),subfeature.seq)
		assert len(subfeature.seq) % 3 == 0,(subfeature.ID,subfeature.strand,len(subfeature.seq),subfeature.seq)
		assert subfeature.seq[:3] in ('CTG','ATG'),(subfeature.ID,subfeature.seq[:3],subfeature.seq)
		assert subfeature.seq[-3:] in ('TGA','TAA','TAG')
		assert '*' not in subfeature.pep[-1],(subfeature.ID,subfeature.pep)
		return True

	def fix_orf(self,subfeature,starts=('ATG','CTG'),stops=('TAA','TGA','TAG'),min_len=6):
		"""
		Finds longest ORF in spliced transcript.

		Args:
			subfeature: class: GffSubPart()
			starts: list/tuple with start codons
			stops: list/tuple with stop codons
			min_len: minimum ORF length
		Returns:
			True if ORF is found
		"""
		orf = subfeature._find_orf()
		if not orf:
			return False
		if max(orf) - min(orf) < min_len:
			return False
		range_map = self._range_map(subfeature)
		if not range_map:
			return False
		#print orf
		#print range_map
		genome_orf = range_map[orf[0]],range_map[orf[1]] #self._map_start_stop(orf,range_map)
		#print genome_orf
		if subfeature.strand == '-':
			genome_orf = (genome_orf[1],genome_orf[0])
		if not genome_orf:
			return False
		change = self._change_cds(subfeature,genome_orf)
		return change

