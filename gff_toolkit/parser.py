__author__ = 'rensholmer'

from gff import Gff
from gffsubpart import GffSubPart,TranslateError
from itertools import groupby
import uuid
import re

class Parser(object):
	"""
	Parser object gff3 formatted annotation files.
	The parse() method return the processed Gff object.
	Has several methods for known incorrectly formatted files: _ratt() and _manual(), which are called by specifying filetype
	"""
	_filetypes = ('standard','ratt','manual','interproscan')
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
				if not line.strip():
					continue
				if line[0] == '#':
					if 'FASTA' in line:
						return
					continue
				if self.filetype == 'interproscan':
					line = re.sub('; ',': ',line)
				parts = line.strip().split('\t')
				if self.limit != None:
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
				yield sub

	def _ratt(self):
		name_index_map = {}
		name_index_remove = set()
		child_map = {}
		
		for transcript in self.gff.getitems(featuretype='mRNA'):
			new_transcript_ID = transcript.attributes['locus_tag'][0]
			if transcript.ID == new_transcript_ID:
				continue

			child_map.setdefault(new_transcript_ID,[])

			name_index_map[new_transcript_ID] = transcript._key
			name_index_remove.add(transcript.ID)
			
			sub_counter = {}
			for sub in self.gff.get_children(transcript):
				if sub == transcript:
					continue
				sub_counter.setdefault(sub.featuretype,0)
				sub_counter[sub.featuretype] += 1
				
				sub.parents = [new_transcript_ID]
				new_sub_ID = '{0}.{1}{2}'.format(new_transcript_ID,sub.featuretype,sub_counter[sub.featuretype])
				
				child_map.setdefault(new_transcript_ID,[]).append(new_sub_ID)
				name_index_map[new_sub_ID] = sub._key
				name_index_remove.add(sub.ID)

				sub.ID = new_sub_ID
				sub.source = 'ratt'

			transcript.ID = new_transcript_ID
			transcript.source = 'ratt'
			transcript.attributes = {key:value for key,value in transcript.attributes.iteritems() if key in ('ID','Parent')}

			gene_ID = transcript.ID.split('.')[0]
			
			if not gene_ID in self.gff.name_index:
				gene_parts = transcript.gff_fields
				gene = GffSubPart(*gene_parts)
				gene.featuretype = 'gene'
				gene.ID = gene_ID
				gene.children = [transcript.ID]
				self.gff.update(gene)
			else:
				for gene in self.gff[gene_ID]:
					gene.children.append(transcript.ID)

			transcript.parents = [gene_ID]

		for new_name,_key in name_index_map.iteritems():
			self.gff.name_index[new_name] = [_key]
		for old_name in name_index_remove:
			self.gff.name_index.pop(old_name)
		for parent_id in child_map:
			for parent in self.gff[parent_id]:
				parent.children = child_map[parent_id]

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
	def _interproscan(self):
		type_dict = {}
		for feature in self.gff:
			for attribute,values in feature.attributes.iteritems():
				feature.attributes[attribute] = [value.strip('"') for value in values]
			if feature.featuretype == 'protein_match':
				#set parent attribute
				parent = feature.attributes['Target'][0]
				parent = parent.split()[0]
				feature.parents = [parent]
				#determine new ID based on parent name and source, this should result in unique IDs, opposed to what comes out of IPS
				source = feature.source
				type_dict.setdefault(parent,{}).setdefault(source,0)
				type_dict[parent][source] += 1
				type_count = type_dict[parent][source]
				ID = '{0}.{1}.{2}'.format(parent,source,type_count)
				self.gff.name_index[ID] = self.gff.name_index.pop(feature.ID)
				feature.ID = ID

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

	def _get_nested_ratt_ids(self,ID):
		IDs = set()
		for lower_ID in self._get_lower_ratt_ids(ID):
			IDs.add(lower_ID)
			yield lower_ID
		for higher_ID in self._get_higher_ratt_ids(ID):
			if higher_ID not in IDs:
				yield higher_ID

	def _get_lower_ratt_ids(self,ID):
		yield ID
		if ID[-4] != '.':
			return
		for sub_ID in self._get_lower_ratt_ids(ID[:-2]):
			yield sub_ID

	def _get_higher_ratt_ids(self,ID):
		yield ID
		counter = int(ID[-1])
		new_ID = '{0}.{1}'.format(ID,counter + 1)
		if new_ID not in self.gff:
			return
		for super_ID in self._get_higher_ratt_ids(new_ID):
			yield super_ID

	def _remove_fragmented_ratt(self):
		remove = set()
		ID_map = {}
		for transcript in self.gff.getitems(featuretype='mRNA'):
			if transcript.ID in remove:
				continue
			if transcript.ID[-4] == '.':
				#get all fragments
				fragments = (ff for f in self._get_nested_ratt_ids(transcript.ID) for ff in self.gff[f])
				#sort fragments by CDS length
				sorted_fragments = sorted(fragments,key = lambda x:len(x.seq))
				#remove the longest fragment from the list so it doesnt get removed
				longest = sorted_fragments.pop()
				#remove all the other fragments	
				remove |= {f.ID for f in sorted_fragments}#set(sorted_fragments)
				#make sure longest fragment keeps the original name
				transcript_ID = [ID for ID in self._get_nested_ratt_ids(longest.ID) if ID[-4] != '.'][0]
				ID_map[longest.ID] = transcript_ID
				#set ID and Parent attribute of nested CDS features
				sub_counter = {}
				for sub in self.gff.get_children(longest):
					if sub == longest:
						continue
					sub.attributes['Parent'] = [transcript_ID]
					sub_counter.setdefault(sub.featuretype,0)
					sub_counter[sub.featuretype] += 1
					sub_ID = '{0}.{1}{2}'.format(transcript_ID,sub.featuretype,sub_counter[sub.featuretype])
					ID_map[sub.ID] = sub_ID
		self.gff.remove(remove)
		for old_ID,new_ID in ID_map.iteritems():
			assert old_ID not in remove,(old_ID,new_ID)
			for feature in self.gff[old_ID]:
				feature.ID = new_ID

	def parse(self):
		self.remove_fragmented_ratt = False
		for subpart in self._readlines():
			#if self.filetype == 'interproscan':
			#	if subpart.featuretype == 'protein_match':
			#		parent = subpart.attributes['Target'][0]
			#		parent = parent.split()[0]
			#		subpart.parents = [parent]
			self.gff.update(subpart)
		if self.filetype == 'ratt':
			self.remove_fragmented_ratt = True
			self.gff.set_children()
			self._ratt()
		elif self.filetype == 'manual':
			self._manual()
			self.gff.set_children()
		elif self.filetype == 'interproscan':
			self._interproscan()
			self.gff.set_children()
		else:
			pass
		if self.filetype != 'ratt':
			self.gff.set_children()
		if self.fasta_file:
			self.gff.add_fasta(self.fasta_file)
		if self.remove_fragmented_ratt:
			self._remove_fragmented_ratt()
		if self.remove_noncoding:
			self._remove_noncoding()
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

