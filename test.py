#!/usr/bin/python
"""
docstring for testing
"""
__author__ = 'rensholmer'

from .parser import Parser
import tempfile
import os

_gff = 'scf1\ttest\tgene\t1\t105\t.\t+\t.\tID=test0\n\
scf1\ttest\tmRNA\t1\t105\t.\t+\t.\tID=test0.t1;Parent=test0\n\
scf1\ttest\tCDS\t1\t20\t.\t+\t.\tID=test0.t1.c1;Parent=test0.t1\n\
scf1\ttest\tCDS\t25\t105\t.\t+\t.\tID=test0.t1.c2;Parent=test0.t1\n\
scf1\ttest\tgene\t2\t105\t.\t+\t.\tID=test1\n\
scf1\ttest\tmRNA\t2\t105\t.\t+\t.\tID=test1.t1;Parent=test1\n\
scf1\ttest\tCDS\t2\t20\t.\t+\t.\tID=test1.t1.c1;Parent=test1.t1\n\
scf1\ttest\tCDS\t25\t105\t.\t+\t.\tID=test1.t1.c2;Parent=test1.t1\n\
scf2\ttest\tgene\t1\t15\t.\t+\t.\tID=test2\n\
scf2\ttest\tmRNA\t1\t15\t.\t+\t.\tID=test2.t1;Parent=test2\n\
scf2\ttest\tCDS\t1\t15\t.\t+\t.\tID=test2.t1.c1;Parent=test2.t1\n\
scf3\ttest\tgene\t1\t15\t.\t-\t.\tID=test3\n\
scf3\ttest\tmRNA\t1\t15\t.\t-\t.\tID=test3.t1;Parent=test3\n\
scf3\ttest\tCDS\t1\t15\t.\t-\t.\tID=test3.t1.c1;Parent=test3.t1\n\
scf4\ttest\tgene\t1\t105\t.\t-\t.\tID=test4\n\
scf4\ttest\tmRNA\t1\t105\t.\t-\t.\tID=test4.t1;Parent=test4\n\
scf4\ttest\tCDS\t1\t20\t.\t-\t.\tID=test4.t1.c1;Parent=test4.t1\n\
scf4\ttest\tCDS\t25\t105\t.\t-\t.\tID=test4.t1.c2;Parent=test4.t1\n\
scf4\ttest\ttest\t104\t105\t.\t-\t.\tID=test4.t1.c2.test;Parent=test4.t1.c2\n\
scf5\ttest\tgene\t1\t105\t.\t+\t.\tID=test5\n\
scf5\ttest\tmRNA\t1\t105\t.\t+\t.\tID=test5.t1;Parent=test5\n\
scf5\ttest\tCDS\t1\t20\t.\t+\t.\tID=test5.t1.c1;Parent=test5.t1\n\
scf5\ttest\tCDS\t25\t105\t.\t+\t.\tID=test5.t1.c2;Parent=test5.t1\n\
scf6\ttest\tgene\t1\t105\t.\t+\t.\tID=test6\n\
scf6\ttest\tmRNA\t1\t105\t.\t+\t.\tID=test6.t1;Parent=test6\n\
scf6\ttest\tCDS\t1\t20\t.\t+\t.\tID=test6.t1.c1;Parent=test6.t1\n\
scf6\ttest\tCDS\t25\t105\t.\t+\t.\tID=test6.t1.c2;Parent=test6.t1\n'

_seq = '>scf1\n\
AAATGGAGGAGGAGGAGGAGGAGGAGGAGGAGG\n\
AGGAGGAGGAGGAGGAGGAGAGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGTAAAGAGAGGAGA\n\
>scf2\n\
ATGGAGAGATGATGA\n\
>scf3\n\
TCATCACTCCTCCAT\n\
>scf4\n\
TCTCCTCTCTTTACTCCTCCTCCTCCTCCTCCTCCTCCTTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCATTT\n\
>scf5\n\
AAATGGAGTGAGAGGAGGAGGAGGAGGAGGAGGAGG\n\
AGGAGGAGGAGGAGGAGGAGAGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGTAAAGAGAGGAGA\n\
>scf6\n\
AAATGGAGGAGAGGAGGAGGAGGAGGAGGAGGAGG\n\
AGGAGGAGGAGGAGGAGGAGAGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGTAAAGAGAGGAGA\n'


def test():
	temp_gff = tempfile.NamedTemporaryFile('w+b', delete=False, prefix=os.getcwd() + '/temp/')
	temp_gff_name = temp_gff.name
	temp_gff.write(_gff)
	temp_gff.close()
	print 'parsing test gff'
	gff = parser(temp_gff_name)
	print 'succes'

	temp_seq = tempfile.NamedTemporaryFile('w', delete=False, prefix=os.getcwd() + '/temp/')
	temp_seq_name = temp_seq.name
	temp_seq.write(_seq)
	temp_seq.close()
	print 'parsing test fasta'
	with open(temp_seq_name, 'rU') as fh:
		gff.add_fasta(fh)
	print 'succes'

	print 'testing nesting'
	for f in gff.get_children('test0'):
		pass
		# print f
	print 'succes'

	print 'testing getitems'
	for i in gff.getitems(level=['CDS','mRNA']):
		for j in gff.getitems(level='CDS'):
			print i.match(j)
		#print i
	print 'succes'
	for i in gff.getitems(level='mRNA'):
		print i.seqid,i.ID
		l = []
		for x in gff.get_children(i,featuretype='CDS'):
			for y in gff.getitems(seqid=x.seqid,level='CDS'):
				if x.match(y):
					l.append(y.parents)
		print set([x for y in l for x in y])
		q = [y.parents for x in gff.get_children(i,featuretype='CDS') for y in gff.getitems(seqid=x.seqid,level='CDS') if x.match(y)]
		#print q
		p = [x for y in q for x in y]
		print set(p)


	quit()





	q = [j.parents for i in gff.getitems(level='CDS') for j in gff.getitems(seqid=i.seqid,level='CDS')]
	print q
	p = [x for y in q for x in y]
	print set(p)
	quit()
	gff.getseq(toplevel='mRNA',sublevel='CDS')

	transcripts = [x for x in gff.getitems(level='mRNA')]
	remove = []
	for t in transcripts:
		print '----'
		print t.strand,t.seqid,t.ID
		#print len(t.seq)
		#print t.seq
		#print 'fixing ORF'
		print t._find_orf()
		if not gff.fix_orf(t):
			print 'orf fix fail',t.ID
			#remove.append(t)

		#print t.seq
	gff.remove(remove,nested=True)
	m = gff.getitems(level='mRNA')
	for n in m:
		print n
	#print gff


# pp.pprint(gff.features)

if __name__ == '__main__':
	test()
