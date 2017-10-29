#!/usr/bin/env python

import sys
from ete2 import Tree
import numpy as np
import re

############## define functions#############

# Given initial state n0, transition matrix mat, time t, return a vector containing probability of each state after evolving.
def evo(n0,t,mat):
	n = np.dot(n0,np.linalg.matrix_power(mat,t))
	return n

# check if node i is ancestral or decendant of node j, if not, return True, so that the two nodes are in independent branches
def checkNode(nodei,nodej,tree):
    if nodei.up == nodej.up:
        return False
    anci = []
    node = nodei
    while node.up:
        anci.append(node)
        node = node.up
    if nodej in anci:
        return False
    ancj = []
    node = nodej
    while node.up:
        ancj.append(node)
        node = node.up
    if nodei in ancj:
        return False
    return True
############################################

if len(sys.argv)<3:
	print "Usage: probCalP.py <genelist> <species_tree> <prefix>"
	sys.exit(2)


GENE = open(sys.argv[1],'r')
sptree = Tree(sys.argv[2],format=3)
#PAIR = open(sys.argv[2],'r')
prefix = sys.argv[3] # for prefix input anything you want to appear at start of output files
OUT = open(prefix+"_results.genes",'w')
# OUTSITE = open(prefix+"_results.sites",'w')

pairlist = [] #list of branch pairs to look at
pairres = []
flag = 0

# initialize JTT matrix and JTT frequency vector
jmat = np.loadtxt("JTT.mat")
jmat.shape = 20,20
jmat = jmat/100
jfreq = np.loadtxt("jtt.freq")

aaL = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

cnt = 0
for gene in GENE:
	gene = gene.rstrip()

	# initialize tree
	TREE = open('../02_Trees/'+gene+'.ntree','r')
	treeline = TREE.readline().rstrip()
	tree = Tree(treeline,format=3)
	if flag==0:
		nodelist = [x.name for x in tree.get_descendants()]
		nnode = len(nodelist)
		for i in range(nnode):
			for j in range(i):
				nodei = tree&nodelist[i]
				nodej = tree&nodelist[j]
				if checkNode(nodei,nodej,tree):
					pairlist.append([nodelist[i],nodelist[j]])
		# store results for each pair of branches
		pairres = [[0,0,sptree.get_distance(x[0],x[1])] for x in pairlist]
		flag = 1

	# read from rate/freq files into numpy arrays for later use
	rateL = np.loadtxt('../03_rate/'+gene+'.rat') # rate of each site
	freqL = np.loadtxt('../04_siteFreq/'+gene+'.freq')
	NODE = open('../01_NodeSeq/'+gene+'.nodes','r')
	nodeD = {}
	for line in NODE:
		col = line.rstrip().split()
		if col[0][:4]=='node':
			nodeD[col[0][4:]] = col[1]
		else:
			nodeD[col[0]] = col[1]
	for pairL in pairlist:
		taxus1,taxus2 = pairL
		print >>OUT,'%s\t%s\t'%(gene,'-'.join(pairL)),
		print '%s\t%s\t'%(gene,'-'.join(pairL)),
		
		# get branch lengths and decide which ancestral sequence should be used
		node1 = tree&taxus1
		node2 = tree&taxus2
		node1a = node1.up
		node2a = node2.up
		ancnode = tree.get_common_ancestor(node1,node2)
		if node1a == node2a:
			sys.exit(2)
		t0 = tree.get_distance(node1a,ancnode)
		t1 = tree.get_distance(node1,node1a)
		t2 = tree.get_distance(node2a,ancnode)
		t3 = tree.get_distance(node2,node2a)
		# get node sequences
		node1seq = nodeD[taxus1]
		node2seq = nodeD[taxus2]
		node1aseq = nodeD[node1a.name]
		node2aseq = nodeD[node2a.name]
		ancseq = nodeD[ancnode.name]

		cntL = []
		probL = []
		
		for pos in range(len(rateL)):
			# print >>OUTSITE,"%s\t%d\t%s\t"%(gene,pos,'-'.join(pairL)),
			# count parallel/convergent change
			n1 = node1seq[pos]
			n2 = node2seq[pos]
			n1a = node1aseq[pos]
			n2a = node2aseq[pos]
			if n1==n2 and n1a!=n1 and n2a!=n2 and n1a==n2a:
				cntL.append(1)
				# print >>OUTSITE,"1\t",
			else:
				cntL.append(0)
				# print >>OUTSITE,"0\t",

			# calculate parallel/convergent probability
			T = [int(np.around(rateL[pos]*x*10000)) for x in [t0,t1,t2,t3]]
			ancaa = ancseq[pos]
			anc = np.array([1 if ancaa==aaL[i] else 0 for i in range(20)])
			freq = freqL[pos,:]
			mul = freq/jfreq
			mul.shape = 1,20
			mulmat = mul.repeat(20,axis=0)
			jfmat = jmat*mulmat
			for i in range(20):
				jfmat[i,i] = 1-(np.sum(jfmat[i,:])-jfmat[i,i])

			n1aprob = evo(anc,T[0],jfmat)
			n2aprob = evo(anc,T[2],jfmat)
			Id = np.identity(20)
			mat1 = evo(Id,T[1],jfmat)
			mat2 = evo(Id,T[3],jfmat)
			for i in range(20):
				mat1[i,i] = 0
				mat2[i,i] = 0
			prob = 0

			n1aprob.shape = 20,1
			n2aprob.shape = 20,1
			prob = np.sum(np.multiply(np.multiply(n1aprob,n2aprob),np.multiply(mat1,mat2)))
			probL.append(prob)
			# print >>OUTSITE,"%f"%(prob)

		print >>OUT,"%.4f\t%.4f"%(sum(cntL),sum(probL))
		print "%.4f\t%.4f"%(sum(cntL),sum(probL))
		pairres[pairlist.index(pairL)][0] += sum(cntL)
		pairres[pairlist.index(pairL)][1] += sum(probL)
		# pairres[pairlist.index(pairL)][2] += t0+t1+t2+t3
		cnt += 1

PAIROUT = open(prefix+"_pairwise_res.txt",'w')
print >>PAIROUT,'\n'.join(['-'.join(pairlist[pairres.index(y)])+'\t'+'\t'.join([str(x) for x in y]) for y in pairres])
