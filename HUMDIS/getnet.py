#usr/bin/python!

import numpy as np
import networkx as nx
import itertools as it
import sys

FILE="data.dat"

def main():
	data=[]
	f=open(FILE,"r")
	i=0
	for line in f:
		data.append(line[0:-1].split(";"))
		i+=1
	
	humdis=nx.Graph()
	humdis2=nx.Graph()
	humgene=nx.Graph()
	humbipart=nx.Graph()
#human disease network
	for line in data[1:]:	
		genes_=line[2].replace(" ","").split(",") #get a list of genes in the new line
		if humdis.has_node(line[1]): #ask if the disease already exists
			for x in genes_:
				if x not in nx.get_node_attributes(humdis,"genes")[line[1]]:#ask if one of the genes  is not in the gene pool of the disease
					nx.get_node_attributes(humdis,"genes")[line[1]].append(x)
		else:
			humdis.add_node(line[1],nid=int(line[0]),genes=genes_,Class=line[5])
		
		if humdis2.has_node(line[0]): #ask if the disease already exists
			for x in genes_:
				if x not in nx.get_node_attributes(humdis2,"genes")[line[0]]:#ask if one of the genes  is not in the gene pool of the disease
					nx.get_node_attributes(humdis2,"genes")[line[0]].append(x)
			nx.set_node_attributes(humdis2,"numofgenes",{line[0]:float(len(nx.get_node_attributes(humdis2,"genes")[line[0]]))})

		else:
			humdis2.add_node(line[0],nid=int(line[0]),genes=genes_,Class=line[5],numofgenes=int(len(genes_)))
	#make edges

	#make edges
	nodebuff=[]
#human gene network
	for line in data[1:]:	
		genes_=line[2].replace(" ","").split(",") #get a list of genes in the new line
		for x in genes_:
			if humgene.has_node(x):
				nx.get_node_attributes(humgene,"disease")[x].append(line[1])
				nx.get_node_attributes(humgene,"nid")[x].append(line[0])
				if line[5] != nx.get_node_attributes(humgene,"Class")[x]:
					nx.set_node_attributes(humgene,"Class",{x:0})
					
			else:
				humgene.add_node(x,num=float(line[0]),numofdis=float(1),disease=[line[1]],nid=[line[0]],Class=line[5])
			nx.set_node_attributes(humgene,"numofdis",{x:float(len(nx.get_node_attributes(humgene,"disease")[x]))})
	
#biparted network
	humbipart.add_nodes_from(humdis2.nodes(),bipartite=0)
	humbipart.add_nodes_from(humgene.nodes(),bipartite=1)

	for node in humgene.nodes(): #got through all genes
		diseases=nx.get_node_attributes(humgene,"disease")[node]
		edges=it.combinations(diseases,2)
		humdis.add_edges_from(edges)
	for node in humgene.nodes(): #got through all genes
		nids=nx.get_node_attributes(humgene,"nid")[node]
		edges=it.combinations(nids,2)
		humdis2.add_edges_from(edges)
	for node in humdis2.nodes(): #got through all genes
		genes=nx.get_node_attributes(humdis2,"genes")[node]
		edges=it.combinations(genes,2)
		humgene.add_edges_from(edges)
	
	nx.set_node_attributes(humdis,"genes",0)
	nx.set_node_attributes(humdis2,"genes",0)
	nx.set_node_attributes(humgene,"nid",0)
	nx.set_node_attributes(humgene,"disease",0)
	nx.write_gexf(humdis,"text.gexf")	
	nx.write_gexf(humdis2,"text2.gexf")	
	nx.write_gexf(humgene,"text3.gexf")	
	
	
if __name__=="__main__":
	main();
