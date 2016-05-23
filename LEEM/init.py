#/usr/bin/python!
#Initailisattion of a directed random erdoes-renyi Network with probability p to connect

#Input Parameter: Number of Nodes S, connecting probability p

#------------Imports---------------------------
import sys
import networkx as nx
import numpy as np
#------------Globals--------------------------
OUTPATH="./init_out.net"
#------------Input parameters------------------
if len(sys.argv)!=3:
	print "Help\nCall the program like that:\n\n	init.py [number of nodes] [connecting probability p]\n"
	sys.exit()

S=int(sys.argv[1])
p=float(sys.argv[2])

#------------Main Program----------------------
def main():
	G=nx.fast_gnp_random_graph(S,p,directed=True) #make a random erdoes-renyi graph with S nodes and prob p
	ran=np.random.rand(S)	# make list with S random numbers between 0 and 1
	ran=ran/sum(ran)	# normalise ran so that the sum over all is 1
	ran=map(float,ran)
	x=dict(zip(np.arange(S),ran)) #create a dictionary with first parameter the index of the node and second the relative population xi
	
	nx.set_node_attributes(G,"x",x)	#set the attribute of the nodes
	nx.write_gexf(G,"init_out.gexf") #write network to file
	f=open("init_parameter_out.data","w")
	f.write(str(S))
	f.write(" ")
	f.write(str(p))
	
	f.close()
if __name__=="__main__":
	main() 
