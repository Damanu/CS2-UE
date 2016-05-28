#usr/bin/python!

#Dynamics simulation of katalytic network
#input: initialized network, ntprint, ntperturb number of steps after which the network gets saved and perturbed afterwards, ntjob total number of steps executed

#------------------imports-----------------
import sys
import numpy as np
import networkx as nx

#---------------Globals---------------------------
setsize=3
smallx=0.002# small x for new node and perturbation (smaller than that)
nt=0
#-----------------inputs-----------------------
if len(sys.argv)!=4:
	print "Usage:\n\n	dynamics.py [ntprint] [ntjob] [ntperturb]"
	sys.exit()

ntprint=int(sys.argv[1])
ntjob=int(sys.argv[2])
ntperturb=int(sys.argv[3])

f=open("init_parameter_out.data","r")
data=f.read()
(S,p)=data.split()
S=int(S)
p=float(p)
f.close()

#---------------Subroutines-----------------------
def timestep(G,x,C):
#	print x[list(x).index(min(x))]*sum(np.squeeze(np.asarray(np.dot(C,x))))
	x=x+np.squeeze(np.asarray(np.dot(C,x)))-x*sum(np.squeeze(np.asarray(np.dot(C,x))))#rate equation
	x=x/sum(abs(x))
	return x

def print_net(G,x,nt):
	x_=map(float,x)	#change all elements to a float
	dic=dict(zip(map(str,np.arange(S)),x_))	#write dictionary for the attributes
	nx.set_node_attributes(G,"x",dic) #set attributes new
	fn=str("./dynamic_out/dynamic_out_nt="+str(nt)+".gexf")
	nx.write_gexf(G,fn) #write gexf file with network

def kill_node(G,x):
	weak=[]
	x_=x
	dic=dict(zip(np.arange(S),x_))
	index=0
	while len(weak)<setsize: #till set is full
		weak.append(min(dic.values())) #append the smallest x
		print len(weak)
#		print dic.values().index(min(dic.values()))
		del dic[dic.values().index(min(dic.values()))] #delete the entry of the smallest x in the dict
		print dic
	index=np.random.randint(len(weak))
	j=0
	while j<S:	#iterate over all nodes, set new edges of new species
		if j!=index and np.random.rand()<p: #is an other node and with certain probability?
			G.add_edge(str(index),str(j))	#add a connection between new node and an other one
		elif G.has_edge(j,index):
			G.remove_edge(j,index)	#remove a connection between new node and an other one
		if j!=index and np.random.rand()<p: #is an other node and with certain probability?
			G.add_edge(str(j),str(index))	#add a connection between new node and an other one
		elif G.has_edge(index,j):
			G.remove_edge(index,j)	#remove a connection between new node and an other one
		j+=1
	x_=perturb_all(x_)
	dic=dict(zip(map(str,np.arange(S)),x_))
	nx.set_node_attributes(G,"x",dic)
	return G


def perturb_all(x_):	#perturb all x by an amount smaller than smallx and rescale to normalization
	i=0
	while i<S:
		x_[i]+=smallx*(np.random.rand()-0.5)
		i+=1
	x_=x_/sum(x_)	
	return x_

#-----------------Main program-------------------
def main():
	G=nx.read_gexf("init_out.gexf")
	x=np.array(nx.get_node_attributes(G,"x").values())
	dic=dict(zip(map(str,np.arange(S)),x))
	C=nx.to_numpy_matrix(G)
	nt=0
	while nt<ntjob:
		C=nx.to_numpy_matrix(G)
		x=np.array(nx.get_node_attributes(G,"x").values())
		x=timestep(G,x,C) #change the attributes according to the rate equation
		nt+=1
		if nt%ntprint==0: #every time i has run ntprint steps through
			print_net(G,x,nt)
		if nt%ntperturb==0: #every time i has run ntperturb steps through
			G=kill_node(G,x)
		
	print x			
if __name__=="__main__":
	main()
