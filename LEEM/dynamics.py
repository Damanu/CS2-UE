#usr/bin/python!

#Dynamics simulation of katalytic network
#input: initialized network, ntprint, ntperturb number of steps after which the network gets saved and perturbed afterwards, ntjob total number of steps executed

#------------------imports-----------------
import sys
import numpy as np
import networkx as nx

#---------------Globals---------------------------
tres=0.01
smallx=0.002

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
	x=x+np.squeeze(np.asarray(np.dot(C,x)))-x*sum(np.squeeze(np.asarray(np.dot(C,x))))
	x=x/sum(x)
	return x
#-----------------Main program-------------------
def main():
	weak=[]

	G=nx.read_gexf("init_out.gexf")
	x=np.array(nx.get_node_attributes(G,"x").values())
	C=nx.to_numpy_matrix(G)
	i=0
	while i<ntjob:
		x=timestep(G,x,C)
		i+=1
		print "i ",i
		if i%ntprint==0: #every time i has run ntprint steps through
			x_=map(float,x)	#change all elements to a float
			dic=dict(zip(map(str,np.arange(S)),x_))	#write dictionary for the attributes
			nx.set_node_attributes(G,"x",dic) #set attributes new
			fn=str("./dynamic_out/dynamic_out_nt="+str(i)+".gexf")
			nx.write_gexf(G,fn) #write gexf file with network
		if i%ntperturb==0: #every time i has run ntperturb steps through
			x_=map(float,x) 
			dic=dict(zip(map(str,np.arange(S)),x_))
			nx.set_node_attributes(G,"x",dic)
			weak=[]
			for index,xi in dic.iteritems():
				if xi < tres:
					weak.append(index)
				#	print index
			delind=np.random.randint(len(weak))
			print delind
#			print G.node
			G.remove_node(str(delind)) #remove one random node of the weak nodes
			G.add_node(str(delind))	#add node with same name
			dic[str(delind)]=smallx #set the relative concentration of the new node to a small value
			j=0
			while j<S:	#iterate over all nodes
				if j!=delind and np.random.rand()<p: #is an other node and with certain probability?
					G.add_edge(delind,j)	#add a connection between new node and an other one
				if j!=delind and np.random.rand()<p: #is an other node and with certain probability?
					G.add_edge(delind,j)	#add a connection between new node and an other one
				j+=1

				
			
					
if __name__=="__main__":
	main()
