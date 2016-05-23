#usr/bin/python!

#Dynamics simulation of katalytic network
#input: initialized network, ntprint number of steps after which the network gets saved and perturbed afterwards, ntjob total number of steps executed

#------------------imports-----------------
import sys
import numpy as np
import networkx as nx

#---------------Globals---------------------------



#-----------------inputs-----------------------
if len(sys.argv)!=3:
	print "Usage:\n\n	dynamics.py [ntprint] [ntjob]"
	sys.exit()

ntprint=int(sys.argv[1])
ntjob=int(sys.argv[2])
f=open("init_parameter_out.data","r")
data=f.read()
(S,p)=data.split()
f.close()

#---------------Subroutines-----------------------
def timestep(G,x,C):
	x=x+np.squeeze(np.asarray(np.dot(C,x)))-x*sum(np.squeeze(np.asarray(np.dot(C,x))))
	x=x/sum(x)
	return x
#-----------------Main program-------------------
def main():
	G=nx.read_gexf("init_out.gexf")
	x=np.array(nx.get_node_attributes(G,"x").values())
	C=nx.to_numpy_matrix(G)
	print x
	print len(C)
	i=0
	while i<ntjob:
		x=timestep(G,x,C)
		if i%ntprint==0:
			print x
			print sum(x)
			map(float,x)
			dic=dict(zip(np.arange(int(S)),x))
			nx.set_node_attributes(G,"x",dic)			
			fn=str("./dynamic_out/dynamic_out_nt="+str(i)+".gexf")
			nx.write_gexf(fn)
		i+=1
if __name__=="__main__":
	main()
