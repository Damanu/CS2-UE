#usr/bin/python!

#Dynamics simulation of katalytic network
#input: initialized network, ntprint, ntperturb number of steps after which the network gets saved and perturbed afterwards, ntjob total number of steps executed

#------------------imports-----------------
import sys
import numpy as np
import networkx as nx


#-----------------inputs-----------------------
f=open("init_parameter_out.data","r")
data=f.read()
(S,p)=data.split()
S=int(S)
p=float(p)
f.close()

if len(sys.argv)!=4:
	print "Usage:\n\n	dynamics.py [ntprint] [ntjob] [ntperturb]"
	sys.exit()

ntprint=int(sys.argv[1])
ntjob=int(sys.argv[2])
ntperturb=int(sys.argv[3])


#---------------Globals---------------------------
setsize=int(S*0.5)
smallx=0.1/S# small x for new node and perturbation (smaller than that)
nt=0
#---------------Subroutines-----------------------
def timestep(G,x,C,nt):
	x_=x
	eig=np.linalg.eig(C)
#	x_=x+(np.squeeze(np.asarray(np.dot(C,x)))-x*sum(np.squeeze(np.asarray(np.dot(C,x)))))#rate equation
	maxew=max(np.real(eig[0]))
	if maxew == 0.:
		print "000000"
		x_=x+(np.squeeze(np.asarray(np.dot(C,x)))-x*sum(np.squeeze(np.asarray(np.dot(C,x)))))#rate equation
	else:
		maxind=list(eig[0]).index(maxew)
		x=np.squeeze(np.asarray(np.real(eig[1].transpose()[maxind])))
		print eig[0]
		print eig[1].transpose()[maxind]

		print x		
#	print maxew
#	print eig[0]
#	print eig[1].transpose()[maxind]
#	x=eig[1].transpose()[maxind]
#	print x		
#	if any(n==maxew for n in list(np.real(eig[0]))):
#		for val2 in eig[1].transpose():
#		#	print np.squeeze(np.asarray(np.real(val2)))
#			if any(n2 < 0 for n2 in np.squeeze(np.asarray(np.real(val2)))) or sum(np.squeeze(np.asarray(np.real(val2))))<=0:
#				print "break"
#			else:
#				x=np.squeeze(np.asarray(np.real(val2)))		
#				#print x
#				break;
					
#	i=0
#	while i<S:
#		j=0
#		while j<S:
#			x[i]+=C[i,j]*x[j]
#			j+=1
#		k=0
#		sum2=0
#		while k<S:
#			j=0
#			while j<S:
#				sum2+=C[k,j]*x[j]
#				j+=1
#			k+=1
#		x[i]-=x[i]*sum2
#		i+=1

	if ntperturb==0:
		print nt,sum(np.squeeze(np.asarray(np.dot(C,x)))-x*sum(np.squeeze(np.asarray(np.dot(C,x)))))
	i=0
	while i<S:
		if x[i]<0:
			x[i]=0
		i+=1
	x=x/sum(x)
	if min(x)<0:
		print "x<0 error"
		sys.exit()
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
		keys=dic.keys()
		np.random.shuffle(keys)
		l=[]	
		for n in keys:
			l.append(dic[n])
		ind=l.index(min(map(float,l)))
		weak.append(int(keys[ind])) #append the smallest x key
#		print dic.values().index(min(dic.values()))
		del dic[keys[ind]] #delete the entry of the smallest x in the dict
	dic=dict(zip(np.arange(S),x_))
	index=str(np.random.choice(weak))
	x_[int(index)]=smallx #change x of new species
#	print weak
#	print "kill: ",index,dic[int(index)]

	j=0
	while j<S:
		if G.has_edge(index,str(j)):
			G.remove_edge(index,str(j))	#remove a connection between new node and an other one
		if G.has_edge(str(j),index):
			G.remove_edge(str(j),index)	#remove a connection between new node and an other one
		j+=1
	j=0	
	while j<S:	#iterate over all nodes, set new edges of new species
		if j!=int(index) and np.random.rand()<p: #is an other node and with certain probability?
			G.add_edge(index,str(j))	#add a connection between new node and an other one
		if j!=int(index) and np.random.rand()<p: #is an other node and with certain probability?
			G.add_edge(str(j),index)	#add a connection between new node and an other one
		j+=1
	x_=perturb_all(x_)
	dic=dict(zip(map(str,np.arange(S)),x_))
	nx.set_node_attributes(G,"x",dic)
	return x,G
def nonzerox(x):
	count=0
	for val in x:
		if val > 10**-6:
			count+=1
	return count

def perturb_all(x_):	#perturb all x by an amount smaller than smallx and rescale to normalization
	i=0
	while i<S:
		x_[i]-=smallx*(np.random.rand())/1.
		if x_[i]<0:
			x_[i]=0
		i+=1
	x_=x_/sum(x_)	
	return x_

#-----------------Main program-------------------
def main():
	G=nx.read_gexf("init_out.gexf")
	x=np.array(nx.get_node_attributes(G,"x").values())
	dic=dict(zip(map(str,np.arange(S)),x))
	nodel=map(str,np.arange(S))
	C=nx.to_numpy_matrix(G,nodelist=nodel).transpose()
	nt=0
	while nt<ntjob:
#		x=np.array(nx.get_node_attributes(G,"x").values())
		x=timestep(G,x,C,nt) #change the attributes according to the rate equation
		nt+=1
		if nt%ntprint==0: #every time i has run ntprint steps through
			print_net(G,x,nt/ntprint)
			print nt/ntprint,nonzerox(x),np.real(max(np.linalg.eig(C)[0]))
#			print "loop: ",x[95],x[97],x[98]
		#	print "x ",x
		if ntperturb!=0: 
			if nt%ntperturb==0: #every time i has run ntperturb steps through
				x,G=kill_node(G,x)
				C=nx.to_numpy_matrix(G,nodelist=nodel).transpose()
#				x=np.array(nx.get_node_attributes(G,"x").values())
		
if __name__=="__main__":
	main()
