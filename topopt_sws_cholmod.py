# A 200 LINE TOPOLOGY OPTIMIZATION CODE BY NIELS AAGE AND VILLADS EGEDE JOHANSEN, JANUARY 2013
# Updated by Niels Aage February 2016
from __future__ import division
import numpy as np

import osqp
from scipy.sparse import csc_matrix
from scipy.optimize import linprog

from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt
import cvxopt ;import cvxopt.cholmod

def main(nelx,nely,volfrac,penal,rmin,ft):
	print("Minimum compliance self-weight problem with LP(Scipy)+AML")
	print("ndes: " + str(nelx) + " x " + str(nely))
	print("volfrac: " + str(volfrac) + ", rmin: " + str(rmin) + ", penal: " + str(penal))
	print("Filter method: " + ["Sensitivity based","Density based"][ft])

	mov_fct=0.5e0*np.ones(nelx*nely)

	# Max and min stiffness
	Emin=1e-9
	Emax=1.0

	# dofs:
	ndof = 2*(nelx+1)*(nely+1)

	# Allocate design variables (as array), initialize and allocate sens.
	x=volfrac * np.ones(nely*nelx,dtype=float)
	xold=x.copy()
	xolder=x.copy()
	xoldest=x.copy()
	xPhys=x.copy()

	g=0 # must be initialized to use the NGuyen/Paulino OC approach (change)
	dc=np.zeros(nely*nelx, dtype=float)

	# FE: Build the index vectors for the for coo matrix format.
	KE=lk()
	edofMat=np.zeros((nelx*nely,8),dtype=int)
	for elx in range(nelx):
		for ely in range(nely):
			el = ely+elx*nely
			n1=(nely+1)*elx+ely
			n2=(nely+1)*(elx+1)+ely
			edofMat[el,:]=np.array([2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3,2*n2, 2*n2+1, 2*n1, 2*n1+1])
	# Construct the index pointers for the coo format
	iK = np.kron(edofMat,np.ones((8,1))).flatten()
	jK = np.kron(edofMat,np.ones((1,8))).flatten()    

	# Filter: Build (and assemble) the index+data vectors for the coo matrix format
	nfilter=int(nelx*nely*((2*(np.ceil(rmin)-1)+1)**2))
	iH = np.zeros(nfilter)
	jH = np.zeros(nfilter)
	sH = np.zeros(nfilter)
	cc=0
	for i in range(nelx):
		for j in range(nely):
			row=i*nely+j
			kk1=int(np.maximum(i-(np.ceil(rmin)-1),0))
			kk2=int(np.minimum(i+np.ceil(rmin),nelx))
			ll1=int(np.maximum(j-(np.ceil(rmin)-1),0))
			ll2=int(np.minimum(j+np.ceil(rmin),nely))
			for k in range(kk1,kk2):
				for l in range(ll1,ll2):
					col=k*nely+l
					fac=rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))
					iH[cc]=row
					jH[cc]=col
					sH[cc]=np.maximum(0.0,fac)
					cc=cc+1
	# Finalize assembly and convert to csc format
	H=coo_matrix((sH,(iH,jH)),shape=(nelx*nely,nelx*nely)).tocsc()	
	Hs=H.sum(1)

	# BC's and support (change)
	dofs=np.arange(2*(nelx+1)*(nely+1))
	ndofy=2*(nely+1)
	fixed = np.union1d(dofs[0:ndofy:2],
		np.array([ndof - 1]))
	free=np.setdiff1d(dofs,fixed)

	# Solution and RHS vectors
	f=np.zeros((ndof,1))
	u=np.zeros((ndof,1))

	# Set load (move)

	# Initialize plot and plot the initial design
	plt.ion() # Ensure that redrawing is possible
	fig,ax = plt.subplots()
	im = ax.imshow(-xPhys.reshape((nelx,nely)).T, cmap='gray',\
	interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
	fig.show()

	loop=0
	change=1
	dv = np.ones(nely*nelx)
	dc = np.ones(nely*nelx)
	ce = np.ones(nely*nelx)
	while change>0.1 and loop<2000:
		loop=loop+1
		# Setup and solve FE problem
		sK=((KE.flatten()[np.newaxis]).T*(Emin+(xPhys)**penal*(Emax-Emin))).flatten(order='F')
		K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
		# Remove constrained dofs from matrix and convert to coo
		K = deleterowcol(K,fixed,fixed).tocoo()
 
	    # Set load (change)
		f_tmp=np.zeros(ndof)
		np.add.at(f_tmp, edofMat[:, 1::2].flatten(),
			np.kron(xPhys, -1e0/5000 * np.ones(4) / 4))
		f[:,0] = f_tmp
 
		# Solve system 
		K = cvxopt.spmatrix(K.data,K.row.astype(np.int),K.col.astype(np.int))
		B = cvxopt.matrix(f[free,0])
		cvxopt.cholmod.linsolve(K,B)
		u[free,0]=np.array(B)[:,0] 

		# Objective and sensitivity (change)
		ce[:] = (np.dot(u[edofMat].reshape(nelx*nely,8),KE)*u[edofMat].reshape(nelx*nely,8) ).sum(1)
		obj=( (Emin+xPhys**penal*(Emax-Emin))*ce ).sum()
		dc[:]=(-penal*xPhys**(penal-1)*(Emax-Emin))*ce
		dc[:] -= u[edofMat[:, 1], 0] * 1e0 / 10000; dc[:] -= u[edofMat[:, 3], 0] * 1e0 / 10000
		dc[:] -= u[edofMat[:, 5], 0] * 1e0 / 10000; dc[:] -= u[edofMat[:, 7], 0] * 1e0 / 10000
		dv[:] = -np.ones(nely*nelx)/float(nelx*nely)/volfrac
		g=1e0-np.sum(x)/float(nelx*nely)/volfrac

		# Sensitivity filtering:
		if ft==0:
			dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
		elif ft==1:
			dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0]
			dv[:] = np.asarray(H*(dv[np.newaxis].T/Hs))[:,0]

#	dc[:]=dc/6.2e7
#	obj=obj/6.2e7

		# Optimality criteria (change)
		xold[:]=x
		(x[:])=lp(nelx,nely,x,volfrac,dc,dv,g,xolder,xoldest,loop,mov_fct)
		xoldest[:]=xolder; xolder[:]=xold
 
		# Filter design variables
		if ft==0:   xPhys[:]=x
		elif ft==1:	xPhys[:]=np.asarray(H*x[np.newaxis].T/Hs)[:,0]
	    
		# Compute the change by the inf. norm
		change=np.linalg.norm(x.reshape(nelx*nely,1)-xold.reshape(nelx*nely,1),np.inf)
		bwf=(sum(y>0.99 for y in x)+sum(y<0.01 for y in x))/nelx/nely

		# Plot to screen
		im.set_array(-xPhys.reshape((nelx,nely)).T)
		fig.canvas.draw()

		# Write iteration history to screen (req. Python 2.6 or newer)
		print("it.: {0} , obj.: {1:.3f} Vol.: {2:.3f}, ch.: {3:.3f}, bw: {4:.3f}".format(\
					loop,obj,g,change,bwf))

	# Make sure the plot stays and that the shell remains	
	plt.show()
	input("Press any key...")
    
#element stiffness matrix
def lk():
	E=1
	nu=0.3
	k=np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
	KE = E/(1-nu**2)*np.array([ [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
	[k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
	[k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
	[k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
	[k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
	[k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
	[k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
	[k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]] ]);
	return (KE)

def lp(nelx,nely,x,volfrac,dc,dv,g,xolder,xoldest,loop,mov_fct):
#
	n=nelx*nely
#
	mov=0.2
#
	if loop > 2:
		osc=(x-xolder)*(xolder-xoldest)/mov
		for i in range(n):
			if osc[i] > 1e-9: mov_fct[i]=mov_fct[i]*1.2
			else: mov_fct[i]=mov_fct[i]*0.7
#
	dx_l=np.maximum(x-mov_fct*mov,np.zeros(n))-x
	dx_u=np.minimum(x+mov_fct*mov,np.ones(n))-x
#
	bds=[]
	for i in range(n):
		bds.append((dx_l[i],dx_u[i]))
#
	res = linprog(dc, A_ub=dv.reshape((1,n)), b_ub=-g, bounds=bds,method='interior-point')
#
	xnew=x+np.maximum(np.minimum(res.x,dx_u),dx_l)
#
	return (xnew)
#
def oc(nelx,nely,x,volfrac,dc,dv,g):
	l1=0
	l2=1e9
	move=0.2
	# reshape to perform vector operations
	xnew=np.zeros(nelx*nely)

	while (l2-l1)/(l1+l2)>1e-3:
		lmid=0.5*(l2+l1)
		xnew[:]= np.maximum(0.0,np.maximum(x-move,np.minimum(1.0,np.minimum(x+move,x*np.sqrt(
			np.maximum(-dc/dv/lmid,1e-3*np.ones(nelx*nely))    )))))
		gt=g+np.sum((dv*(xnew-x)))
		if gt<0 :
			l1=lmid
		else:
			l2=lmid
	return (xnew,gt)
    
def deleterowcol(A, delrow, delcol):
	# Assumes that matrix is in symmetric csc form !
	m = A.shape[0]
	keep = np.delete (np.arange(0, m), delrow)
	A = A[keep, :]
	keep = np.delete (np.arange(0, m), delcol)
	A = A[:, keep]
	return A    

if __name__ == "__main__":
	# Default input parameters
	nelx=180
	nely=60
	volfrac=0.4
	rmin=5.4
	penal=3.0
	ft=1 # ft==0 -> sens, ft==1 -> dens

	import sys
	if len(sys.argv)>1: nelx   =int(sys.argv[1])
	if len(sys.argv)>2: nely   =int(sys.argv[2])
	if len(sys.argv)>3: volfrac=float(sys.argv[3])
	if len(sys.argv)>4: rmin   =float(sys.argv[4])
	if len(sys.argv)>5: penal  =float(sys.argv[5])
	if len(sys.argv)>6: ft     =int(sys.argv[6])

	main(nelx,nely,volfrac,penal,rmin,ft)
