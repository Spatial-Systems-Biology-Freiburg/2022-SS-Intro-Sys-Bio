from scipy import *
import numpy as np

def couplingMatrix(xmax=1,ymax=1,bndcondition="zeroflux",celltype="quadratic"):
	if ((xmax==1) or (ymax==1)):
		ysten=r_[0,0]
		xsten=r_[-1,1]
	elif (xmax>0) and (ymax>0):
		if (celltype=="quadratic"):
			ysten=r_[-1,1,0,0]
			xsten=r_[0,0,-1,1]
		elif (celltype=="hexagonal"):
			ysten=r_[-1,1,0,0,-1,1]
			xsten=r_[0,0,-1,1,1,-1]
		else:
			print("Unknown cell type!")
	else:
		print("Impossible dimensions!")
	if (bndcondition=="zeroflux"):
		bflag=1
	elif (bndcondition=="periodic"):
		bflag=0
	else:
		print("Unknown boundary condition!")
	return cmatrix(ymax,xmax,ysten,xsten,bflag,1)
	
def cmatrix(ymax,xmax,ysten,xsten,bflag,cflag):
	n=ymax*xmax
	D=zeros((n,n))
	idx=reshape(np.arange(n),(ymax,xmax), order='F')
	for s in np.arange(len(ysten)):
		sidx=roll(roll(idx,xsten[s],axis=1),ysten[s],axis=0)
		ind=ravel_multi_index((idx,sidx),(n,n))
		if bflag:
			if (ysten[s]>0):
				ind=delete(ind,0,0)
			elif (ysten[s]<0):
				ind=delete(ind,-1,0)
			if (xsten[s]>0):
				ind=delete(ind,0,1)
			elif (xsten[s]<0):
				ind=delete(ind,-1,1)
		ind=unravel_index(ind,(n,n))
		D[ind]=1
	if cflag:
		D=D-diag(sum(D,axis=1))
	return D.astype(int)
	
def IJKth(s,y,x,ymax,NVar):
	xp=tile(x*ymax*NVar,(len(y),1))
	yp=tile(y*NVar,(len(x),1)).T
	idx=xp+yp+s-1
	return idx.flatten('F').astype(int)