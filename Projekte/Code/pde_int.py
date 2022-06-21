import numpy as np

def couplingMatrix(xmax=1,ymax=1,bndcondition="zeroflux",celltype="quadratic"):
	if ((xmax==1) or (ymax==1)):
		ysten=[0,0]
		xsten=[-1,1]
	elif (xmax>0) and (ymax>0):
		if (celltype=="quadratic"):
			ysten=[-1,1,0,0]
			xsten=[0,0,-1,1]
		elif (celltype=="hexagonal"):
			ysten=[-1,1,0,0,-1,1]
			xsten=[0,0,-1,1,1,-1]
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
	return __cmatrix(ymax,xmax,ysten,xsten,bflag,1)


def __cmatrix(ymax,xmax,ysten,xsten,bflag,cflag):
	n=ymax*xmax
	D=np.zeros((n,n))
	idx=np.reshape(np.arange(n),(ymax,xmax), order='F')
	for s in np.arange(len(ysten)):
		sidx=np.roll(np.roll(idx,xsten[s],axis=1),ysten[s],axis=0)
		ind=np.ravel_multi_index((idx,sidx),(n,n))
		if bflag:
			if (ysten[s]>0):
				ind=np.delete(ind,0,0)
			elif (ysten[s]<0):
				ind=np.delete(ind,-1,0)
			if (xsten[s]>0):
				ind=np.delete(ind,0,1)
			elif (xsten[s]<0):
				ind=np.delete(ind,-1,1)
		ind=np.unravel_index(ind,(n,n))
		D[ind]=1
	if cflag:
		D=D-np.diag(np.sum(D,axis=1))
	return D.astype(int)
	
def IJKth(s,y,x,ymax,NVar):
	xp=np.tile(x*ymax*NVar,(len(y),1))
	yp=np.tile(y*NVar,(len(x),1)).T
	idx=xp+yp+s-1
	return idx.flatten('F').astype(int)