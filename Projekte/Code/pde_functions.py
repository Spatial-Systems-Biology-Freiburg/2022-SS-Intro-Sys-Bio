from scipy import sparse
import numpy as np
import time

def jpat(D,ind):
	NVar = ind[1]-ind[0]
	sysmat = np.ones((NVar, NVar))
	diffmat = sparse.eye(NVar)
	Js = sparse.kron(sparse.eye(len(D)), sysmat)
	sD = sparse.csr_matrix(D)
	sD.data.fill(1)
	Jd = sparse.kron(sD, diffmat)
	J = Js + Jd
	J.data.fill(1)
	return J

def sd(t,y,D,ind,k):
	dydt = np.zeros(y.shape)
	TTG1 = y[ind]
	GL3  = y[ind+1]
	AC   = y[ind+2]
	dydt[ind]   = k[0]-k[1]*TTG1 - TTG1*GL3 + k[2]*np.dot(D,TTG1)
	dydt[ind+1] = k[3]*(AC*AC) - GL3 - TTG1*GL3
	dydt[ind+2] = TTG1*GL3 - AC
	return dydt


def full_model(t,y,D,ind,k,start_time=None):
	if start_time != None:
		print("[{: >8.4f}s] Solving ...".format(time.time()-start_time), end="\r")
	dydt = np.zeros(y.shape)
	TTG1 = y[ind]
	GL1  = y[ind+1]
	GL3  = y[ind+2]
	TRY  = y[ind+3]
	CPC  = y[ind+4]
	AC1  = y[ind+5]
	AC2  = y[ind+6]
	dydt[ind]   = k[0] - TTG1*(k[1] + k[2]*GL3) + (k[1]*k[3])*np.dot(D,TTG1)
	dydt[ind+1] = k[4] + k[5]*AC2 - GL1*(k[6] + k[7]*GL3)
	dydt[ind+2] = k[8] + (k[22]*k[10]*AC2*AC2)/(k[22]+AC2*AC2) \
                - GL3*(k[11] + k[7]*GL1 + k[13]*CPC) \
                + (k[23]*k[9]*AC1*AC1)/(k[23]+AC1*AC1) \
                - k[2]*TTG1*GL3 - k[12]*TRY*GL3
	dydt[ind+3] = k[14]*AC1*AC1 - TRY*k[15] - TRY*GL3*k[12] + k[15]*k[16]*np.dot(D,TRY)
	dydt[ind+4] = k[17]*AC2*AC2 - k[18]*CPC - k[13]*CPC*GL3 + k[18]*k[19]*np.dot(D,CPC)
	dydt[ind+5] = k[2]*GL3*TTG1 - k[20]*AC1
	dydt[ind+6] = k[7]*GL3*GL1 - k[21]*AC2
	return dydt


def jac_full_model(t,y,k):
	dy_dtdp = np.zeros((y.size, y.size))
	TTG1 = y[0]
	GL1  = y[1]
	GL3  = y[2]
	TRY  = y[3]
	CPC  = y[4]
	AC1  = y[5]
	AC2  = y[6]
	# TODO there is a error somewhere in here!
	dy_dtdp = np.array([
		[-(k[1]+k[2]*GL3), 0, -TTG1*k[2], 0, 0, 0, 0],
		[0, -(k[6]+k[7]*GL3), -GL1*k[7], 0, 0, 0, k[5]],
		[
			-k[2]*GL3,
			-k[7]*GL3, 
			-(k[11]+k[7]*GL1+k[13]*CPC)-k[2]*TTG1-k[12]*TRY,
			-k[12]*GL3,
			-k[13]*GL3,
			2*(k[23]*k[ 9]*AC1)/(k[23]+AC1*AC1) - 2*AC1*(k[23]*k[ 9]*AC1*AC1)/(k[23]+AC1*AC1)**2,
			2*(k[22]*k[10]*AC2)/(k[22]+AC2*AC2) - 2*AC2*(k[22]*k[10]*AC2*AC2)/(k[22]+AC2*AC2)**2
		],
		[0, 0, TRY*k[12], GL3*k[12], 0, 2*k[14]*AC1, 0],
		[0, 0, -k[13]*CPC, 0, -k[13]*GL3, 0, 2*k[17]*AC2],
		[k[2]*GL3, 0, k[2]*TTG1, 0, 0, -k[20], 0],
		[0, k[7]*GL3, k[7]*GL1, 0, 0, 0, -k[21]]
	])
	return dy_dtdp


def MYC1_model(t,y,D,ind,k,start_time=None):
	if start_time != None:
		print("[{: >8.4f}s] Solving ...".format(time.time()-start_time), end="\r")
	dydt = np.zeros(y.shape)
	Ac = y[ind]
	An = y[ind+1]
	Ic = y[ind+2]
	In = y[ind+3]
	T  = y[ind+4]

	dydt[ind]   = An*(k[0] + k[1]*T) - Ac*(k[0] + k[2]) + np.dot(D, Ac)
	dydt[ind+1] = k[3] + (k[4]*An*An)/(k[5] + In) - An*(k[0] + k[2] + k[1]*T) + k[0]*Ac
	dydt[ind+2] = k[0]*In - Ic*(k[0] + k[6] + k[7]*T) + k[8]*np.dot(D,Ic)
	dydt[ind+3] = k[9]*An*An + Ic*(k[0] + k[7]*T) - In*(k[6] + k[0])
	dydt[ind+4] = k[10] - T*(k[11] + k[7]*Ic + k[1]*An)
	return dydt


def jac_MYC1_model(t,y,k):
	Ac = y[0]
	An = y[1]
	Ic = y[2]
	In = y[3]
	T  = y[4]

	dy_dtdp = np.array([
		[-k[0]-k[2],k[0] + T*k[1],0,0,An*k[1]],
		[k[0], (2*An*k[4])/(In + k[5])-k[2]-T*k[1]-k[0],0,-(An**2*k[4])/(In + k[5])**2, -An*k[1]],
		[0,0,-k[0]-k[6]-T*k[7],k[0],-Ic*k[7]],
		[0,2*An*k[9],k[0] + T*k[7],-k[0]-k[6],Ic*k[7]],
		[0,-T*k[1],-T*k[7],0, -k[11]-An*k[1]-Ic*k[7]]
	])
	return dy_dtdp