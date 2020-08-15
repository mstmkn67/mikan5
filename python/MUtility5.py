#MUtility.py

from numpy import zeros,array,transpose,dot,trace,abs
from numpy.linalg import eig,det
from UDFManager import *
from copy import copy

def getSubTensorA(udf,index):
	tA=udf.get("particle_type[%d].resistance_tensor.A"%(index))
	A=zeros((3,3),float)
	A[0][0],A[1][1],A[2][2]=tA[0],tA[1],tA[2]
	A[1][0]=A[0][1]=tA[3]
	A[2][0]=A[0][2]=tA[4]
	A[2][1]=A[1][2]=tA[5]
	return A

def getSubTensorB(udf,index):
	tB=udf.get("particle_type[%d].resistance_tensor.B"%(index))
	B=zeros((3,3),float)
	B[0][0],B[1][1],B[2][2]=tB[0],tB[1],tB[2]
	B[1][0]=B[0][1]=tB[3]
	B[2][0]=B[0][2]=tB[4]
	B[2][1]=B[1][2]=tB[5]
	return B

def getSubTensorC(udf,index):
	tC=udf.get("particle_type[%d].resistance_tensor.C"%(index))
	C=zeros((3,3),float)
	C[0][0],C[1][1],C[2][2]=tC[0],tC[1],tC[2]
	C[1][0]=C[0][1]=tC[3]
	C[2][0]=C[0][2]=tC[4]
	C[2][1]=C[1][2]=tC[5]
	return C

def getSubTensorG(udf,index):
	tG=zeros((3,3,3),float)
	G=array(udf.get("particle_type[%d].resistance_tensor.G"%(index)))
	tG[0][0][0],tG[0][0][1],tG[0][0][2]=G[0][0],G[0][1],G[0][2]#G.xx
	tG[1][1][0],tG[1][1][1],tG[1][1][2]=G[1][0],G[1][1],G[1][2]#G.yy
	tG[2][2][0],tG[2][2][1],tG[2][2][2]=G[2][0],G[2][1],G[2][2]#G.zz
	tG[1][0][0]=tG[0][1][0]=G[3][0]
	tG[1][0][1]=tG[0][1][1]=G[3][1]
	tG[1][0][2]=tG[0][1][2]=G[3][2] #G.yx
	tG[2][0][0]=tG[0][2][0]=G[4][0]
	tG[2][0][1]=tG[0][2][1]=G[4][1]
	tG[2][0][2]=tG[0][2][2]=G[4][2] #G.zx
	tG[2][1][0]=tG[1][2][0]=G[5][0]
	tG[2][1][1]=tG[1][2][1]=G[5][1]
	tG[2][1][2]=tG[1][2][2]=G[5][2] #G.zy
	return tG

def getSubTensorH(udf,index):
	tH=zeros((3,3,3),float)
	H=array(udf.get("particle_type[%d].resistance_tensor.H"%(index)))
	tH[0][0][0],tH[0][0][1],tH[0][0][2]=H[0][0],H[0][1],H[0][2]#H.xx
	tH[1][1][0],tH[1][1][1],tH[1][1][2]=H[1][0],H[1][1],H[1][2]#H.yy
	tH[2][2][0],tH[2][2][1],tH[2][2][2]=H[2][0],H[2][1],H[2][2]#H.zz
	tH[1][0][0]=tH[0][1][0]=H[3][0]
	tH[1][0][1]=tH[0][1][1]=H[3][1]
	tH[1][0][2]=tH[0][1][2]=H[3][2] #H.yx
	tH[2][0][0]=tH[0][2][0]=H[4][0]
	tH[2][0][1]=tH[0][2][1]=H[4][1]
	tH[2][0][2]=tH[0][2][2]=H[4][2] #H.zx
	tH[2][1][0]=tH[1][2][0]=H[5][0]
	tH[2][1][1]=tH[1][2][1]=H[5][1]
	tH[2][1][2]=tH[1][2][2]=H[5][2] #H.zy
	return tH

def getSubTensorK(udf,index):
	K=zeros((3,3,3,3),float)
	tK=udf.get("particle_type[%d].resistance_tensor.K.xx"%(index))
	K[0][0][0][0]=tK[0]
	tK=udf.get("particle_type[%d].resistance_tensor.K.yy"%(index))
	K[1][1][0][0]=K[0][0][1][1]=tK[0]
	K[1][1][1][1]=tK[1]
	tK=udf.get("particle_type[%d].resistance_tensor.K.zz"%(index))
	K[2][2][0][0]=K[0][0][2][2]=tK[0]
	K[2][2][1][1]=K[1][1][2][2]=tK[1]
	K[2][2][2][2]=tK[2]
	tK=udf.get("particle_type[%d].resistance_tensor.K.yx"%(index))
	K[1][0][0][0]=K[0][0][1][0]=K[0][1][0][0]=K[0][0][0][1]=tK[0]
	K[1][0][1][1]=K[1][1][1][0]=K[0][1][1][1]=K[1][1][0][1]=tK[1]
	K[1][0][2][2]=K[2][2][1][0]=K[0][1][2][2]=K[2][2][0][1]=tK[2]
	K[1][0][1][0]=K[0][1][1][0]=K[1][0][0][1]=K[0][1][0][1]=tK[3]
	tK=udf.get("particle_type[%d].resistance_tensor.K.zx"%(index))
	K[2][0][0][0]=K[0][0][2][0]=K[0][2][0][0]=K[0][0][0][2]=tK[0]
	K[2][0][1][1]=K[1][1][2][0]=K[0][2][1][1]=K[1][1][0][2]=tK[1]
	K[2][0][2][2]=K[2][2][2][0]=K[0][2][2][2]=K[2][2][0][2]=tK[2]
	K[2][0][1][0]=K[0][2][1][0]=K[2][0][0][1]=K[0][2][0][1]=tK[3]
	K[1][0][2][0]=K[0][1][2][0]=K[1][0][0][2]=K[0][1][0][2]=tK[3]
	K[2][0][2][0]=K[0][2][2][0]=K[2][0][0][2]=K[0][2][0][2]=tK[4]
	tK=udf.get("particle_type[%d].resistance_tensor.K.zy"%(index))
	K[2][1][0][0]=K[0][0][2][1]=K[1][2][0][0]=K[0][0][1][2]=tK[0]
	K[2][1][1][1]=K[1][1][2][1]=K[1][2][1][1]=K[1][1][1][2]=tK[1]
	K[2][1][2][2]=K[2][2][2][1]=K[1][2][2][2]=K[2][2][1][2]=tK[2]
	K[2][1][1][0]=K[1][2][1][0]=K[2][1][0][1]=K[1][2][0][1]=tK[3]
	K[1][0][2][1]=K[0][1][2][1]=K[1][0][1][2]=K[0][1][1][2]=tK[3]
	K[2][1][2][0]=K[1][2][2][0]=K[2][1][0][2]=K[1][2][0][2]=tK[4]
	K[2][0][2][1]=K[0][2][2][1]=K[2][0][1][2]=K[0][2][1][2]=tK[4]
	K[2][1][2][1]=K[1][2][2][1]=K[2][1][1][2]=K[1][2][1][2]=tK[5]
	return K

#type is the string of ParticleType[i]
def getSubTensorTracelessK(udf,type):
	K=getSubTensorK(udf,type)
	deltaK=zeros((3,3),float)
	Kdelta=zeros((3,3),float)
	deltaKdelta=0
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					deltaK[i][j]=deltaK[i][j]+delta(k,l)*K[l][k][i][j]
					Kdelta[i][j]=Kdelta[i][j]+K[i][j][k][l]*delta(l,k)
					deltaKdelta=deltaKdelta+delta(i,j)*K[j][i][k][l]*delta(l,k)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					K[i][j][k][l]=K[i][j][k][l]-delta(i,j)*deltaK[k][l]/3.-Kdelta[i][j]*delta(k,l)/3.+delta(i,j)*deltaKdelta*delta(k,l)/9.
	return K

#get resistance tensor of particle_type[index] in udf
def getResistanceTensor(udf,index):
	return [getSubTensorA(udf,index),getSubTensorB(udf,index),getSubTensorC(udf,index),
	        getSubTensorG(udf,index),getSubTensorH(udf,index),getSubTensorK(udf,index)]

#cronecker's delta
def delta(i,j):
	if i==j:
		return 1.0
	else:
		return 0.0

def deltaTensor():
	return array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])

#Levi-chivita's symbol
def epsilon(i,j,k):
	if i==j or j==k or k==i:
		return 0
	if i==0:
		if j==1:
			return 1.
		else:
			return -1.
	if i==1:
		if j==0:
			return -1.
		else:
			return 1.
	if j==1:
		return -1.
	else:
		return 1.

def epsilonTensor():
	return array([[[0.,0.,0.],[0.,0.,1.],[0.,-1.,0.]],
	              [[0.,0.,-1.],[0.,0.,0.],[1.,0.,0.]],
	              [[0.,1.,0.],[-1.,0.,0.],[0.,0.,0.]]])

#inverse of 3x3 array
#For avoiding the error using inverse() in Numerical python
def inverse3x3(a):
	det=a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[1][2]*a[2][0]+a[1][0]*a[2][1]*a[0][2] \
	  -(a[0][2]*a[1][1]*a[2][0]+a[0][1]*a[1][0]*a[2][2]+a[1][2]*a[2][1]*a[0][0])
	return array([[(a[1][1]*a[2][2]-a[1][2]*a[2][1])/det,
	               (a[0][2]*a[2][1]-a[0][1]*a[2][2])/det,
	               (a[0][1]*a[1][2]-a[1][1]*a[0][2])/det],
	              [(a[1][2]*a[2][0]-a[1][0]*a[2][2])/det,
	               (a[0][0]*a[2][2]-a[0][2]*a[2][0])/det,
	               (a[0][2]*a[1][0]-a[0][0]*a[1][2])/det],
	              [(a[1][0]*a[2][1]-a[1][1]*a[2][0])/det,
	               (a[0][1]*a[2][0]-a[0][0]*a[2][1])/det,
	               (a[0][0]*a[1][1]-a[0][1]*a[1][0])/det]])

#Q is the second rank tensor
#Shuould be |det(Q)|=1
#A,G and K are the second, third and forth tensor.
def convertTensor3x3(Q,A):
	Ap=zeros((3,3),float)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					Ap[i][j]=Ap[i][j]+Q[i][k]*Q[j][l]*A[k][l]
	return Ap

def convertTensor3x3x3(Q,G):
	Gp=zeros((3,3,3),float)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					for m in range(3):
						for n in range(3):
							Gp[i][j][k]=Gp[i][j][k]+Q[i][l]*Q[j][m]*Q[k][n]*G[l][m][n]
	return Gp

def convertTensor3x3x3x3(Q,K):
	Kp=zeros((3,3,3,3),float)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					for m in range(3):
						for n in range(3):
							for o in range(3):
								for p in range(3):
									Kp[i][j][k][l]=Kp[i][j][k][l]+Q[i][m]*Q[j][n]*Q[k][o]*Q[l][p]*K[m][n][o][p]
	return Kp

#average of tensor for all direction
def averageForDirection3x3(t):
	rr=zeros((3,3),float)
	for i in range(3):
		for j in range(3):
			for i2 in range(3):
				for j2 in range(3):
					rr[i][j]=rr[i][j]+t[i2][j2]*delta(i,j)*delta(i2,j2)/3.
	return rr

def averageForDirection3x3x3(t):
	rr=zeros((3,3,3),float)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for i2 in range(3):
					for j2 in range(3):
						for k2 in range(3):
							rr[i][j][k]=rr[i][j][k]+t[i2][j2][k2]*epsilon(i,j,k)*epsilon(i2,j2,k2)/6.
	return rr

def averageForDirection3x3x3x3(t):
	rr=zeros((3,3,3,3),float)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					for i2 in range(3):
						for j2 in range(3):
							for k2 in range(3):
								for l2 in range(3):
									s=delta(i2,j2)*delta(k2,l2)*(4.*delta(i,j)*delta(k,l)-delta(i,k)*delta(j,l)-delta(i,l)*delta(j,k))
									s=s+delta(i2,k2)*delta(j2,l2)*(-delta(i,j)*delta(k,l)+4*delta(i,k)*delta(j,l)-delta(i,l)*delta(j,k))
									s=s+delta(i2,l2)*delta(j2,k2)*(-delta(i,j)*delta(k,l)-delta(i,k)*delta(j,l)+4.*delta(i,l)*delta(j,k))
									rr[i][j][k][l]=rr[i][j][k][l]+t[i2][j2][k2][l2]*s/30.
	return rr

##############################
#r:translationVector(r=R1-R0)#
#R0:CoordinationOfOrigin     #
#R1:AfterCoordination        #
##############################
#Z=(A,B,C,G,H,K)
def translationResistanceTensor(Z,r):
	A,B,C,G,H,K=Z[0],Z[1],Z[2],Z[3],Z[4],Z[5]
	A1,B1,C1,G1,H1,K1=copy(A),copy(B),copy(C),copy(G),copy(H),copy(K)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				G1[i][j][k]=G1[i][j][k]-0.5*(r[i]*A[j][k]+A[i][k]*r[j])
				for l in range(3):
					B1[i][j]=B1[i][j]-epsilon(i,k,l)*r[k]*A[l][j]
					C1[i][j]=C1[i][j]+epsilon(j,k,l)*B[i][k]*r[l]-epsilon(i,k,l)*r[k]*B[j][l]
					K1[i][j][k][l]=K1[i][j][k][l]-0.5*(G[i][j][k]*r[l]+G[i][j][l]*r[k]+r[i]*G[l][k][j]+r[j]*G[i][k][l])
					K1[i][j][k][l]=K1[i][j][k][l]+0.25*(r[i]*A[j][k]*r[l]+r[i]*A[j][l]*r[k]+r[j]*A[i][k]*r[l]+r[j]*A[i][l]*r[k])
					for m in range(3):
						for n in range(3):
							C1[i][j]=C1[i][j]-epsilon(i,k,l)*epsilon(j,m,n)*r[k]*r[n]*A[l][m]
	#using a value of above B1
	for i in range(3):
		for j in range(3):
			for k in range(3):
				H1[i][j][k]=H1[i][j][k]-0.5*(r[j]*B1[k][i]+r[i]*B1[k][j])
				for l in range(3):
					for m in range(3):
						H1[i][j][k]=H1[i][j][k]+epsilon(k,l,m)*G[i][j][l]*r[m]
	return A1,B1,C1,G1,H1,K1

#z=(a,b,c,g,h,k,)
def translationMobilityTensor(z,r):
	a0,b0,c0,g0,h0,k0=copy(z[0]),copy(z[1]),copy(z[2]),copy(z[3]),copy(z[4]),copy(z[5])
	a1,b1,g1=copy(a0),copy(b0),copy(g0)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				g1[i][j][k]=g1[i][j][k]-0.5*(r[i]*delta(j,k)+delta(i,k)*r[j])
				for l in range(3):
					b1[i][j]=b1[i][j]+epsilon(j,k,l)*c0[i][k]*r[l]
					a1[i][j]=a1[i][j]+epsilon(j,k,l)*b0[k][i]*r[l]
					a1[i][j]=a1[i][j]-epsilon(i,k,l)*r[k]*b0[l][j]
					for m in range(3):
						g1[i][j][k]=g1[i][j][k]+epsilon(k,l,m)*h0[i][j][l]*r[m]
						for n in range(3):
							a1[i][j]=a1[i][j]-epsilon(i,k,l)*epsilon(j,m,n)*r[k]*c0[l][m]*r[n]
	return a1,b1,c0,g1,h0,k0

#convert from ABC to abc using Block matrix method
def convertABCToabc(A,B,C):
	eps=10.0**-15
	x=zeros((6,6),float)
	for i in range(3):
		for j in range(3):
			x[i  ][j  ]=A[i][j]
			x[i+3][j  ]=B[i][j]
			x[i  ][j+3]=B[j][i]
			x[i+3][j+3]=C[i][j]
	w,x=eig(x)
	if abs(w[0])<eps or abs(w[1])<eps or abs(w[2])<eps or abs(w[3])<eps or abs(w[4])<eps or abs(w[5])<eps:
		y=zeros((6,6),float)
		for aa in range(6):
			if w[aa]>eps:
				for i in range(6):
					for j in range(6):
						y[i][j]+=x[aa][i]*x[aa][j]/w[aa]
		a,b,c=zeros((3,3),float),zeros((3,3),float),zeros((3,3),float)
		for i in range(3):
			for j in range(3):
				a[i][j]=y[i][j]
				b[i][j]=y[i+3][j]
				c[i][j]=y[i+3][j+3]
	else:
		Bt=transpose(B)
		A_=inverse3x3(A)
		C_=inverse3x3(C)
		a=inverse3x3(A-dot(dot(Bt,C_),B))
		b=transpose(-dot(dot(C_,B),a))
		c=inverse3x3(C-dot(dot(B,A_),Bt))
	return a,b,c

#convert from resistance tensor to mobility tensor
def convertResistanceToMobilityTensor(Z):
	A,B,C,G,H,KK=Z[0],Z[1],Z[2],Z[3],Z[4],Z[5]
	a,b,c=convertABCToabc(A,B,C)
	g,h=zeros((3,3,3),float),zeros((3,3,3),float)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					g[i][j][k]=g[i][j][k]+G[i][j][l]*a[l][k]+H[i][j][l]*b[l][k]
					h[i][j][k]=h[i][j][k]+G[i][j][l]*b[k][l]+H[i][j][l]*c[l][k]
	kk=-copy(KK)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					for m in range(3):
						kk[i][j][k][l]=kk[i][j][k][l]+G[i][j][m]*g[l][k][m]+H[i][j][m]*h[l][k][m]
	return a,b,c,g,h,kk

def convertTensorMobilityToResistanceTensor(Z):
	return convertResistanceToMobilityTensor(Z)

#b and c are the sub mobility tensors
#return values is relative position from particle's center
def calcCenterOfMobility(b,c):
	eps=10.0**-15
	if abs(det(c))<eps:
		return zeros(3)
	x=zeros((3),float)
	t=b-transpose(b)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				x[i]=x[i]+epsilon(i,j,k)*t[k][j]
	r=0.5*dot(inverse3x3(trace(c)*deltaTensor()-c),x)
	return r

#A and B are the sub resistance tensors
#return values is relative position from particle's center
def calcCenterOfResistance(A,B):
	eps=10.0**-15
	if abs(det(A))<eps:
		return zeros(3)
	x=zeros((3),float)
	t=B-transpose(B)
	for i in range(3):
		for j in range(3):
			for k in range(3):
				x[i]=x[i]+epsilon(i,j,k)*t[k][j]
	r=0.5*dot(inverse3x3(trace(A)*deltaTensor()-A),x)
	return r


