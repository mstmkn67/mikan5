from MUtility import *
from LinearAlgebra import *
from math import *
#Z:resistance tensor
def linearViscosity(Z,kBT=1.0,omega1=10e-5,omega2=10e0,division=100):
	eps=10e-20
	z=convertResistanceToMobilityTensor(Z)#mobility tensors at center of resistance
	translate=calcCenterOfMobility(z[1],z[2])     #center of mobility
	a,b,c,g,h,k,p,q,w,x=translationMobilityTensor(z,translate)#mobility tensors at center of mobility
	evalues,evectors=eigenvectors(c)
	c,h,k=convertTensor3x3(evectors,c),convertTensor3x3x3(evectors,h),convertTensor3x3x3x3(evectors,k)
	c,trc=[c[0][0],c[1][1],c[2][2]],c[0][0]+c[1][1]+c[2][2]
	alpha,lam,D=zeros(5,Float),zeros(5,Float),0
	if abs(c[0]-c[1])>eps or abs(c[1]-c[2])>eps:
	#if abs((c[0]-c[1])/c[0])>eps or abs((c[1]-c[2])/c[1])>eps:
		if  trc*trc-3*(c[0]*c[1]+c[1]*c[2]+c[2]*c[1])>0.0:
			D=sqrt(trc*trc-3*(c[0]*c[1]+c[1]*c[2]+c[2]*c[1]))
		if abs((c[0]-c[1])/c[0])<eps:
			v1,v2=array([c[1]-c[2],c[0]-c[1]+D,-c[0]+c[2]-D]),array([c[1]-c[2],c[0]-c[1]-D,-c[0]+c[2]+D])
		else:
			v1,v2=array([-c[0]+c[2]+D,c[1]-c[2]-D,c[0]-c[1]]),array([-c[0]+c[2]-D,c[1]-c[2]+D,c[0]-c[1]])
		beta=array([2*(h[1][2][0]-h[2][1][0]),2*(h[2][0][1]-h[0][2][1]),2*(h[0][1][2]-h[1][0][2])])
		alpha[0],alpha[1]=dot(beta,v1)*dot(beta,v1)/dot(v1,v1),dot(beta,v2)*dot(beta,v2)/dot(v2,v2)
		alpha[2]=2*(h[1][2][1]-h[2][1][1]+h[2][0][0]-h[0][0][2])*(h[1][2][1]-h[2][1][1]+h[2][0][0]-h[0][0][2])
		alpha[3]=2*(h[2][0][2]-h[0][2][2]+h[0][1][1]-h[1][1][0])*(h[2][0][2]-h[0][2][2]+h[0][1][1]-h[1][1][0])
		alpha[4]=2*(h[0][1][0]-h[1][0][0]+h[1][2][2]-h[2][2][1])*(h[0][1][0]-h[1][0][0]+h[1][2][2]-h[2][2][1])
	lam[0],lam[1]=trc+D,trc-D
	lam[2],lam[3],lam[4]=0.5*(trc+3*c[2]),0.5*(trc+3*c[0]),0.5*(trc+3*c[1])
	#print "rotational characteristic times:\n",lam/kBT
	omega,etaReal,etaImag=[],[],[]
	lo1,lo2=log10(omega1),log10(omega2)
	for n in arange(lo1,lo2,(lo2-lo1)/division):
		o=10**n
		omega.append(o)
		eta=0.0
		for i in range(3):
			for j in range(3):
				eta=eta-0.1*(k[i][j][i][j]-k[i][i][j][j]/3.)
		for i in range(5):
			eta=eta+0.1*kBT*alpha[i]/(2*lam[i]*kBT+o*1J)
		etaReal.append(eta.real)
		etaImag.append(-eta.imag)
	return lam/kBT,omega,etaReal,etaImag
