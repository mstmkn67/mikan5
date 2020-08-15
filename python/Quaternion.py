#Quaternion.py
#functions for quaternion
from math import *
#From Euler angle(phi,theta,psi) to quaternion
def eulerToFour(angle):
	phi=angle[0]
	theta=angle[1]
	psi=angle[2]
	xi=sin(theta/2)*sin((psi-phi)/2)
	eta=sin(theta/2)*cos((psi-phi)/2)
	zeta=cos(theta/2)*sin((psi+phi)/2)
	chi=cos(theta/2)*cos((psi+phi)/2)
	return [xi,eta,zeta,chi]

#rotation 3d vector r using quaternion q
def rotate(q,r):
	xi,eta,zeta,chi=q[0],q[1],q[2],q[3]
	x=r[0]*(-xi*xi+eta*eta-zeta*zeta+chi*chi)+r[1]*2*(zeta*chi-xi*eta)+r[2]*2*(eta*zeta+xi*chi)
	y=r[0]*(-2*(xi*eta+zeta*chi))+r[1]*(xi*xi-eta*eta-zeta*zeta+chi*chi)+r[2]*2*(eta*chi-xi*zeta)
	z=r[0]*2*(eta*zeta-xi*chi)+r[1]*(-2*(xi*zeta+eta*chi))+r[2]*(-xi*xi-eta*eta+zeta*zeta+chi*chi)
	return [x,y,z]

#rotation tensor using quaternion q
def rotateTensor(q):
	xi,eta,zeta,chi=q[0],q[1],q[2],q[3]
	return [[-xi*xi+eta*eta-zeta*zeta+chi*chi,2*(zeta*chi-xi*eta),2*(eta*zeta+xi*chi)],
	        [-2*(xi*eta+zeta*chi),xi*xi-eta*eta-zeta*zeta+chi*chi,2*(eta*chi-xi*zeta)],
	        [2*(eta*zeta-xi*chi),-2*(xi*zeta+eta*chi),-xi*xi-eta*eta+zeta*zeta+chi*chi]]


#inverse rotation 3d vector r using quaternion q
def rotate_(q,r):
	xi,eta,zeta,chi=q[0],q[1],q[2],q[3]
	x=r[0]*(-xi*xi+eta*eta-zeta*zeta+chi*chi)+r[1]*(-2*(xi*eta+zeta*chi))+r[2]*2*(eta*zeta-xi*chi)
	y=r[0]*2*(zeta*chi-xi*eta)+r[1]*(xi*xi-eta*eta-zeta*zeta+chi*chi)+r[2]*(-2*(xi*zeta+eta*chi))
	z=r[0]*2*(eta*zeta+xi*chi)+r[1]*2*(eta*chi-xi*zeta)+r[2]*(-xi*xi-eta*eta+zeta*zeta+chi*chi)
	return [x,y,z]

#inverse rotation tensor using quaternion q
def rotate_Tensor(q):
	xi,eta,zeta,chi=q[0],q[1],q[2],q[3]
	return [[-xi*xi+eta*eta-zeta*zeta+chi*chi,-2*(xi*eta+zeta*chi),2*(eta*zeta-xi*chi)],
	        [2*(zeta*chi-xi*eta),xi*xi-eta*eta-zeta*zeta+chi*chi,-2*(xi*zeta+eta*chi)],
	        [2*(eta*zeta+xi*chi),2*(eta*chi-xi*zeta),-xi*xi-eta*eta+zeta*zeta+chi*chi]]

#x vector in zero Euler angle from quaternion q
def x_vector(q):
	xi,eta,zeta,chi=q[0],q[1],q[2],q[3]
	xx=-xi*xi+eta*eta-zeta*zeta+chi*chi
	yy=2*(-xi*eta+zeta*chi)
	zz=2*(eta*zeta+xi*chi)
	return [xx,yy,zz]

#y vector in zero Euler angle from quaternion q
def y_vector(q):
	xi,eta,zeta,chi=q[0],q[1],q[2],q[3]
	xx=-2*(zeta*chi+xi*eta)
	yy=xi*xi-eta*eta-zeta*zeta+chi*chi
	zz=2*(-xi*zeta+eta*chi)
	return [xx,yy,zz]

#z vector in zero Euler angle from quaternion q
def z_vector(q):
	xi,eta,zeta,chi=q[0],q[1],q[2],q[3]
	xx=2*(eta*zeta-xi*chi)
	yy=-2*(eta*chi+xi*zeta)
	zz=-xi*xi-eta*eta+zeta*zeta+chi*chi
	return [xx,yy,zz]

