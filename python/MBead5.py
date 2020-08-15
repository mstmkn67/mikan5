from math import *

def draw_beads(p,udf):
	for i in p:
		udf.sphere(i,[1,0,0,1,1])

def rect(x,y,z):
	p=[]
	for j in range(y):
		for i in range(x):
			p.append([2.0*i,2.0*j,0.0])
			if (z-1)!=0:
				p.append([2.0*i,2.0*j,2.0*(z-1)])
	for k in range(1,z-1):
		for i in range(x):
			p.append([2.0*i,0.0,2.0*k])
			if (y-1)!=0:
				p.append([2.0*i,2.0*(y-1),2.0*k])
	for k in range(1,z-1):
		for j in range(1,y-1):
			p.append([0.0,2.0*j,2.0*k])
			if (x-1)!=0:
				p.append([2.0*(x-1),2.0*j,2.0*k])
	return p

#helix
#a:radius
#b:pitch
#curvature:a/(a^2+b^2)
#torsion:b/(a^2+b^2)
#n:number of beads
#ch:chirality(1 or -1)
def helix(a,b,n,ch=1):
	p=[]
	c=sqrt(a*a+b*b)
	for i in range(n):
		s=2.0*i
		x=a*cos(s/c)
		y=ch*a*sin(s/c)
		z=b*s/c
		p.append([x,y,z])
	return p

#spherical shell
#l:radius
def sphere(l):
	N=int(pi*l*0.5)+1#number of rings
	p=[]
	for i in range(0,N+1):
		thetan=i*pi/(N)
		zn=l*cos(thetan)
		rn=l*sin(thetan)
		M=int(pi*rn)+1
		for j in range(0,M):
			thetam=j*2*pi/(M)
			xm=rn*cos(thetam)
			ym=rn*sin(thetam)
			p.append([xm,ym,zn])
	return p

#cylinderical shell
#l:length,r:radius
def cylinder(l,r):
	N=int(l*0.5)+1#number of rings
	p=[]
	for i in range(0,N):
		zn=float(i)*l/(N-1)
		M=int(pi*r)+1
		for j in range(0,M):
			thetam=j*2*pi/(M)
			xm=r*cos(thetam)
			ym=r*sin(thetam)
			p.append([xm,ym,zn])
	return p

#disk l:number of rings
def disk(l):
	aida=sqrt(3.0)
	N=l#number of rings
	p=[]
	for i in range(0,N+1):
		#rn=2*i
		rn=aida*i
		M=int(pi*rn)+1
		for j in range(0,M):
			thetam=j*2*pi/(M)
			xm=rn*cos(thetam)
			ym=rn*sin(thetam)
			p.append([xm,ym,0])
	return p

#p=rect(10,10,10)
#p=helix(3,3,100,-1)
#p=sphere(15)
#p=cylinder(20,5)
#p=disk(5)
#draw_beads(p,_udf_)
