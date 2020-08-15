from Quaternion import rotate_,x_vector,y_vector,z_vector
from numpy import sqrt,array
#############################################################
#BEM               pos_vertex[id]=pos :dictionary 
def draw_BEM_model_concolor(udf,pos,angle,pos_vertex,index_face,attr):
	p={}
	for i in pos_vertex:
		p[i]=rotate_(angle,pos_vertex[i])
		p[i]=[p[i][0]+pos[0],p[i][1]+pos[1],p[i][2]+pos[2]]
	for i in index_face:
		udf.polygon([p[i[0]],p[i[1]],p[i[2]]],attr)
#  each face is drew by each color
def draw_BEM_model_color(udf,pos,angle,pos_vertex,index_face,attr):
	p={}
	for i in pos_vertex:
		p[i]=rotate_(angle,pos_vertex[i])
		p[i]=[p[i][0]+pos[0],p[i][1]+pos[1],p[i][2]+pos[2]]
	#draw face
	for i in range(len(index_face)):
		udf.polygon([p[index_face[i][0]],p[index_face[i][1]],p[index_face[i][2]]],attr[i])
##########################################################################
#beads
def draw_beads_model_concolor(udf,pos,angle,pos_bead,attr):
	for i in pos_bead:
		a=rotate_(angle,i)
		a=[a[0]+pos[0],a[1]+pos[1],a[2]+pos[2]]
		udf.sphere(a,attr+[1.0])
def draw_beads_model_color(udf,pos,angle,pos_bead,attr):
	for i in range(len(pos_bead)):
		a=rotate_(angle,pos_bead[i])
		a=[a[0]+pos[0],a[1]+pos[1],a[2]+pos[2]]
		udf.sphere(a,attr[i]+[1.0])
##########################################################################
#superposition
def draw_superposition_model(udf,pos,angle,element,color):
	n=len(element)
	#color
	co=[]
	if len(color)==0:
		for i in range(n):
			co.append([1.0,0.0,0.0,1.0])
	elif len(color)==1:
		for i in range(n):
			co.append(color[0])
	elif len(element)!=len(color):
		print("Numbers of element[] and color[] are different!")
		return
	else:
		co=color
	for i in range(n):
		ele_pos=rotate_(angle,element[i][0])
		ele_pos=[ele_pos[0]+pos[0],ele_pos[1]+pos[1],ele_pos[2]+pos[2]]
		if element[i][1]=="sphere":
			a=element[i][2][0]
			udf.sphere(ele_pos,co[i]+[a])
		elif element[i][1]=="prolate_spheroid":
			d,l,s=element[i][3][0],element[i][3][1],element[i][3][2]
			d=rotate_(angle,d)
			rr=sqrt(l*l-s*s)
			p1=[ele_pos[0]+rr*d[0],ele_pos[1]+rr*d[1],ele_pos[2]+rr*d[2]]
			p2=[ele_pos[0]-rr*d[0],ele_pos[1]-rr*d[1],ele_pos[2]-rr*d[2]]
			udf.ellipsoid2(p1,p2,co[i]+[l])
		elif element[i][1]=="oblate_spheroid":
			d,l,s=element[i][4][0],element[i][4][1],element[i][4][2]
			d=rotate_(angle,d)
			c1,c2=d[2]*d[2]+d[1]*d[1],d[2]*d[2]+d[0]*d[0]
			if c1>c2:
				d2=[0,d[2]/sqrt(c1),-d[1]/sqrt(c1)]
			else:
				d2=[-d[2]/sqrt(c2),0,d[0]/sqrt(c2)]
			udf.ellipsoid1(ele_pos,co[i]+[s,l,l,d[0],d[1],d[2],d2[0],d2[1],d2[2]])
		elif element[i][1]=="disk":
			d=element[i][5][0]
			d=rotate_(angle,d)
			a=element[i][5][1]
			udf.disk(ele_pos,co[i]+[a,d[0],d[1],d[2]])
		elif element[i][1]=="needle":
			d,l,s=element[i][6][0],element[i][6][1],element[i][6][2]
			d=rotate_(angle,d)
			rr=sqrt(l*l-s*s)
			p1=[ele_pos[0]+rr*d[0],ele_pos[1]+rr*d[1],ele_pos[2]+rr*d[2]]
			p2=[ele_pos[0]-rr*d[0],ele_pos[1]-rr*d[1],ele_pos[2]-rr*d[2]]
			udf.ellipsoid2(p1,p2,co[i]+[l])
		elif element[i][1]=="point":
			udf.point(ele_pos,co[i])
		elif element[i][1]=="line":
			d,h=element[i][7][0],element[i][7][1]
			d=rotate_(angle,d)
			r1=[ h*d[0]+ele_pos[0], h*d[1]+ele_pos[1], h*d[2]+ele_pos[2]]
			r2=[-h*d[0]+ele_pos[0],-h*d[1]+ele_pos[1],-h*d[2]+ele_pos[2]]
			udf.line(r1,r2,co[i])
#dice
def drawDice(udf,q,rc,l):
	a=0.5*l
	ex,ey,ez=array(x_vector(q)),array(y_vector(q)),array(z_vector(q))
	r=[a*(ex-ey-ez),a*(ex+ey-ez),a*(ex+ey+ez),a*(ex-ey+ez),
	  a*(-ex-ey-ez),a*(-ex+ey-ez),a*(-ex+ey+ez),a*(-ex-ey+ez)]
	r2=[[]]*8
	for i in range(8):
		r2[i]=(r[i]+rc).tolist()
	udf.polygon([r2[0],r2[1],r2[2],r2[3]],[1,1,1,1])
	udf.polygon([r2[1],r2[5],r2[6],r2[2]],[1,1,1,1])	
	udf.polygon([r2[3],r2[2],r2[6],r2[7]],[1,1,1,1])
	udf.polygon([r2[0],r2[1],r2[5],r2[4]],[1,1,1,1])
	udf.polygon([r2[0],r2[4],r2[7],r2[3]],[1,1,1,1])
	udf.polygon([r2[4],r2[5],r2[6],r2[7]],[1,1,1,1])
	e=1.01
	for i in range(8):
		r2[i]=(e*r[i]+rc).tolist()
	udf.polyline([r2[0],r2[1],r2[2],r2[3],r2[0]],[0,0,0,1])
	udf.polyline([r2[4],r2[5],r2[6],r2[7],r2[4]],[0,0,0,1])
	udf.line(r2[1],r2[5],[0,0,0,1]);udf.line(r2[2],r2[6],[0,0,0,1])
	udf.line(r2[0],r2[4],[0,0,0,1]);udf.line(r2[3],r2[7],[0,0,0,1])
	a,b,c,d=a*1.01,0.15*l,l/4.0,l/3.0
	#1
	udf.disk((a*ex+rc).tolist(),[1,0,0,1,b,ex[0],ex[1],ex[2]])
	#2
	udf.disk((a*ey+c*ex+c*ez+rc).tolist(),[0,0,0,1,b,ey[0],ey[1],ey[2]])
	udf.disk((a*ey-c*ex-c*ez+rc).tolist(),[0,0,0,1,b,ey[0],ey[1],ey[2]])
	#3
	udf.disk((a*ez+rc).tolist(),[0,0,0,1,b,ez[0],ez[1],ez[2]])
	udf.disk((a*ez+c*ex+c*ey+rc).tolist(),[0,0,0,1,b,ez[0],ez[1],ez[2]])
	udf.disk((a*ez-c*ex-c*ey+rc).tolist(),[0,0,0,1,b,ez[0],ez[1],ez[2]])
	#4
	udf.disk((-a*ez+c*ex+c*ey+rc).tolist(),[0,0,0,1,b,ez[0],ez[1],ez[2]])
	udf.disk((-a*ez+c*ex-c*ey+rc).tolist(),[0,0,0,1,b,ez[0],ez[1],ez[2]])
	udf.disk((-a*ez-c*ex+c*ey+rc).tolist(),[0,0,0,1,b,ez[0],ez[1],ez[2]])
	udf.disk((-a*ez-c*ex-c*ey+rc).tolist(),[0,0,0,1,b,ez[0],ez[1],ez[2]])
	#5
	udf.disk((-a*ey+rc).tolist(),[0,0,0,1,b,ey[0],ey[1],ey[2]])
	udf.disk((-a*ey+c*ex+c*ez+rc).tolist(),[0,0,0,1,b,ey[0],ey[1],ey[2]])
	udf.disk((-a*ey-c*ex+c*ez+rc).tolist(),[0,0,0,1,b,ey[0],ey[1],ey[2]])
	udf.disk((-a*ey+c*ex-c*ez+rc).tolist(),[0,0,0,1,b,ey[0],ey[1],ey[2]])
	udf.disk((-a*ey-c*ex-c*ez+rc).tolist(),[0,0,0,1,b,ey[0],ey[1],ey[2]])
	#6
	udf.disk((-a*ex-c*ey+d*ez+rc).tolist(),[0,0,0,1,b,ex[0],ex[1],ex[2]])
	udf.disk((-a*ex-c*ey+rc).tolist(),[0,0,0,1,b,ex[0],ex[1],ex[2]])
	udf.disk((-a*ex-c*ey-d*ez+rc).tolist(),[0,0,0,1,b,ex[0],ex[1],ex[2]])
	udf.disk((-a*ex+c*ey+d*ez+rc).tolist(),[0,0,0,1,b,ex[0],ex[1],ex[2]])
	udf.disk((-a*ex+c*ey+rc).tolist(),[0,0,0,1,b,ex[0],ex[1],ex[2]])
	udf.disk((-a*ex+c*ey-d*ez+rc).tolist(),[0,0,0,1,b,ex[0],ex[1],ex[2]])

def draw_frame(udf,smin,max,color=[1,1,1,1]):
	size=[max[0]-smin[0],max[1]-smin[1],max[2]-smin[2]]
	r0=[smin[0]        ,smin[1]        ,smin[2]]
	r1=[smin[0]+size[0],smin[1]        ,smin[2]]
	r2=[smin[0]+size[0],smin[1]+size[1],smin[2]]
	r3=[smin[0]        ,smin[1]+size[1],smin[2]]
	r4=[smin[0]        ,smin[1]        ,smin[2]+size[2]]
	r5=[smin[0]+size[0],smin[1]        ,smin[2]+size[2]]
	r6=[smin[0]+size[0],smin[1]+size[1],smin[2]+size[2]]
	r7=[smin[0]        ,smin[1]+size[1],smin[2]+size[2]]
	udf.polyline([r0,r1,r2,r3,r0],color)
	udf.polyline([r4,r5,r6,r7,r4],color)
	udf.line(r0,r4,color)
	udf.line(r1,r5,color)
	udf.line(r2,r6,color)
	udf.line(r3,r7,color)

def draw_xyz(udf,length=1.0):
	udf.line([0,0,0],[length,0,0],3)
	udf.line([0,0,0],[0,length,0],2)
	udf.line([0,0,0],[0,0,length],1)
