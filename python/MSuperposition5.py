def sphere(radius):
	return [[[0,0,0],"sphere",[radius]]]

def prolate_spheroid(director,long_axis,short_axis):
	return [[[0,0,0],"prolate_spheroid",[director,long_axis,short_axis]]]

def oblate_spheroid(director,long_axis,short_axis):
	return [[[0,0,0],"oblate_spheroid",[director,long_axis,short_axis]]]

def disk(director,radius):
	return [[[0,0,0],"disk",[director,radius]]]

def needle(director,long_axis,short_axis):
	return [[[0,0,0],"needle",[director,long_axis,short_axis]]]

def dumbbell(radius0,radius1,half_length):
	e0=[[0,0, half_length],"sphere",[radius0]]
	e1=[[0,0,-half_length],"sphere",[radius1]]
	e2=[[0.0,0.0,0.0],"line",[[0,0,1],half_length]]
	return [e0,e1,e2]

def propeller(radius0,radius1,twist_angle,half_length):
	import numpy as np
	t=0.5*twist_angle/180*np.pi
	u1,u2=np.array([1,0,0]),np.array([0,1,0])
	director0= np.sin(t)*u1+np.cos(t)*u2
	director1=-np.sin(t)*u1+np.cos(t)*u2
	e0=[[0,0,-half_length],"disk",[director0.tolist(),radius0]]
	e1=[[0,0, half_length],"disk",[director1.tolist(),radius1]]
	e2=[[0.0,0.0,0.0],"line",[[0,0,1],half_length]]
	return [e0,e1,e2]
