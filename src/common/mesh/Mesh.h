#ifndef _MESH_H_
#define _MESH_H_

#include "../Vector3d.h"
#include <vector>
#include <cmath>
#include <string>
using namespace std;

class Vertex{
public:
	Vertex(){}
	Vertex(const Vector3d& r){position=r;}
	virtual ~Vertex(){}
	Vector3d position;
};

class Face{
public:
	Face() {}
	Face(Vertex* v0,Vertex* v1,Vertex* v2){vertex[0]=v0;vertex[1]=v1;vertex[2]=v2;}
	virtual ~Face(){}
	
	Vertex* vertex[3];
	Vector3d position;
	Vector3d normal;
	double area;
	string type;
	virtual void calc_position_area_normal(){
		const Vector3d& r0=vertex[0]->position;
		const Vector3d& r1=vertex[1]->position;
		const Vector3d& r2=vertex[2]->position;
		position=(r0+r1+r2)/3.0;
		normal=((r1-r0)^(r2-r0));
		area=0.5*normal.length();
		normal/=normal.length();
	}
};

class Mesh{
public:
	Mesh(){};
	virtual ~Mesh(){};
	virtual void calc_positions_areas_normals(){
		vector<Face*>::iterator i=face.begin();
		for(;i!=face.end();i++){
			(*i)->calc_position_area_normal();
		}
	}
	vector<Vertex*> vertex;
	vector<Face*> face;
};

#endif // _MESH_H_
