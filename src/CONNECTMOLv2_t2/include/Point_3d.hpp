#ifndef __POINT_3D_HPP
#define __POINT_3D_HPP 
#include <iostream>
#include <cstdio>
#include <cmath>
#include "constants.hpp"

using namespace std;
/*
 * A 3d vector class
 */

class Point_3d{
private:
	void rotate(double *,double *,double *,double ); 
	float x,y,z;
	
public:

	Point_3d(float xin=0.0, float yin=0.0, float zin=0.0):x(xin),y(yin),z(zin)
	{
	}

	//returns dot product
	double dot_with(const Point_3d &p) const{
		return(x*p.x+y*p.y+z*p.z);
	}

	//returns cross product of 2 vectors
	Point_3d cross_with(const Point_3d &p) const{
		return(Point_3d(y*p.z-p.y*z,p.x*z-x*p.z, x*p.y-p.x*y));
	}
		
	//subtracts with vector p and returns the resultant vector
	Point_3d operator-(const Point_3d &p) const{
		return(Point_3d(x-p.x,y-p.y,z-p.z));
	}

	//adds vector p and returns the resultant
	Point_3d operator+(const Point_3d &p) const{
		return(Point_3d(x+p.x,y+p.y,z+p.z));
	}
	
	//returns sqrt(x*x+y*y+z*z) 
	double mod() const{
		return(sqrt(x*x+y*y+z*z));
	}

	//returns x*x+y*y+z*z 
	double mod_square() const{
		return(x*x+y*y+z*z);
	}

	//multiplies x,y,z, by s
	void scale_with(double s){
		x*=s;y*=s;z*=s;
	}
		
	//prints x,y,z values
	void showPosition() const{
		printf("x: %lf, y: %lf, z: %lf\n", x, y, z);
	}
	
	//returns angle between 2 vectors in radians
	double angle(const Point_3d &p) const{
		return acos( this->dot_with(p) / (sqrt (this->mod_square() * p.mod_square()) ) ); 
	}
	
	//returns the distance between 2 vector
	double distance(const Point_3d &p) const{
		return (((*this)-p).mod()); 
	}

	//rotates the point with an angle(in degree) around an arbitrary axis, axis is from pt1
	//to pt2.
	void rotate(const Point_3d &pt1, const Point_3d &pt2, double an);

	void setX(float xin) { x=xin; } 
	void setY(float yin) { y=yin; } 
	void setZ(float zin) { z=zin; }
        void setXYZ(float xin,float yin, float zin)  { x=xin;y=yin;z=zin; }
	float getX() const { return x;}
	float getY() const { return y;}
	float getZ() const { return z;}
	
};

inline void Point_3d::rotate(const Point_3d &pt1, const Point_3d &pt2, double an){
	
	double b[3],c[3],d[3];

	b[0] = pt1.x;	
	b[1] = pt1.y;	
	b[2] = pt1.z;

	c[0] = pt2.x;	
	c[1] = pt2.y;	
	c[2] = pt2.z;

	d[0] = x;
	d[1] = y;
	d[2] = z;	
	
	rotate(b,c,d,an);
	
	x = d[0];
	y = d[1];
	z = d[2];	
}

//rotates 3d coordinate d with angle x(in degrees) about the axis formed by points b and c
//Surojit's module
inline void Point_3d::rotate(double *b,double *c,double *d,double an){
	int i;
	double c1[3],d1[3];
	//double pt_b[3],pt_c[3],pt_d[3];
	double rc,rc1;
	double cos_c1_x,cos_c1_xy,sin_c1_x,sin_c1_xy;
	double d_x1,d_x,d_y,d_z,d_y1,d_z1;
	double d_z3,d_x4,d_y4,d_x3;//d_y3, d_z4;
	double phi;

	phi = (PI * an)/180;

	//STEP ONE
	for(i=0;i<3;i++){
		c1[i] = c[i] - b[i];
		d1[i] = d[i] - b[i];
	}


	//STEP TWO
	rc = sqrt((c1[0] * c1[0]) + (c1[1] * c1[1]) + (c1[2] * c1[2]));
	rc1 = sqrt((c1[0] * c1[0]) + (c1[1] * c1[1]));
	if(rc1==0){
//		cos_c1_x= 1/sqrt(2);
//		sin_c1_x= 1/sqrt(2);
		cos_c1_x= 1;
		sin_c1_x= 0;
		
	}
	else{
		cos_c1_x = c1[0]/rc1;
		sin_c1_x = c1[1]/rc1;
	}
	cos_c1_xy = rc1/rc;
	sin_c1_xy = c1[2]/rc;
	d_x1 = d1[1] * sin_c1_x + d1[0] * cos_c1_x;
	d_y = d1[1] * cos_c1_x - d1[0] * sin_c1_x;
	d_x = d1[2] * sin_c1_xy + d_x1 * cos_c1_xy;
	d_z = d1[2] * cos_c1_xy - d_x1 * sin_c1_xy;
	d_y1 = d_z * sin(phi) + d_y * cos(phi);
	d_z1 = d_z * cos(phi) - d_y * sin(phi);


	//STEP THREE
	d_x3 = d_z1 *(-sin_c1_xy) + d_x * cos_c1_xy;
	d_z3 = d_z1 * cos_c1_xy - d_x * (-sin_c1_xy);
	d_x4 = d_y1 *(- sin_c1_x) + d_x3 * cos_c1_x;
	d_y4 = d_y1 * cos_c1_x - d_x3 * (-sin_c1_x);

	d[0]=d_x4;
	d[1]=d_y4;
	d[2]=d_z3;
	for(i=0;i<3;i++){
		d[i]=d[i]+b[i];
	}
}

#endif
