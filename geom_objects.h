#ifndef GEOM_OBJECTS_H
#define GEOM_OBJECTS_H 1

#include <iostream>
#include <math.h>
#include "const.h"

using namespace std;

enum Boolean {False,True};


class Vector{
           protected:
             double x,y,z;
             double norm;
	   public:
	     Vector(double x0=0., double y0=0., double z0=0.);
	     Vector(const Vector &v);
	     double operator *(Vector v);
	     Vector operator ^(Vector v);
	     Vector operator *(double t);
	     Vector operator +(Vector v);
	     Vector operator -(Vector v);
	     Vector operator -();
	     friend Vector operator *(double t,Vector v);
	     void print(){cout<<x<<"*"<<y<<"*"<<z<<endl;};
	     
	     void rotation(Vector axis,double cos_angle,double sin_angle);
	     void rotation(Vector axis,double angle)
	       {rotation(axis,cos(angle),sin(angle));};
	     void normalize(){
		norm=sqrt(x*x+y*y+z*z);
		if(norm==0.) return;
		*this=(*this)*(1./norm);
		norm=1.;};
	     double get_norm(){return(norm=sqrt(x*x+y*y+z*z));};
	     double get_x(){return(x);};
	     double get_y(){return(y);};
	     double get_z(){return(z);};	
	     void set_x(double xx){x=xx;};
	     void set_y(double yy){y=yy;};
	     void set_z(double zz){z=zz;};
	     void set(double xx,double yy,double zz){
		x=xx,y=yy,z=zz;
		norm=sqrt(x*x+y*y+z*z);};
	   };

 

class Line{
        protected:
         Vector passing_point;
	 Vector direction;
	 Vector point;
	 double parameter;

       public:
	 Vector get_passing_point(){return(passing_point);};
	 Vector get_direction(){
		direction.normalize();
		return direction;};
	 Vector get_point(){return(point);}
	 double get_parameter(){return(parameter);};
	 void set_passing_point(Vector pass){passing_point=pass;};
	 void set_direction(Vector dir){
	   direction=dir;
	   direction.normalize();
	 };
	 void set_point(Vector pt){point=pt;};
	 void set_parameter(double par){parameter=par;};

	 Line(){};
	 Line(Line & l)
	   {passing_point=l.get_passing_point();
	    direction=l.get_direction();
	    point=l.get_point();
	    parameter=l.get_parameter();};
	 Line(Vector pt,Vector direct)
	   {passing_point=pt;
	    direct.normalize();
	    direction=direct;
            point=passing_point;
            parameter=0.;
	  };
	 void translate(Vector new_pt){passing_point=new_pt;};
	 void rotation(Vector axis,double cos_angle,double sin_angle)
	   {direction.rotation(axis,cos_angle,sin_angle);};
	 void rotation(Vector axis,double angle)
	   {direction.rotation(axis,angle);};
	 void rotation(Vector new_direct)
	   { new_direct.normalize();
	     direction=new_direct;};
	 double intersect(Line l)
	   { return((passing_point-l.get_passing_point())*
		    (l.get_direction()-direction)/
		    ((l.get_direction()-direction)*
		     (l.get_direction()-direction)));
		  };

	 double map(Vector pt);
	 Vector map(double par);
	 void translate(double par){
	   Vector t;
	   t = map(par);
	   translate(t);};
       };


class Surface{
        protected: 
	 Vector local_xaxis,local_yaxis;
         Vector location;
	 Vector orientation;
	 Vector point;
	 double bound1,bound2;
	 double a[3][3];
	public:
	 Vector get_local_xaxis(void){return(local_xaxis);};
	 Vector get_local_yaxis(void){return(local_yaxis);};
	 void set_transform_matrix(void);
	 void set_local_xaxis(Vector axis){
	   local_xaxis = axis;
	   local_xaxis.normalize();
	   local_yaxis = orientation^local_xaxis;
	   set_transform_matrix();
	 };
	 void set_local_yaxis(Vector axis){
	   local_yaxis=axis;
	   local_yaxis.normalize();
	   local_xaxis = local_yaxis^orientation;
	   set_transform_matrix();
	 };
	 Vector get_location(void){return(location);};
	 Vector get_orientation(void){return(orientation);};
	 Vector get_point(void){return(point);};
	 void get_bound(double &b1, double &b2){b1=bound1;b2=bound2;};
	 void set_location(Vector loc){location=loc;};
	 void set_orientation(Vector orient){
	   orientation=orient;
	   orientation.normalize();
	 };
	 void set_point(Vector pt){point=pt;};
	 void set_bound(double b1, double b2){bound1=b1;bound2=b2;};

	 Surface():local_xaxis(0,0,1.){};
	 Surface(Surface& s)
	   {location=s.get_location();
	    orientation=s.get_orientation();
	    point=s.get_point();
	    local_xaxis=s.get_local_xaxis();
	    local_yaxis=s.get_local_yaxis();
	    s.get_bound(bound1,bound2);};
	 Surface(Vector pt,Vector orient):local_xaxis(0.,0.,1.)
	    {location=pt;
	    orient.normalize();
	    orientation=orient;
	    local_yaxis=orientation^local_xaxis;
            point=location;
	  };
	 Surface(Vector pt, Vector orient, Vector lx){
	    location=pt;
	    orient.normalize();
	    orientation=orient;
	    local_xaxis=lx;
	    local_yaxis=orientation^local_xaxis;};

	 void translate(Vector new_pt){location=new_pt;};
	 void rotation(Vector axis,double angle)
	   {orientation.rotation(axis,angle);};
	 void rotation(Vector axis,double cos_angle,double sin_angle)
	   {orientation.rotation(axis,cos_angle,sin_angle);};
      	 void rotation(Vector new_orient)
	   { new_orient.normalize();
	     orientation=new_orient;};
	 Vector get_local(Vector pt);
	 Vector get_global(Vector pt); 

	 virtual void deriv(double &d1, double &d2, Vector pt)=0;
			// pt is in local system.
	 virtual Boolean intersect(Line l,Vector &intersection, Vector &normal)
			=0;
		//intersection and normal returned in local system.
	 Boolean is_inside(Vector local_pt);//the parameter is in local system.
	 Boolean is_inside();//determine if point is inside.
       };

class Plane:public Surface{
       public:
         Plane():Surface(){};
	 Plane(Plane& p):Surface(p){};
         Plane(Vector pt,Vector orient):Surface(pt,orient){};
         Plane(Vector pt,Vector orient,Vector lx):Surface(pt,orient,lx){};
	 virtual void deriv(double &d1, double &d2, Vector pt){
		d1=0; d2=0;};
	 virtual Boolean intersect(Line l,Vector &intersection, 
				Vector &normal); 
       };


class Sphere:public Surface{
       protected:
         Vector center;
	 double radius;
       public:
	 Vector get_center(){return(center);};
	 double get_radius(){return(radius);};
	 void set_center(Vector c){center=c;};
	 void set_radius(double r){radius=r;};
	 void set_orientation(){
		orientation=center-location;
		orientation.normalize();
		};
         Sphere(Sphere& sph):Surface(sph)
	   {center=sph.get_center();
	    radius=sph.get_radius();};
	 Sphere(Vector loc,Vector ctr,double rad):Surface(loc,ctr-loc)
	   {center=ctr;
	    radius=rad;};
	 Sphere():Surface(){};
	 Sphere(Vector ctr,double rad):Surface(ctr-(rad/ctr.get_norm())
					       *ctr,ctr)
	  {center=ctr;
	   radius=rad;};
	 Sphere(Vector loc,Vector ctr):Surface(loc,ctr-loc),center(ctr){
	   radius=sqrt((ctr-loc)*(ctr-loc));};
	 virtual void deriv(double &d1, double &d2, Vector pt);
         virtual Boolean intersect(Line l,Vector &intersection, 
			Vector &normal);
       };



#endif

extern Vector Global_Xaxis;
extern Vector Global_Yaxis;
extern Vector Global_Zaxis;
extern Vector Origin; 




