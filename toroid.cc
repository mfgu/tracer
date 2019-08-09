#include <math.h>
#include <iostream>
#include "toroid.h"

const double TINY=1.e-8;

void Toroid::deriv(double &d1, double &d2, Vector pt)
{
  double x,y,z;
  x=pt.get_x();
  y=pt.get_y();
  z=pt.get_z();
  d2=y/(horizontal_radius-z);
  d1=x/(horizontal_radius-z);
  d1/=(1.-(horizontal_radius-vertical_radius)/sqrt(y*y
		+(horizontal_radius-z)*(horizontal_radius-z)));
}

Boolean Toroid::intersect(Line l,Vector &intersection, Vector &normal)
{
  double a,b,c,t,t0,r;
  double x,y,z,x0,y0,z0,dir_x,dir_y,dir_z;
  Vector temp;
  t0=0.;
  r=vertical_radius/horizontal_radius;
  temp=l.get_passing_point();
  temp=get_local(temp);
  x0=temp.get_x();
  y0=temp.get_y();
  z0=temp.get_z();
  temp=l.get_direction()+location;
  temp=get_local(temp);
  temp.normalize();
  dir_x=temp.get_x();
  dir_y=temp.get_y();
  dir_z=temp.get_z();
  a=dir_x*dir_x+dir_y*dir_y*r+dir_z*dir_z;
  b=-2*(dir_x*x0+dir_y*y0*r-dir_z*(vertical_radius-z0));
  c=(vertical_radius-z0)*(vertical_radius-z0);
  c+=x0*x0+r*y0*y0-vertical_radius*vertical_radius;
  c=b*b-4*a*c;
  if(c<0.){
    normal=Origin;
    return(False);
  }
  t=(b+sqrt(c))/2/a;
  
  while(my_abs(t-t0)>t*TINY){
    t0=t;
    y=dir_y*t0+y0;
    z=dir_z*t0+z0;
    b=-2.*(dir_x*x0+dir_y*y0-(horizontal_radius-z0)*dir_z);
    c=(horizontal_radius-z0)*(horizontal_radius-z0);
    c+=x0*x0+y0*y0-horizontal_radius*horizontal_radius;
    a=(horizontal_radius-z)*(horizontal_radius-z)+y*y;
    a=sqrt(a);
    a=horizontal_radius-a;
    c+=2*(horizontal_radius-vertical_radius)*a;
    c=b*b-4*c;
    if(c<0.){
      normal=Origin;
      return(False);
    }
    t=(b+sqrt(c))/2.;
  }
  x=dir_x*t+x0;
  y=dir_y*t+y0;
  z=dir_z*t+z0;
  r=sqrt(vertical_radius*vertical_radius-x*x);
  r=horizontal_radius-vertical_radius+r;
  r=y/r;
  x0=0.;
  y0=(horizontal_radius-vertical_radius)*r;
  z0=(horizontal_radius-vertical_radius)*sqrt(1.-r*r);
  z0=horizontal_radius-z0;	
  normal.set(x0-x,y0-y,z0-z);
  normal.normalize();
  intersection.set(x,y,z);
  return(is_inside(intersection));
}

 

