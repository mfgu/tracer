#include <math.h>
#include "geom_objects.h"
#include <iostream>

Vector Global_Xaxis(1.,0.,0.);
Vector Global_Yaxis(0.,1.,0.);
Vector Global_Zaxis(0.,0.,1.);
Vector Origin(0.,0.,0.);

void prn(Vector p)
{cout<<"("<<p.get_x()<<"*"<<p.get_y()<<"*"<<p.get_z()<<")"<<endl;
}


Vector::Vector(double x0, double y0, double z0)
{
  x=x0;
  y=y0;
  z=z0;
  norm=sqrt(x*x+y*y+z*z);
}


Vector::Vector(const Vector &v0)
{
  Vector v;
  v = v0;
  x=v.get_x();
  y=v.get_y();
  z=v.get_z();
  norm=sqrt(x*x+y*y+z*z);
}

double Vector::operator *(Vector v)
{
  return(x*v.get_x()+y*v.get_y()+z*v.get_z());
}



Vector Vector::operator ^(Vector v)
{
  Vector temp;
  double xx,yy,zz;
  xx=y*(v.get_z())-z*(v.get_y());
  yy=z*(v.get_x())-x*(v.get_z());
  zz=x*(v.get_y())-y*(v.get_x());
  temp.set(xx,yy,zz);
  temp.get_norm();
  return temp;
}



Vector Vector::operator *(double t)
{
  Vector temp;
  temp.set(x*t,y*t,z*t);
  temp.get_norm();
  return temp;
}



Vector Vector::operator +(Vector v)
{
  Vector temp;
  temp.set(x+v.get_x(),y+v.get_y(),z+v.get_z());
  temp.get_norm();
  return temp;
}


Vector Vector::operator -(Vector v)
{
  Vector temp;
  temp.set(x-v.get_x(),y-v.get_y(),z-v.get_z());
  temp.get_norm();
  return temp;
}


Vector Vector::operator -()
{
 Vector temp;
 temp.set(-x,-y,-z);
 temp.get_norm();
 return temp;
}




Vector operator *(double t,Vector v)
{
  Vector temp;

  temp = v*t;
  return temp;
}



void Vector::rotation(Vector axis,double cos_angle,double sin_angle)
{
  Vector temp1,tempx,tempy,projection;
  if(this->get_norm()==0) return;
  if(axis.get_norm()==0) return ;
  axis.normalize();
  projection=(axis*(*this))*axis;
  tempx=temp1=*this-projection;

  if(temp1.get_norm()==0) return;
  tempx.normalize();
  tempy=axis^tempx;

  temp1=temp1.get_norm()*cos_angle*tempx+temp1.get_norm()*sin_angle*tempy;
  *this=projection+temp1;
  };




double Line::map(Vector pt)
{
  return(parameter=(pt-passing_point)*direction);
}

Vector Line::map(double par)
{
  parameter=par;
  return(point=passing_point+direction*par);
}


Boolean Surface::is_inside(Vector local_pt)
{
  if(my_abs(local_pt.get_x())<=bound1/2. && my_abs(local_pt.get_y())<=bound2/2.) 
	return(True);
  else return(False);
}

void Surface::set_transform_matrix(void) {
  a[0][0] = local_xaxis * Global_Xaxis;
  a[0][1] = local_yaxis * Global_Xaxis;
  a[0][2] = orientation * Global_Xaxis;
  a[1][0] = local_xaxis * Global_Yaxis;
  a[1][1] = local_yaxis * Global_Yaxis;
  a[1][2] = orientation * Global_Yaxis;
  a[2][0] = local_xaxis * Global_Zaxis;
  a[2][1] = local_yaxis * Global_Zaxis;
  a[2][2] = orientation * Global_Zaxis;
}

Boolean Surface::is_inside()
{
  return(is_inside(get_local(point)));
}

Vector Surface::get_global(Vector pt) {
  Vector temp;
  double x, y, z;

  x = pt.get_x();
  y = pt.get_y();
  z = pt.get_z();

  temp.set_x(a[0][0]*x + a[0][1]*y + a[0][2]*z);
  temp.set_y(a[1][0]*x + a[1][1]*y + a[1][2]*z);
  temp.set_z(a[2][0]*x + a[2][1]*y + a[2][2]*z);
  
  return temp + location;
}

Vector Surface::get_local(Vector pt) {
  Vector temp;
  double x, y, z;
  
  temp = pt - location;
  x = temp.get_x();
  y = temp.get_y();
  z = temp.get_z();

  temp.set_x(a[0][0]*x + a[1][0]*y + a[2][0]*z);
  temp.set_y(a[0][1]*x + a[1][1]*y + a[2][1]*z);
  temp.set_z(a[0][2]*x + a[1][2]*y + a[2][2]*z);
  
  return temp;
} 

Boolean Plane::intersect(Line l,Vector &intersection,Vector &normal)
{
  l.set_parameter(orientation*(location-l.get_passing_point())
		  /(l.get_direction()*orientation));
  point=l.map(l.get_parameter());
  normal.set(0.,0.,1);
  intersection=get_local(point); 
  return(is_inside(intersection));
}

Boolean Sphere::intersect(Line l, Vector &intersection, Vector &normal)
{
  double b,c;
  Vector temp;
  temp=l.get_passing_point()-center;
  b=l.get_direction()*temp*2;
  c=temp*temp-radius*radius;
  if((c=b*b-4*c)<0) 
    {normal=Origin;
    return(False);}
  else              
    {b=-b/2+sqrt(c)/2;
    point=l.map(b);
    temp=point;
    normal=center-temp;
    normal=get_local(location+normal);
    normal.normalize();
    intersection=get_local(temp);
    return(is_inside(intersection));
    }
}


void Sphere::deriv(double &d1, double &d2, Vector pt)
{
  double x,y,z;
  x=pt.get_x();
  y=pt.get_y();
  z=pt.get_z();
  d1=x/(radius-z);
  d2=y/(radius-z);
}

