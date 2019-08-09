
#include "opt_devices.h"
#include <iostream>

extern void prn(Vector);

double Mirror::efficiency()
{
  double in_angle,cos_ang,f1,f2,y;
  double wavelength,eff;
  double critical_angle;
  cos_ang=(incoming.get_direction())*(outgoing.get_direction());
  in_angle=acos(cos_ang)/2.; 
  wavelength=incoming.get_wavelength();
  f1f2(wavelength,f1,f2);
  critical_angle=sqrt(4.069e-5*2.*f1)*wavelength/12.4/A2cm;
  y=f2/f1;
  eff=reflectr(in_angle/critical_angle,y);
  return(eff);
} 


void Mirror::make_outgoing()
{
  Vector d_in,d_out;
  double inx,iny,inz; 
  double outx,outy,outz;
  double deriv1,deriv2;
  double a,b,c,delta,t1,t2;
  Vector intersection,normal;
  d_in=incoming.get_direction();
  inside=surf->intersect(incoming,intersection,normal);
  surf->deriv(deriv1,deriv2,intersection);
  d_in=d_in+surf->get_location();
  d_in=surf->get_local(d_in);
  inx=d_in.get_x();
  iny=d_in.get_y();
  inz=d_in.get_z();
  a=1.+deriv1*deriv1+deriv2*deriv2;
  t1=inx+inz*deriv1;
  t2=iny+inz*deriv2;
  b=2*(t1*deriv1+t2*deriv2);
  c=t1*t1+t2*t2-1;
  delta=b*b-4*a*c;
  outz=(-b-sqrt(delta))/(2*a);
  outx=-(t1+outz*deriv1);
  outy=-(t2+outz*deriv2);
  d_out.set(-outx,-outy,-outz);
  if((d_out*normal)*(d_in*normal)>0.){
  	outz=(-b+sqrt(delta))/(2*a);
	outx=-(t1+outz*deriv1);
	outy=-(t2+outz*deriv2);
	d_out.set(-outx,-outy,-outz);
  }
  outgoing.set_passing_point(surf->get_global(intersection));
  outgoing.set_direction(surf->get_global(d_out)-surf->get_location());
  outgoing.set_wavelength(incoming.get_wavelength());
}	
 
double Grating::efficiency()
{ 
  double in_angle,out_angle,cos_ang,alpha,beta,f1,f2,y;
  double g,p,Q,wavelength,eff;
  double critical_angle;
  Vector dir_in,dir_out,intersection,normal;
  dir_in=incoming.get_direction();
  dir_out=outgoing.get_direction();
  inside=surf->intersect(incoming,intersection,normal);
  normal=surf->get_global(normal);
  normal=normal-surf->get_location();
  normal.normalize();
  alpha=my_abs(asin(normal*dir_in));
  beta=my_abs(asin(normal*dir_out));
  in_angle=alpha+param.blazing_angle;
  out_angle=beta-param.blazing_angle;
  wavelength=incoming.get_wavelength();
  f1f2(wavelength,f1,f2);
  critical_angle=sqrt(4.069e-5*2.*f1)*wavelength/12.4/A2cm;
  y=f2/f1;
  p=sin(in_angle)+sin(out_angle);
  g=sin(alpha)/sin(in_angle);
  Q=Pi*g*param.m*(cos(in_angle)-cos(out_angle))/(cos(alpha)-cos(beta));
  Q=my_abs(Q);
  if(Q>1e-4)
	eff=(g*g*p*p*sin(Q)*sin(Q)/(sin(alpha)*sin(beta)*Q*Q*4.));
  else
	eff=sin(alpha)/sin(beta);
  eff*=sqrt(reflectr(in_angle/critical_angle,y) * 
	    reflectr(out_angle/critical_angle,y));

  return(eff);
} 

void Grating::make_outgoing()
{
  Vector d_in,d_out;
  double inx,iny,inz; 
  double outx,outy,outz;
  double deriv1,deriv2,deriv_n1,deriv_n2;
  double a,b,c,delta,t1,t2,sign;
  Vector intersection,normal;
  d_in=incoming.get_direction();
  inside=surf->intersect(incoming,intersection,normal);
  surf->deriv(deriv1,deriv2,intersection);
  deriv_n(deriv_n1,deriv_n2,param,intersection);
  d_in=d_in+surf->get_location();
  d_in=surf->get_local(d_in);
  inx=d_in.get_x();
  iny=d_in.get_y();
  inz=d_in.get_z();
  a=1.+deriv1*deriv1+deriv2*deriv2;
  t1=inx+inz*deriv1+deriv_n1*param.m*incoming.get_wavelength();
  t2=iny+inz*deriv2+deriv_n2*param.m*incoming.get_wavelength();
  b=2*(t1*deriv1+t2*deriv2);
  c=t1*t1+t2*t2-1;
  if(param.is_reflective) sign=1; 
  else sign=-1;
  delta=b*b-4*a*c;
  outz=(-b-sign*sqrt(delta))/(2*a);
  outx=-(t1+outz*deriv1);
  outy=-(t2+outz*deriv2);
  d_out.set(-outx,-outy,-outz);
  if((d_out*normal)*(d_in*normal)*sign>0.){
  	outz=(-b+sign*sqrt(delta))/(2*a);
	outx=-(t1+outz*deriv1);
	outy=-(t2+outz*deriv2);
	d_out.set(-outx,-outy,-outz);
  }
  outgoing.set_passing_point(surf->get_global(intersection));
  outgoing.set_direction(surf->get_global(d_out)-surf->get_location());
  outgoing.set_wavelength(incoming.get_wavelength());
}	
 

void Grating::deriv_n(double &d1, 
		      double &d2, RulingParams p, Vector pt)
{
 double x,y,z,yy,t1,t2,deriv1,deriv2;
 x=pt.get_x();
 y=pt.get_y();
 z=pt.get_z();
 yy=y-z*p.tan;
 surf->deriv(deriv1,deriv2,pt);
 t1=-deriv1*p.tan;
 t2=1-deriv2*p.tan;
 d2=(1.+2.*p.b2*yy+3*p.b3*yy*yy+4*p.b4*yy*yy*yy)*t2/p.spacing_at_center; 
 d1=(1.+2.*p.b2*yy+3*p.b3*yy*yy+4*p.b4*yy*yy*yy)*t1/p.spacing_at_center; 
}


void Detector::make_outgoing()
{
 Vector intersection,normal;
 inside=surf->intersect(incoming,intersection,normal);
 outgoing=incoming;
 outgoing.set_passing_point(surf->get_global(intersection));
}



