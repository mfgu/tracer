#include "spectrometer.h"
#include <fstream>
#include <stdio.h>

struct SmRGSParams{
  double mirror_grating_separation;
  double mirror_grazing_angle;
  double grating_grazing_angle;
  };

class SmRGS :public Spectrometer{
 public:
  Sphere mirror_surface;
  Plane grating_blank,det_surface;
  SmRGSParams more_param;
  Mirror sphere_mirror;
  Grating plane_grating;
  Detector plane_detector;
  SmRGS():Spectrometer(){};
  virtual void read_input(void);
  virtual void optimize();
  virtual void write_params(void);
  virtual void setup_bench();
  virtual double resolving_power(double w);
}; 



double SmRGS::resolving_power(double w)
{
  double alpha, beta, delta_alpha, delta_beta, delta_w;
  RulingParams p;
  p=plane_grating.get_param();
  alpha=more_param.grating_grazing_angle-p.blazing_angle;
  beta=more_param.grating_grazing_angle+p.blazing_angle;
  delta_alpha=param.source_size/param.source_distance;
  delta_beta=param.det_resolution/param.focal_length;
  delta_w=p.spacing_at_center/my_abs(p.m);
  delta_w*=sqrt(sin(alpha)*sin(alpha)*delta_alpha*delta_alpha
	+sin(beta)*sin(beta)*delta_beta*delta_beta);
  return(w/delta_w);
}

void SmRGS::read_input(void)
{
  ifstream infile;
  char dummy[80];
  double x,y;
  RulingParams p;
  infile.open("mirrg.in",ios::in);
  if(!infile) cout<<"can't open file "<<"torg.in"<<endl;
  infile.getline(dummy,sizeof(dummy));
  (infile>>param.principal_resolving_power
	>>param.total_length).getline(dummy,sizeof(dummy));
  param.total_length*=100;

  infile.getline(dummy,sizeof(dummy));
  (infile>>p.spacing_at_center
	>>param.principal_wavelength).getline(dummy,sizeof(dummy));
  p.spacing_at_center=0.1
	/p.spacing_at_center;
  param.principal_wavelength*=A2cm;

  infile.getline(dummy,sizeof(dummy));
  (infile>>param.source_size
	>>param.det_resolution).getline(dummy,sizeof(dummy));
  param.source_size*=MICRON;
  param.det_resolution*=MICRON;

  infile.getline(dummy,sizeof(dummy));
  (infile>>p.m>>param.wavelength_min).getline(dummy,sizeof(dummy));
  param.wavelength_min*=A2cm;

  infile.getline(dummy,sizeof(dummy));
  (infile>>x>>y).getline(dummy,sizeof(dummy));
  mirror_surface.set_bound(x,y);
  infile.getline(dummy,sizeof(dummy));
  (infile>>x>>y).getline(dummy,sizeof(dummy));
  grating_blank.set_bound(x,y);
  infile.getline(dummy,sizeof(dummy));
  (infile>>x>>y).getline(dummy,sizeof(dummy));
  det_surface.set_bound(x,y);

  infile.getline(dummy,sizeof(dummy));
  (infile>>more_param.mirror_grating_separation).getline(dummy,sizeof(dummy));

  infile.getline(dummy,sizeof(dummy));
  (infile>>param.grasp_file
	>>param.calib_file).getline(dummy,sizeof(dummy));
  infile.getline(dummy,sizeof(dummy));
  (infile>>param.respower_file
	>>param.image_file).getline(dummy,sizeof(dummy));
  infile.getline(dummy,sizeof(dummy));
  (infile>>param.param_file).getline(dummy,sizeof(dummy));

  p.is_reflective=True;
  p.tan=0.;
  plane_grating.set_param(p);
}


void SmRGS::optimize()
{
  double critical_graz,inc_ang;
  double grazm,graz,blaz,d,s;
  double graz1min,graz1max,grazmin,grazmax,blazmin,blazmax;
  double bestgraz1,bestgraz,bestblaz,bestd;
  double alpha,beta,A,B,C,D;
  double f1,f2;
  double y;
  double eff=0,besteff=0;
  double spacing,blazing_wavelength;
  double msize_x,msize_y,gsize_x,gsize_y,dsize_x,dsize_y;
  int n,m;
  RulingParams p;

  p=plane_grating.get_param();
  mirror_surface.get_bound(msize_x,msize_y);
  grating_blank.get_bound(gsize_x,gsize_y);
  det_surface.get_bound(dsize_x,dsize_y);
  spacing=p.spacing_at_center;
  blazing_wavelength=param.principal_wavelength;
  s=more_param.mirror_grating_separation;
  m=-p.m;

  f1f2(blazing_wavelength,f1,f2);
  critical_graz=sqrt(4.069e-5*2*f1)*blazing_wavelength/12.4/A2cm;
  y=f2/f1;

  cout<<"optimizing...wait..."<<endl;
  
  inc_ang= 10e-4;
  graz1min=5e-4;
  graz1max=2.5*critical_graz;
//  cout<<"critical angle: "<<critical_graz<<endl;
  for(grazm=graz1min;grazm<graz1max;grazm+=inc_ang){
    grazmin=my_max(5e-4,
		   sqrt(m*blazing_wavelength/spacing/2.));
    grazmax=2.5*critical_graz;
    if(grazmax<grazmin) grazmax=5.*grazmin;
    for(graz=grazmin;graz<grazmax;graz+=inc_ang){
      blaz=(m*blazing_wavelength/spacing)/sin(graz)/2.;
      if(blaz>=1.) continue;
      blaz=asin(blaz);
      if(blaz>graz) continue;
      alpha=graz-blaz; 
      beta=graz+blaz;
      C=(beta*param.principal_resolving_power)/(2*graz*blaz);
      A=C*param.source_size*alpha/beta;
      B=C*param.det_resolution;
      D=param.total_length-s;
      if(D<(A+B)) continue;
      double Q=(-A+B-D)*(-A+B-D)-4*A*D;
      if(Q<0.0) continue;
      d=(A-B+D-sqrt(Q))/2;
      if(d<50) continue;
      eff=reflectr(grazm/critical_graz,y)*msize_y*grazm;
      eff*=sin(alpha)/sin(beta)/d;
      if(gsize_y*alpha<msize_y*grazm)
	eff*=alpha*gsize_y/(msize_y*grazm);
      eff*=reflectr(graz/critical_graz,y);
      if(eff>besteff){
	besteff=eff;
	bestgraz1=grazm;
	bestgraz=graz;
	bestblaz=blaz;
	bestd=d;
      }
    }
  }
  if(eff==0){
    cout<<"can't find proper parameters for the spectrometer."<<endl;
    exit(1);
  }
  more_param.mirror_grazing_angle=grazm=bestgraz1;
  more_param.grating_grazing_angle=graz=bestgraz;
  blaz=p.blazing_angle=bestblaz;
  param.source_distance=d=bestd;
  eff=besteff*dsize_x/param.total_length;
  param.focal_length=param.total_length-d-s;
//cout<<"the optimum grasp is: "<<eff<<endl;

 double a,b,c,rad,t,q,u,v;
 double C2,C3,C4,H1,H2,H3,B01,B10,B20,B02,B11,B30,B03,B12,B21;
 double cos_grazm,sin_grazm,cos_alpha,sin_alpha,cos_beta,sin_beta;
 double T1,T2;
 double L;
 L=param.focal_length;
 alpha=graz-blaz;
 beta=graz+blaz;
 cos_grazm=cos(grazm);
 sin_grazm=sin(grazm);
 cos_alpha=cos(alpha);
 sin_alpha=sin(alpha);
 cos_beta=cos(beta);
 sin_beta=sin(beta);
 param.band_pass=spacing*beta*dsize_y/L/m;
 param.wavelength_max=param.wavelength_min+param.band_pass;
 p.b2=cos(beta)/L;
 param.range_y=0.99*grazm*msize_y/d;
 param.range_x=0.99*dsize_x/(L+d+s);
 Vector vect(0,d,0);
 Line temp1(vect,-Global_Yaxis);
 param.principal_incoming=*(new Ray(temp1,blazing_wavelength));
 mirror_surface.set_location(Origin);
 temp1.translate(mirror_surface.get_location());
 temp1.rotation(Global_Zaxis,grazm+Pi/2);
 a=sin_alpha*sin_alpha*(1./s+1./L)
   +(cos_alpha-cos_beta)*(cos_alpha-cos_beta)/L;
 b=4*d*sin_grazm*cos_grazm*sin_alpha*sin_alpha/s;
 c=4*d*d*cos_grazm*cos_grazm*sin_alpha*sin_alpha;
 t=4*d*s*cos_grazm;
 q=-2*sin_grazm*cos_grazm*(d+s);
 u=sin_grazm*sin_grazm*(1./d+1./s);
 v=-2*sin_grazm;
 A=a*q*q+b*q+c*u;
 B=2*a*t*q+b*t+c*v;
 C=a*t*t;
 if((q=B*B-4*A*C)<0){
   cout<<"can't find the proper mirror"<<endl;
   exit(1);
 }
 mirror_surface.set_radius(rad=(-B+sqrt(q))/A/2);
 mirror_surface.set_center(temp1.map(rad));
 mirror_surface.set_orientation();
 mirror_surface.set_local_xaxis(Global_Zaxis);
//cout<<"radius of the mirror**"<<rad<<endl;
 sphere_mirror.surf=&mirror_surface;
 sphere_mirror.incoming=param.principal_incoming;
 sphere_mirror.make_outgoing();

 temp1 = sphere_mirror.get_outgoing();
 grating_blank.set_location(temp1.map(s));
//cout<<"grat loc***";
//prn(grating_blank.get_location()); 
 temp1.translate(grating_blank.get_location());
 temp1.rotation(Global_Zaxis,alpha+Pi/2);
 grating_blank.set_orientation(temp1.get_direction());
//cout<<"grat orient****";
//prn(grating_blank.get_orientation());
 grating_blank.set_local_xaxis(Global_Zaxis);
 T1=sin_grazm*sin_grazm/d-sin_grazm/rad;
 T2=sin_grazm*sin_grazm/s-sin_grazm/rad;
 B01=-2*rad*d*cos_grazm*sin_alpha;
 B10=4*d*s*cos_grazm-(d+s)*rad*2*cos_grazm*sin_grazm;
 B02=0.;
 B20=(d-s)*2*cos_grazm*cos_grazm
   -d*rad*sin_grazm*T2
     +s*rad*sin_grazm*T1;
 B11=2*d*(sin_alpha*sin_grazm-2*cos_alpha*cos_grazm)
   -2*(rad*d/s)*sin_grazm*sin_grazm*sin_alpha
     +2*rad*cos_grazm*sin(alpha+grazm);
 B30=-cos_grazm+(s/d+d/s*cos_grazm)*sin_grazm*sin_grazm;
 B30+=rad*s*T1*cos_grazm*sin_grazm/d;
 B30+=d*rad*T2*cos_grazm*sin_grazm/s;
 B21=d*sin_alpha*sin_grazm*cos_grazm/s;
 B21-=rad*sin(grazm+alpha)*T1;
 B21+=cos_alpha*(1-d*sin_grazm/rad);
 B21+=d*sin(grazm+alpha)/rad;
 B21+=2*cos(grazm+alpha);
 B21-=d*rad*(sin_grazm*sin_grazm*cos_alpha
			 -2*sin_grazm*cos_grazm*sin_alpha)*sin_grazm/s/s;
 B12=-d*rad*(2*sin_grazm*sin_alpha*cos_alpha
			 -cos_grazm*sin_alpha*sin_alpha)*sin_grazm/s/s;
 B12-=2*rad*sin_alpha*sin_alpha*sin_grazm*cos_grazm/s;
 B12+=d*sin_alpha*sin_alpha*cos_grazm/s;
 B03=rad*d*(2*sin_alpha*cos_alpha*cos(grazm-alpha)
			-cos_alpha*sin_alpha*sin_alpha*sin_grazm
			-2*cos_grazm*sin_alpha*cos_alpha*cos_alpha)/s/s;
 H1=-B01/B10;
 H2=-B02/B10-(B20/B10)*H1*H1-(B11/B10)*H1;
 H3=2*B01*B20*H2/B10/B10-B30*H1*H1*H1/B10-B21*H1*H1/B10-B12*H1/B10
   -B03/B10+(B02+B20*H1*H1+B11*H1)*B11/B10/B10;
//cout<<"B01**"<<B01<<"**B10**"<<B10<<endl;
//cout<<"B02**"<<B02<<"**B20**"<<B20<<"**B11**"<<B11<<endl;
//cout<<"H1**"<<H1<<"**H2**"<<H2<<"**H3**"<<H3<<endl;
 C2=sin_alpha*sin_alpha/s/2+sin_grazm*sin_alpha*H1/s
   +sin_grazm*(sin_grazm*(1./d+1./s)-2./rad)*H1*H1/2
     +sin_beta*sin_beta/L/2;
 C3=cos_alpha*sin_alpha*sin_alpha/s/s/2;
 C3+=sin_grazm*sin_alpha*H2/s;
 C3+=(2*sin_grazm*cos_alpha*cos_alpha*sin_alpha
       -cos_grazm*sin_alpha*sin_alpha)*H1/s/s/2;
 C3+=cos_grazm*(sin_grazm*sin_grazm*(1./d/d-1./s/s)
		  -(sin_grazm/rad)*(1./d-1./s))*H1*H1*H1/2;
 C3+=sin_grazm*(sin_grazm*(1./d+1./s)-2./rad)*H1*H2;
 C3+=((sin_grazm*sin_grazm*cos_alpha-2*sin_grazm*cos_grazm*sin_alpha)/s
	+cos_grazm*sin_alpha/rad)*H1*H1/s/2;
 C3-=cos_beta*sin_beta*sin_beta/L/L/2;
 C4=(-1+6*cos_alpha*cos_alpha-5*cos_alpha*cos_alpha*cos_alpha*cos_alpha)
   /s/s/s/8;
 C4+=(-1+6*cos_beta*cos_beta-5*cos_beta*cos_beta*cos_beta*cos_beta)
   /L/L/L/8;
 C4+=sin_grazm*sin_alpha*H3/s;
 C4+=(2*sin_grazm*cos_alpha*sin_alpha-cos_grazm*sin_alpha*sin_alpha)
   *H2/s/s/2;
 C4+=(sin_alpha*sin_grazm*(1-1.5*sin_alpha*sin_alpha)
      -cos_grazm*sin_alpha*sin_alpha*cos_alpha)*H1/s/s/s;
 C4+=(T1+T2)*(H2*H2+2*H1*H3)/2;
 C4+=((sin_grazm*sin_grazm*cos_alpha-2*sin_grazm*cos_grazm*sin_alpha)/s/s
      +cos_grazm*sin_alpha/s/rad)*H1*H2;
 C4+=((sin_alpha*sin_alpha-4*sin_alpha*cos_alpha*sin_grazm*cos_grazm
       +(1-1.5*sin_alpha*sin_alpha)*sin_grazm*sin_grazm)/s
      +0.5*(2*sin_alpha*cos_alpha*cos_grazm
	    -sin_grazm*sin_alpha*sin_alpha)/rad)*H1*H1/s/s/2;
 C4+=1.5*(T1/d-T2/s)*cos_grazm*H1*H1*H2;
 C4-=(cos_grazm*cos_alpha*sin_grazm*sin_grazm/s
      +cos(grazm-alpha)*sin_grazm*sin_grazm/s
      -2*cos_grazm*cos_grazm*sin_grazm*sin_alpha/s
      +cos_grazm*cos_grazm*sin_alpha/rad
      -cos(grazm-alpha)*sin_grazm/rad)*H1*H1*H1/s/s/2;
 C4+=((4*cos_grazm*cos_grazm*T1/d/d-T1*T1/d
       +T1/rad/rad)/8
      +(4*cos_grazm*cos_grazm*T2/s/s-T2*T2/s
	+T2/rad/rad)/8)*H1*H1*H1*H1;
//cout<<"C2**"<<C2<<"**C3**"<<C3<<"**C4**"<<C4<<endl;
 p.b3=-C3/(cos_alpha-cos_beta);
 p.b4=C4/(cos_alpha-cos_beta);
 plane_grating.set_param(p);
 plane_grating.surf=&grating_blank;
//cout<<"b2**"<<p.b2<<"**b3**"<<p.b3<<"**b4**"<<p.b4<<endl;
 plane_grating.incoming=sphere_mirror.outgoing;
 plane_grating.make_outgoing();

 temp1 = plane_grating.get_outgoing(); 
 det_surface.set_location(temp1.map(L));
//cout<<"det loc***";
//prn(det_surface.get_location());
 det_surface.set_orientation(-temp1.get_direction());
//cout<<"det orient***";
//prn(det_surface.get_orientation());
 det_surface.set_local_xaxis(Global_Zaxis);
 plane_detector.surf=&det_surface;
}

void SmRGS::write_params(void)
{
 ofstream outfile;
 RulingParams p;
 p=plane_grating.get_param();
 outfile.open(param.param_file,ios::out);
 if(!outfile) {printf("can't open output file for params.\n"); 
	exit(1);}
 outfile<<"the distance between the source and the grating: "
   <<param.source_distance<<" cm"<<endl;
 outfile<<"the distance between the grating and the detecter: "
   <<param.focal_length<<" cm"<<endl;
 outfile<<"the blazing angle of the grating: "
   <<p.blazing_angle/DGR<<" Degree"<<endl;
 outfile<<"the grazing angle on the mirror: "
   <<more_param.mirror_grazing_angle/DGR<<" Degree"<<endl;
 outfile<<"the grazing angle on the grating: "
   <<more_param.grating_grazing_angle/DGR<<" Degree"<<endl;
 outfile<<"             the resolving power: "
	<<param.principal_resolving_power<<endl;
 outfile<<"                   the band pass: "<<int(param.wavelength_min/A2cm)
   <<"A"<<"-"<<int(param.wavelength_max/A2cm)<<"A"<<endl;
 outfile<<"        the radius of the mirror: "
   <<mirror_surface.get_radius()/100<<" m"<<endl;
 outfile<<"           the ruling parameters: "<<endl;
 outfile<<"                              b2: "<<p.b2<<endl;
 outfile<<"                              b3: "<<p.b3<<endl;
 outfile<<"                              b4: "<<p.b4<<endl;
 outfile.close();
}


void SmRGS::setup_bench()
{
  Vnode *anode;
  char e0[] = "mirror";
  char e1[] = "grating";
  char e2[] = "detector";
  
  if((anode=new Node<Mirror>(&sphere_mirror))==NULL)
     mem_error(e0);
  bench.add_node(anode);
  if((anode=new Node<Grating>(&plane_grating))==NULL) 
     mem_error(e1);
  bench.add_node(anode);
  if((anode=new Node<Detector>(&plane_detector))==NULL) 
     mem_error(e2);
  bench.add_node(anode);
}
 
int main()
{
  SmRGS spect;
  spect.do_trace();
  return 0;
}


