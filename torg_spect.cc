#include "spectrometer.h"
#include "toroid.h"
#include <fstream>
#include <stdio.h>

struct ToroidRGSParam{
  double grating_grazing_angle;
};

class ToroidRGS:public Spectrometer{
 public:
  Toroid grating_blank;
  Plane det_surface;
  ToroidRGSParam more_param;
  Grating toroid_grating;
  Detector plane_detector;
  ToroidRGS():Spectrometer(){};
  virtual void read_input(void);
  virtual void optimize();
  virtual void write_params(void);
  virtual void setup_bench();
  virtual double resolving_power(double w);
}; 



double ToroidRGS::resolving_power(double w)
{
  double alpha, beta, delta_alpha, delta_beta, delta_w;
  RulingParams p;
  p=toroid_grating.get_param();
  alpha=more_param.grating_grazing_angle-p.blazing_angle;
  beta=more_param.grating_grazing_angle+p.blazing_angle;
  delta_alpha=param.source_size/param.source_distance;
  delta_beta=param.det_resolution/param.focal_length;
  delta_w=p.spacing_at_center/my_abs(p.m);
  delta_w*=sqrt(sin(alpha)*sin(alpha)*delta_alpha*delta_alpha
	+sin(beta)*sin(beta)*delta_beta*delta_beta);
  return(w/delta_w);
}

void ToroidRGS::read_input(void)
{
  ifstream infile;
  char dummy[80];
  double x,y,rad2;
  RulingParams p;
  infile.open("torg.in",ios::in);
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
  grating_blank.set_bound(x,y);
  infile.getline(dummy,sizeof(dummy));
  (infile>>x>>y).getline(dummy,sizeof(dummy));
  det_surface.set_bound(x,y);

  infile.getline(dummy,sizeof(dummy));
  (infile>>p.tan>>rad2).getline(dummy,sizeof(dummy));
  p.tan=tan(p.tan*DGR);
  rad2*=100;
  grating_blank.set_vertical_radius(rad2);
 
  infile.getline(dummy,sizeof(dummy));
  (infile>>param.grasp_file
   >>param.calib_file).getline(dummy,sizeof(dummy));
  infile.getline(dummy,sizeof(dummy));
  (infile>>param.respower_file
   >>param.image_file).getline(dummy,sizeof(dummy));
  infile.getline(dummy,sizeof(dummy));
  (infile>>param.param_file
   >>param.dump_file).getline(dummy,sizeof(dummy));

  p.is_reflective=True;
  toroid_grating.set_param(p);
}


void ToroidRGS::optimize()
{
  double critical_graz,inc_ang;
  double graz,blaz,d,L;
  double grazmin,grazmax;
  double bestgraz,bestd;
  double alpha,beta,A,B,C,D,Q,d0,d1;
  double f1,f2;
  double y;
  double eff0,eff=0,besteff=0;
  double besteff1=0, bestd1, bestgraz1, bestt1;
  double mind, maxd;
  double spacing,blazing_wavelength;
  double gsize_x,gsize_y,dsize_x,dsize_y;
  int n,m;
  RulingParams p;
  double sin_alpha,cos_alpha,sin_beta,cos_beta;
  double C20,C30,C40,C02,C12,C22,C04;
  double M20,M30,M40,M02,M12,M22,M04,F12;
  double T1,T2,T3,T4;
  double rad1,rad2;
  FILE *outfile;

  outfile = fopen(param.dump_file, "w");
  p=toroid_grating.get_param();
  grating_blank.get_bound(gsize_x,gsize_y);
  det_surface.get_bound(dsize_x,dsize_y);
  spacing=p.spacing_at_center;
  blazing_wavelength=param.principal_wavelength;
  m=-p.m;

  f1f2(blazing_wavelength,f1,f2);
  critical_graz=sqrt(4.069e-5*2*f1)*blazing_wavelength/12.3984/A2cm;
  y=f2/f1;

  //  cout<<"optimizing...wait..."<<endl;
  
  inc_ang= 1e-4;
  grazmin=my_max(1e-4,sqrt(m*blazing_wavelength/spacing/2.));
  grazmax=4.0*critical_graz;
  if(grazmax<grazmin) grazmax=10.0*grazmin;
  //  cout<<"critical angle "<<critical_graz<<endl;
  besteff1 = 0;
  if (param.total_length < 0) {
    mind = 30.0;
    maxd = -param.total_length;
  } else {
    mind = param.total_length;
    maxd = param.total_length;
  }
  for (param.total_length = mind; 
       param.total_length <= maxd;
       param.total_length += 5.0) {
    besteff = 0;
    eff = 0;
    for(graz=grazmin;graz<grazmax;graz+=inc_ang){
      blaz=(m*blazing_wavelength/spacing)/sin(graz)/2.;
      blaz=asin(blaz);
      alpha=graz-blaz; 
      beta=graz+blaz;
      C = param.principal_resolving_power/(2.0*sin(graz)*sin(blaz));
      D = sin(alpha)*param.source_size;
      A = C*D/param.total_length;
      A = A*A;
      D = sin(beta)*param.det_resolution;
      B = C*D/param.total_length;
      B = B*B;
      d1 = 0.5;
      Q = A*(1-d1)*(1-d1) + B*d1*d1;
      D = d1*d1*(1-d1)*(1-d1);
      //    cout<<A<<","<<B<<","<<C<<","<<D<<","<<Q<<endl;
      if (Q > D) continue;
      d0 = 5.0/param.total_length;
      Q = A*(1-d0)*(1-d0)+B*d0*d0;
      D = d0*d0*(1-d0)*(1-d0);
      if (Q < D) continue;
      while (d1-d0 > 1E-5) {
	d = 0.5*(d0 + d1);
	Q = A*(1-d)*(1-d) + B*d*d;
	D = d*d*(1-d)*(1-d);
	if (Q > D) {
	  d0 = d;
	} else {
	  d1 = d;
	}
      }
      d *= param.total_length;
      eff=sin(alpha)/sin(beta)/d/param.total_length;
      eff*=alpha*gsize_y*dsize_x;
      eff0 = reflectr(graz/critical_graz,y);
      eff *= eff0;
      if(eff>besteff){
	besteff=eff;
	bestgraz=graz;
	bestd=d;
      }
    }
    //  cout<<A<<","<<B<<","<<C<<","<<D<<endl;
    if (besteff == 0) continue;
    if (besteff > besteff1) {
      besteff1 = besteff;
      bestgraz1 = bestgraz;
      bestd1 = bestd;
      bestt1 = param.total_length;
    }
    more_param.grating_grazing_angle=graz=bestgraz;
    param.source_distance=d=bestd;
    eff=besteff;
    p.blazing_angle=blaz=asin((m*blazing_wavelength/spacing)/sin(graz)/2.);
    param.focal_length=L=param.total_length-d;
    alpha=graz-blaz;
    beta=graz+blaz;
    sin_alpha=sin(alpha);
    cos_alpha=cos(alpha);
    sin_beta=sin(beta);
    cos_beta=cos(beta);
    T1=sin_alpha*sin_alpha;
    T2=cos_alpha-cos_beta;
    rad1=(T1+T2*T2)/L+T1/d;
    rad1=(sin_alpha+sin_beta-T2*cos_beta/sin_beta)/rad1;
    fprintf(outfile, 
	    "%10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n",
	    param.total_length, param.source_distance, 
	    param.focal_length, rad1, eff,
	    alpha/DGR, beta/DGR);
  }
  fclose(outfile);

  more_param.grating_grazing_angle=graz=bestgraz1;
  param.source_distance=d=bestd1;
  param.total_length = bestt1;
  eff=besteff1;
  p.blazing_angle=blaz=asin((m*blazing_wavelength/spacing)/sin(graz)/2.);
  param.focal_length=L=param.total_length-d;

  rad2=grating_blank.get_vertical_radius();
  alpha=graz-blaz;
  beta=graz+blaz;
  sin_alpha=sin(alpha);
  cos_alpha=cos(alpha);
  sin_beta=sin(beta);
  cos_beta=cos(beta);
  param.band_pass=spacing*beta*dsize_y/L/m;
  param.wavelength_max=param.wavelength_min+param.band_pass;
  param.range_y=0.99*alpha*gsize_y;
  param.range_x=0.99*dsize_x*d/(L+d);
  if (param.range_x > gsize_x) param.range_x = 0.99*gsize_x;

  Vector v(0,d,0);
  Line temp1(v,-Global_Yaxis);
  param.principal_incoming=*(new Ray(temp1,blazing_wavelength));
  grating_blank.set_location(Origin);
  temp1.translate(grating_blank.get_location());
  temp1.rotation(Global_Zaxis,Pi/2.+alpha);
  grating_blank.set_orientation(temp1.get_direction());
    
  T1=sin_alpha*sin_alpha;
  T2=cos_alpha-cos_beta;
  rad1=(T1+T2*T2)/L+T1/d;
  rad1=(sin_alpha+sin_beta-T2*cos_beta/sin_beta)/rad1;
  grating_blank.set_center(temp1.map(rad1));
  grating_blank.set_horizontal_radius(rad1);
  M20=cos_beta*(1./sin_beta/rad1/2-1./L);
  T1=sin_alpha*sin_alpha/d-sin_alpha/rad1;
  T2=sin_beta*sin_beta/L-sin_beta/rad1;
  T3=1./d-sin_alpha/rad1;
  T4=1./L-sin_beta/rad1;
  C20=T1/2.+T2/2;
  C30=T1*cos_alpha/d/2;
  C30-=T2*cos_beta/L/2;
  C40=(T1*4.*cos_alpha*cos_alpha/d/d-T1*T1/d+T3/rad1/rad1)/8.;
  C40+=(T2*4.*cos_beta*cos_beta/L/L-T2*T2/L+T4/rad1/rad1)/8.;
  T2=cos_alpha-cos_beta;
  M30=-C30/T2;
  M40=-C40/T2;
  p.b2=rad1*M20+p.tan/2.;
  if(rad2 <=0.){
    rad2=d*L*(sin_alpha+sin_beta+p.tan)/(d+L);
    grating_blank.set_vertical_radius(rad2);
  }
  if(rad2>rad1){
    rad2=rad1;
    grating_blank.set_vertical_radius(rad2);
  }
  T3=1./d-sin_alpha/rad2;
  T4=1./L-sin_beta/rad2;
  C02=T3+T4;
  C12=T3*cos_alpha/d/2-T4*cos_beta/L/2;
  M12=-p.b2*p.tan/rad1/rad2;
  F12=C12+T2*M12;
  //    cout<<"C02:"<<C02<<"	"<<"C12:"<<C12<<endl;

  p.b3=M30*rad1*rad1+p.b2*p.tan;
  p.b4=M40*rad1*rad1*rad1+1.5*p.b3*p.tan
    -0.25*p.b2*p.tan*p.tan+0.125*p.tan;    
  p.b2=-p.b2/rad1;
  p.b3=p.b3/rad1/rad1;
  p.b4=-p.b4/rad1/rad1/rad1;

  toroid_grating.set_param(p);
  grating_blank.set_local_xaxis(Global_Zaxis);
  toroid_grating.set_surface(&grating_blank);
  
  toroid_grating.set_incoming(param.principal_incoming);
  toroid_grating.make_outgoing();
  
  temp1 = toroid_grating.get_outgoing();
  det_surface.set_location(temp1.map(L));
  det_surface.set_orientation(-temp1.get_direction());
  det_surface.set_local_xaxis(Global_Zaxis);
  plane_detector.set_surface(&det_surface);
}

void ToroidRGS::write_params(void)
{
 ofstream outfile;
 RulingParams p;
 p=toroid_grating.get_param();
 outfile.open(param.param_file,ios::out);
 if(!outfile) {printf("can't open output file for params.\n"); 
	exit(1);}

 outfile<<"the resolving power: "
   <<param.principal_resolving_power<<endl;
 outfile<<"the total length of the spectrometer: "
   <<param.total_length<<" cm"<<endl;
 outfile<<"the distance between the source and the grating: "
   <<param.source_distance<<" cm"<<endl;
 outfile<<"the distance between the grating and the detecter: "
   <<param.focal_length<<" cm"<<endl;
 outfile<<"the blazing angle of the grating: "
   <<p.blazing_angle/DGR<<" Degree"<<endl;
 outfile<<"the grazing angle on the grating: "
   <<more_param.grating_grazing_angle/DGR<<" Degree"<<endl;
 outfile<<"             the resolving power: "
	<<param.principal_resolving_power<<endl;
 outfile<<"                   the band pass: "<<int(param.wavelength_min/A2cm)
   <<"A"<<"-"<<int(param.wavelength_max/A2cm)<<"A"<<endl;
 outfile<<"     the horizontal radius of the grating: "
   <<grating_blank.get_horizontal_radius()/100<<" m"<<endl;
 outfile<<"     the vertical radius of the grating: "
   <<grating_blank.get_vertical_radius()/100<<" m"<<endl;
 outfile<<"           the ruling parameters: "<<endl;
 outfile<<"                              b2: "<<p.b2<<endl;
 outfile<<"                              b3: "<<p.b3<<endl;
 outfile<<"                              b4: "<<p.b4<<endl;
 outfile.close();
}


void ToroidRGS::setup_bench()
{
  Vnode *anode;
  char e0[] = "grating";
  char e1[] = "detector";
  
  if((anode=new Node<Grating>(&toroid_grating))==NULL) 
     mem_error(e0);
  bench.add_node(anode);
  if((anode=new Node<Detector>(&plane_detector))==NULL) 
     mem_error(e1);
  bench.add_node(anode);
}
 
int main()
{
  ToroidRGS spect;
  spect.do_trace();
  return 0;
}


