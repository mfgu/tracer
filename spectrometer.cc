#include "spectrometer.h"
#include <fstream>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>


double Spectrometer::grasp(double wavelength)
{
  double eff, d;
  param.principal_incoming.set_wavelength(wavelength);
  eff=bench.efficiency(param.principal_incoming);
  d = param.source_distance*param.source_distance;
  eff*=(param.range_x * param.range_y)/d;
  return(eff);
}
 

void Spectrometer::write_grasp()
{
  fstream file;
  int i;
  double eff,wavelength;
  file.open(param.grasp_file,ios::trunc|ios::in|ios::out);
  if(!file) cout<<"can't open the file "<<param.grasp_file<<endl;
  for(wavelength=param.wavelength_min;wavelength<param.wavelength_max;
	wavelength+=0.05*A2cm)
    {
      eff=grasp(wavelength);
      file<<wavelength/A2cm<<"  "<<eff<<endl;
    }  
  file.close();
}


void Spectrometer::write_resolving_power() {
  fstream file;
  int i;
  file.open(param.respower_file,ios::trunc|ios::in|ios::out);
  if(!file) cout<<"can't open the file "<<param.respower_file<<endl;
  double wavelength,res,alpha,bata;
  for(wavelength=param.wavelength_min;wavelength<param.wavelength_max;
	wavelength+=0.05*A2cm)
    {
      res=resolving_power(wavelength);
      file<<wavelength/A2cm<<"  "<<res<<endl;
    }
  file.close();
}

void Spectrometer::write_calibration() {
  fstream file;
  int i;
  double wavelength;
  double x,y;
  file.open(param.calib_file,ios::trunc|ios::in|ios::out);
  if(!file) cout<<"can't open the file "<<param.calib_file<<endl;
  for(wavelength=param.wavelength_min;wavelength<param.wavelength_max;
	wavelength+=0.05*A2cm)
    {
      param.principal_incoming.set_wavelength(wavelength);
      bench.incoming=param.principal_incoming;
      bench.make_outgoing();
      if(!bench.inside) continue;
      bench.get_image(x,y);
      file<<wavelength/A2cm<<"  "<<y<<endl;
    }
  file.close();
}

void Spectrometer::do_trace()
{
  fstream file;
  int n,m,i;
  double wavelength=0.;
  double x,y,x0,y0,x1,y1;
  double ux, uz, ur, sig;
  Vector v, p0, p1;
  char buffer[80];
  const gsl_rng_type *T;
  gsl_rng *r;
  
  read_input();
  optimize();
  setup_bench();
  write_calibration();
  write_grasp();
  write_resolving_power();

  write_params();
  
  //return;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  sig = param.source_size/2.35;
  while(True){
    cout<<"input the wavelength("<<int(param.wavelength_min/A2cm)<<"A-"<<
      int(param.wavelength_max/A2cm)<<"A,input 0 to terminate): ";
    cin>>wavelength;
    sprintf(buffer,"%s.%d",param.image_file,int(wavelength));
    wavelength*=A2cm;
    if(wavelength==0) break;
    if(wavelength>param.wavelength_max||wavelength<param.wavelength_min)
      {
	cout<<"wavelength out of the band pass."<<endl;
	continue;
      }
        
    param.principal_incoming.set_wavelength(wavelength);
    bench.incoming=param.principal_incoming; 
    bench.make_outgoing();
    bench.get_image(x0,y0);

    file.open(buffer,ios::trunc|ios::in|ios::out);
    if(!file) cout<<"can't open the file "<<buffer<<endl;
    cout<<"the grasp at wavelength "<<wavelength/A2cm
		<<"A is "<<grasp(wavelength)<<endl;

    for (n = 0; n < MaxN; n++) {
      ur = gsl_ran_gaussian(r, sig);
      p0.set(ur, param.source_distance, 0.0);
      ux = gsl_ran_flat(r, -0.5*param.range_y, 0.5*param.range_y);
      uz = gsl_ran_flat(r, -0.5*param.range_x, 0.5*param.range_x);
      p1.set(ux, 0.0, uz);
      v = p1 - p0;      
      v.normalize();
      bench.incoming.set_passing_point(p0);
      bench.incoming.set_direction(v);
      bench.make_outgoing();
      if (!bench.inside) continue;
      bench.get_image(x, y);      
      file<<(y-y0)/MICRON<<"  "<<(x-x0)<<endl;
    }
    file.close();
  }
  //  cout<<"done"<<endl;
}


