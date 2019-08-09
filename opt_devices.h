#ifndef OPT_DEVICES_H
#define OPT_DEVICES_H 1


#include "geom_objects.h"
#include "const.h"
#include "util.h"
#include <iostream>
using namespace std;

struct RulingParams{
  double spacing_at_center;
  double b2,b3,b4;
  int m;
  Boolean is_reflective;
  double blazing_angle,tan;
}; 


class Ray:public Line{
 public:
  double wavelength;
 public:
  double get_wavelength(){return(wavelength);};
  void set_wavelength(double wvlen){wavelength=wvlen;};
  
  Ray():Line() {
    wavelength=0;
  };
    
  Ray(Ray &r):Line(r) {
    wavelength=r.get_wavelength();
  };
    
  Ray(Line l,double wvlen):Line(l),wavelength(wvlen) {
  };

  Ray(Vector pt,Vector direct,double wvlen):Line(pt,direct) {
    wavelength=0;
  };	  
};


class OptElement{
       public:
	 Surface *surf;
	 Ray incoming,outgoing;
	 Boolean inside;
       public:	 
	 OptElement():inside(True),surf(NULL){};
	 OptElement(Ray in):incoming(in),inside(True){};
	 OptElement(Ray in, Surface *sur):incoming(in),
		inside(True),surf(sur){};
	 OptElement(Surface *sur):surf(sur),inside(True){};
         OptElement(OptElement &ele)
	   {incoming=ele.incoming;
	    outgoing=ele.outgoing;
	    inside=ele.inside;
	    surf=ele.surf;
	  };
	 void set_incoming(Ray in){incoming=in;};
         void set_outgoing(Ray out){outgoing=out;};
	 Ray & get_incoming(){return(incoming);};
	 Ray & get_outgoing(){return(outgoing);};
	 void set_surface(Surface *sur){surf=sur;};
	 Surface *get_surface(){return surf;};
	 Boolean get_inside(){return(inside);};

	 virtual double efficiency(){return(1.);};
	 virtual double efficiency(Ray in){
	   incoming=in;
	   make_outgoing();
	   return efficiency();
	 };
	 virtual void make_outgoing(){outgoing=incoming;};
       };

class Mirror:public OptElement{
        public:	 
         Mirror():OptElement(){};
	 Mirror(Surface *sur):OptElement(sur){};
	 Mirror(Mirror &a_mirr):OptElement(a_mirr) {};
	 Mirror(Ray in,Surface *sur):OptElement(in,sur){};
	 virtual void make_outgoing();
	 virtual double efficiency();
       };

class Grating:public OptElement{
       private:
	 RulingParams param;
       public:
	 Grating():OptElement(){};
	 Grating(RulingParams p, Surface *sur):OptElement(sur){param=p;};
	 Grating(Grating &g):OptElement(g){};
	 Grating(Ray in,Surface *sur):OptElement(in,sur){};
	 RulingParams get_param(){return param;};
	 void set_param(RulingParams p){param=p;};
	 void deriv_n(double &d1, double &d2, RulingParams p, Vector pt);
	 virtual void get_spacing(Vector p){};
	 virtual void make_outgoing();
	 virtual double efficiency();
       };

class Detector:public OptElement{
       public:
         Detector():OptElement(){};
	 Detector(Ray in):OptElement(in){};
	 Detector(Surface *sur):OptElement(sur){};
	 Detector(Ray in,Surface *sur):OptElement(in,sur){};
	 virtual void make_outgoing();

	      };

#endif


