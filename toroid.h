#ifndef TOROID_H
#define TOROID_H

#include "geom_objects.h"

class Toroid:public Surface{
 protected:
  Vector center;
  double horizontal_radius,vertical_radius;
  
 public:
  Vector get_center(){return(center);};
  double get_horizontal_radius(){return(horizontal_radius);};
  double get_vertical_radius(){return(vertical_radius);};
  void set_center(Vector c){center=c;};
  void set_horizontal_radius(double r){horizontal_radius=r;};
  void set_vertical_radius(double r){vertical_radius=r;};

  Toroid(Toroid &tor):Surface(tor){
    center=tor.get_center();
    horizontal_radius=tor.get_horizontal_radius();
    vertical_radius=tor.get_vertical_radius();
  }
  Toroid(Vector loc,Vector ctr,double r1,double r2):Surface(loc,ctr-loc){
    center=ctr;
    horizontal_radius=r1,vertical_radius=r2;
  }
  Toroid():Surface(){};
  virtual void deriv(double &d1, double &d2, Vector pt);
  virtual Boolean intersect(Line l,Vector &intersection,Vector &normal);
};

  
#endif

