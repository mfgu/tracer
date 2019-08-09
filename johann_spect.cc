#include "tracer.h"
#include <fstream>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct _JohannParam_ {
  double R, d2, wc, hc, lamc;
  double D, rminor;
  double wd, hd;
  double lam;
  int nrays;
  char imfile[FILENL];
} JohannParam;

PARAM params[] = {
  {PDBL, 1, "CrystalRadius", NULL},
  {PDBL, 1, "Crystal2d", NULL},
  {PDBL, 1, "CrystalWidth", NULL},
  {PDBL, 1, "CrystalHeight", NULL},
  {PDBL, 1, "LambdaCenter", NULL},
  {PDBL, 1, "SourceDistance", NULL},
  {PDBL, 1, "TorusMinor", NULL},
  {PDBL, 1, "DetectorWidth", NULL},
  {PDBL, 1, "DetectorHeight", NULL},
  {PDBL, 1, "Lambda", NULL},
  {PINT, 1, "NumRays", NULL},
  {PSTR, 1, "ImageFile", NULL},
  {PNUL, 1, "", NULL}
};

void set_params(JohannParam *p) {
  int i;

  i = 0;
  params[i++].addr = &(p->R);
  params[i++].addr = &(p->d2);
  params[i++].addr = &(p->wc);
  params[i++].addr = &(p->hc);
  params[i++].addr = &(p->lamc);
  params[i++].addr = &(p->D);
  params[i++].addr = &(p->rminor);
  params[i++].addr = &(p->wd);
  params[i++].addr = &(p->hd);
  params[i++].addr = &(p->lam);
  params[i++].addr = &(p->nrays);
  params[i++].addr = p->imfile;
}

int main(int argc, char **argv) {
  JohannParam p;
  Sphere xtal_surf;
  Plane det_surf, tp;
  Mirror xtal;
  Detector det;
  Vector v, s0, v0, v1;  
  Ray r;
  double phi0, phi1, phi, x, y, lamd, phid, theta, thetac;
  int n;
  gsl_rng *rnd;
  FILE *f;
  
  if (argc != 2) {
    printf("Usage: johann_spect <parfile>\n");
    exit(1);
  }

  set_params(&p);
  read_param(argv[1], params);
  
  xtal_surf.set_bound(p.hc, p.wc);
  det_surf.set_bound(p.hd, p.wd);
  
  xtal_surf.set_location(Origin);
  xtal_surf.set_radius(p.R);
  xtal_surf.set_center(Vector(p.R, 0, 0));
  xtal_surf.set_orientation();
  xtal_surf.set_local_xaxis(Global_Zaxis);
  
  thetac = 0.5*Pi - asin(p.lamc/p.d2);  
  Line t(Origin, Global_Xaxis);
  t.rotation(Global_Zaxis, thetac);
  s0 = t.map(p.D);
  t.rotation(Global_Zaxis, -2*thetac);
  v = t.map(p.R*cos(thetac));
  det_surf.set_location(v);
  det_surf.set_orientation(-t.get_direction());
  det_surf.set_local_xaxis(Global_Zaxis);
  
  xtal.surf = &xtal_surf;
  det.surf = &det_surf;
    
  phi1 = p.rminor/p.D;
  phi0 = -phi1;
  phi1 += 1.5*Pi;
  phi0 += 1.5*Pi;

  gsl_rng_env_setup();
  rnd = gsl_rng_alloc(gsl_rng_default);

  theta = 0.5*Pi - asin(p.lam/p.d2);
  
  n = 0;
  f = fopen(p.imfile, "w");
  fprintf(f, "# %12.5E %12.5E\n", p.lam, theta);
  while (1) {
    x = gsl_ran_flat(rnd, -0.5*p.hc, 0.5*p.hc);
    //x = 2.0;
    y = gsl_ran_flat(rnd, -0.5*p.wc, 0.5*p.wc);
    phi = gsl_ran_flat(rnd, phi0, phi1);
    //phi = 1.5*Pi;
    v = Vector(p.R-sqrt(p.R*p.R-x*x-y*y), -y, x);
    tp.set_location(v);
    v0 = Vector(p.R, 0, 0) - v;
    v1 = s0 - v;    
    tp.set_orientation(v0);
    tp.set_local_xaxis(v0^v1);
    v0 = Vector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    v0 = tp.get_global(v0) - v;
    r.set_passing_point(v);
    r.set_direction(-v0);
    r.set_wavelength(p.lam);
    xtal.set_incoming(r);
    xtal.make_outgoing();
    det.set_incoming(xtal.get_outgoing());
    det.make_outgoing();
    v1 = det.surf->get_point();
    r.set_passing_point(v1);
    v0 = Vector(p.R*(1-cos(thetac-theta)), p.R*sin(thetac-theta), 0.0);
    r.set_direction(v0 - v1);
    xtal.set_incoming(r);
    v0.set_x(v0.get_x() - p.R);
    v0.normalize();
    lamd = (r.get_direction()*v0)*p.d2;
    xtal.make_outgoing();
    v0 = xtal.get_outgoing().get_direction();
    phid = atan2(-v0.get_y(), v0.get_z());
    phid += 0.5*Pi;
    v = det.surf->get_local(v1);
    fprintf(f, "%15.8f %15.8f %15.8E %15.8E %15.8E %10.3E %10.3E\n", 
	    v.get_y(), v.get_x(), lamd, phid, phi-1.5*Pi, x, y);
    n++;
    if (n % 50000 == 0) {
      printf("processing ray %d\n", n);
      fflush(f);
    }
    if (n >= p.nrays) break;
  }

  fclose(f);
}
