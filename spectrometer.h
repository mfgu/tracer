#ifndef SPECTROMTER_H
#define SPECTROMTER_H

#include "tracer.h"
#include "const.h"

#define MaxN 10000
#define EfficiencyUnit (1e-5)

struct SpectParams{
  double source_size;
  Ray principal_incoming;
  double principal_wavelength;
  double wavelength_min,wavelength_max,band_pass;
  double total_length;
  double focal_length;
  double source_distance;
  double principal_resolving_power;
  double range_x,range_y;
  double det_resolution,detector_x,detector_y;
  char grasp_file[40],calib_file[40];
  char respower_file[40],image_file[40];
  char param_file[40];
  char dump_file[40];
};

class Spectrometer{
 public:
  OptBench bench;
  SpectParams param;
  Spectrometer(){}; 
  virtual void read_input(void)=0;
  virtual void optimize()=0;
  virtual void write_params(void)=0;
  virtual void setup_bench()=0;
  virtual double resolving_power(double wavelength)=0;
  double grasp(double wavelength);
  void do_trace();
  void write_grasp();
  void write_resolving_power();
  void write_calibration();
 };


extern void prn(Vector);


#endif























