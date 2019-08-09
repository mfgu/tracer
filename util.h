#ifndef UTIL_H
#define UTIL_H 1

#include "const.h"

struct Data{
  double f1;
  double f2;
  double wavelength;
};

typedef struct _PARAM_ {
  int type;
  int nelem;
  char name[PARAMNL];
  void *addr;
} PARAM;

void f1f2(double wavelength, double &f1, double &f2);
double reflectr(double x, double y);
int search_param(char *pname, PARAM *params);
void read_param(char *fn, PARAM *params);
void dump_param(char *fn, PARAM *params);

#endif
