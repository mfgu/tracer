#include "util.h"
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>

const Data Table[]={
  {78.50,3.67,0.4},
  {78.51,4.85,0.5},
  {78.22,6.50,0.6},
  {77.08,8.84,0.7},
  {74.44,11.18,0.8},

  {73.41,5.54,1.3},
  {74.79,7.66,1.5},
  {75.08,8.60,1.7},
  {75.29,9.68,1.8},
  {75.41,12.43,2.1},

  {75.26,14.18,2.3},
  {74.86,16.23,2.5},
  {74.12,18.65,2.7},
  {70.60,24.85,3.4},
  {62.93,30.25,4.2},

  {55.57,31.42,4.7},
  {36.18,33.16,5.4},
  {45.05,10.87,6.1},
  {51.19,13.15,7.1},
  {52.64,15.68,8.3},

  {52.38,18.80,9.9},
  {52.00,19.85,10.4},
  {50.23,22.80,11.9},
  {49.63,23.32,12.3},
  {47.71,24.80,13.3},
  
  {45.57,26.00,14.6},
  {43.10,27.19,16.0},
  {40.06,27.92,17.6},
  {38.79,28.15,18.3},
  {36.73,28.20,19.5},

  {33.60,27.86,21.6},
  {31.10,27.60,23.6},
  {30.28,27.47,24.2},
  {26.69,26.29,27.4},
  {23.22,24.57,31.6},

  {18.13,21.05,39.8},
  {15.91,18.58,44.8},
  {13.50,13.46,56.3},
  {13.01,10.40,64.4},
  {13.08,9.07,67.6},

  {16.27,5.51,82.1},
  {18.48,5.73,93.4},
  {20.33,7.56,108.8},
  {20.29,8.22,114.3},
  {20.03,11.18,135.5},
};

int search_param(char *pname, PARAM *params) {
  int i;
  
  i = 0;
  while (1) {
    if (params[i].type == PNUL) break;
    if (strcasecmp(pname, params[i].name) == 0) return i;
    i++;
  }
  return -1;
}

void dump_param(char *fn, PARAM *params) {
  FILE *f;
  int i, j;
  
  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    exit(1);
  }

  i = 0;
  while (1) {
    if (params[i].type == PNUL) break;
    fprintf(f, "%s", params[i].name);
    for (j = 0; j < params[i].nelem; j++) {
      switch (params[i].type) {
      case PINT:
	fprintf(f, "\t%d", *(((int *) params[i].addr) + j));
	break;
      case PDBL:
	fprintf(f, "\t%g", *(((double *) params[i].addr) + j));
	break;
      case PSTR:
	fprintf(f, "\t%s", ((char *) params[i].addr) + j*FILENL);
	break;
      default:
	break;
      }
    }
    fprintf(f, "\n");      
    i++;
  }
  fclose(f);
}

void read_param(char *fn, PARAM *params) {
  FILE *f;
  char pname[PARAMNL];
  int *nset;
  int i, j, n, np;

  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open parameter file %s\n", fn);
    exit(1);
  }

  np = 0;
  while (params[np].type != PNUL) np++;
  nset = (int *) malloc(sizeof(int)*np);
  for (i = 0; i < np; i++) nset[i] = 0;
  while (1) {
    n = fscanf(f, "%s", pname);
    if (n == EOF) break;
    if (n != 1) goto NEXTLINE;    
    if (!(isalpha(pname[0]))) goto NEXTLINE;
    i = search_param(pname, params);
    if (i < 0) goto NEXTLINE;
    nset[i] = 0;
    for (j = 0; j < params[i].nelem; j++) {
      switch (params[i].type) {
      case PINT:	
	n = fscanf(f, "%d", ((int *) params[i].addr)+j);
	break;
      case PDBL:
	n = fscanf(f, "%lf", ((double *) params[i].addr)+j);
	break;
      case PSTR:
	n = fscanf(f, "%s", ((char *) params[i].addr) + j*FILENL);
	break;
      default:
	n = -1;
	break;
      }
      if (n == 1) nset[i]++;
    }
  NEXTLINE:
    while (1) {
      j = fgetc(f);
      if (j == '\n' || j == EOF) break;
    }
  }
  fclose(f);

  for (i = 0; i < np; i++) {
    if (nset[i] != params[i].nelem) {
      printf("Parameter %s is not set properly, expecting %d values, %d given\n",
	     params[i].name, params[i].nelem, nset[i]);
      free(nset);
      exit(1);
    }
  }
  free(nset);
}

void f1f2(double w,double &f1,double &f2)
{
  int i;
  w=w/A2cm;
  for(i=0;i<TableSize&&Table[i].wavelength<w;i++);
  f1=Table[i].f1+(Table[i].f1-Table[i-1].f1)*(w-Table[i].wavelength)\
    /(Table[i].wavelength-Table[i-1].wavelength);
  f2=Table[i].f2+(Table[i].f2-Table[i-1].f2)*(w-Table[i].wavelength)\
    /(Table[i].wavelength-Table[i-1].wavelength);
}


double reflectr(double x,double y)
{
  double t1,t2;
  t1=sqrt((x*x-1)*(x*x-1)+y*y);
  t2=t1-(x*x-1);
  t1=sqrt(t1+(x*x-1));
  double t=((1.414*x-t1)*(1.414*x-t1)+t2)/((1.414*x+t1)*(1.414*x+t1)+t2);
  return(sqrt(t));
}
