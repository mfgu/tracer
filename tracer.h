#ifndef TRACER_H
#define TRACER_H 1
#include <stdlib.h>
#include <iostream>
#include "geom_objects.h"
#include "opt_devices.h"

#define LAST -1

extern void prn(Vector);
extern void mem_error(char *s);


class Vnode{
      public:
	Vnode *next;
	Vnode():next(NULL){};
	virtual void make_outgoing()=0;
	virtual Ray &get_outgoing()=0;
	virtual void set_incoming(Ray &r)=0;
	virtual Boolean is_inside()=0;
        virtual Vector get_local()=0;
	virtual double efficiency()=0;
      };



class List{
 public:
  Vnode *head;
  Vnode *tail;
  List():head(NULL),tail(NULL){};
  List(Vnode *n):head(n)
    {
      tail=head;
      head->next=tail;
    }
  
  void add_node(Vnode *anode,int order=LAST);
};


template <class T> 
class Node:public Vnode{
       public:
         T *dev;
	 Node(T *d):Vnode(){dev = d;};
	 virtual void make_outgoing()
	   {
	    dev->make_outgoing();
	    };
	 virtual Boolean is_inside(){return dev->inside;};
	 virtual Ray &get_outgoing(){return dev->outgoing;};
	 virtual void set_incoming(Ray &r){dev->incoming=r;};
	 virtual Vector get_local(){
	   Vector temp;
	   temp = (dev->surf)->get_local((dev->surf)->get_point());
	   return temp;
	 };
 	 virtual double efficiency(){
	   return dev->efficiency();
	 };
       };



class OptBench:public OptElement,public List{
        public:
	   OptBench():
	     OptElement(),List(){};
	   void setup(Vnode *h){head=h;};
	   virtual double efficiency(Ray in);
	   virtual void make_outgoing();
	   virtual void get_image(double &x, double &y)
		{x=(tail->get_local()).get_x();
		 y=(tail->get_local()).get_y();
		 };
	 };

	       
#endif
























