#include "tracer.h"

extern void prn(Vector);



void OptBench::make_outgoing()
{
  Vnode *anode;
  Ray temp;
  anode=head;
  anode->set_incoming(incoming);
  anode->make_outgoing();
  while(anode!=tail){
    temp=anode->get_outgoing();
    anode=anode->next;
    anode->set_incoming(temp);
    anode->make_outgoing();
    if(!(inside=anode->is_inside())) return;
  }
  outgoing=anode->get_outgoing();
}


double OptBench::efficiency(Ray in)
{
  Vnode *node;
  double eff=1.;
  double t;

  incoming=in;
  make_outgoing();
  node=head;
  while(node!=tail){
	eff *= (node->efficiency());
	node=node->next;
     }
  t = tail->efficiency();
  eff*=t;

  return(eff);
}


void List::add_node(Vnode *anode,int order)
{ 
  int i;
  Vnode *temp;
  temp=head;

  if(head==NULL)
    {
      head=anode;
      tail=head;
      head->next=tail;
      return;
    }

  if(order==LAST)
    {
      tail->next=anode;
      anode->next=tail->next;
      tail=anode;
      return;
    } 
  for(i=1;i<order&&temp!=tail;i++) temp=temp->next;
  temp->next=anode;
  anode->next=temp->next;
  if(temp==tail) tail=anode;

}




void mem_error(char *s)
{
  cout<<"can't allocate memory for:"<<s<<endl;
}




















