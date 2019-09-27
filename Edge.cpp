#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>

using namespace std;

#include "Edge.h"

Edge::Edge(string edgename, Node begin, Node end){
  name=edgename;
  if(begin.getX()>end.getX()){
    head=end;
    tail=begin;
  }else{
    head=begin;
    tail=end;
  }
  float dx= pow(tail.getX()-head.getX(),2);
  float dy=pow(tail.getY()-head.getY(),2);
  //cout<<dx<<" "<<dy<< "lol"<<endl;
  length= sqrt(dx+dy);
  //cout<<length<<" kk"<<endl;

  //Node midP(edgename,midXc,midYc);
  //midpoint=&midP;
  XComponentNorm=(tail.getX()-head.getX())/length;
  YComponentNorm=(tail.getY()-head.getY())/length;
//  Node norm(edgename,midXc,midYc);
  //norme=&norm;
  radius=sqrt(pow(XComponentNorm,2)+ pow(YComponentNorm,2));
}
float Edge::getLength()const{
    return length;
}
Node Edge::getHead()const{
  return head;
}
Node Edge:: getTail()const{
  return tail;
}
string Edge::getEdgeName()const{
  return name;
}
std::ostream& operator<<(std::ostream &strm, const Edge &a) {
  return strm <<a.getEdgeName();
  //return strm << "Edge: " << a.getEdgeName();
  //return strm<<" head: "<<a.getHead().getIndex()<<" tail: " <<a.getTail().getIndex()<<" with length: "<<a.getLength();
}
float Edge:: getXComponentNorm()const{
  return XComponentNorm;
}
float Edge::getYComponentNorm()const{
  return YComponentNorm;
}
float Edge::getradius()const{
  return radius;
}
Node Edge::getmidpoint()const{
  float midXc=head.getX() + abs(tail.getX()-head.getX())/2.0;
  float midYc= head.getY() + abs(tail.getY()-head.getY())/2.0;
  Node midpoint(name, midXc, midYc);
  return midpoint;
}

// Node *Edge::getNorme()const{
//   return norme;
// }
// float * Edge::getNorme()const{
//   static float result[2];
//   result[0]=getXComponentNorm();
//   result[1]=getYComponentNorm();
//   return result;
// }

Edge::~Edge(){}
