#include <iostream>
#include <string>
#include <stdio.h>
#include "Node.h"
 using namespace std;
Node::Node(){
  index="";
  x=0.0;
  y=0.0;

}
Node::Node(std::string indexname,float xCoordinate, float yCoordinate){
     index = indexname;
     x = xCoordinate;
     y = yCoordinate;
}

void Node::setX(float xcord){
  x=xcord;
}

void Node::setY(float ycord){
  y=ycord;
}

float Node::getX()const{
  return x;
}

float Node::getY()const{
  return y;
}

bool Node::issame(Node &point)const{
  if(x==point.getX() && y==point.getY()){
    return true;
  }
  return false;
}

void Node::setIndex(std::string indexname){
  index=indexname;
}
std::string Node::getIndex()const{
  return index;
}
std::ostream& operator<<(std::ostream &strm, const Node &a) {
  //return strm << "Node: " << a.getIndex()<<" x: "<<a.getX()<<" y: " <<a.getY();
  return strm <<a.getIndex();
}
Node::~Node(){}
