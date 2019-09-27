//#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <bits/stdc++.h>
#include <map>
#include "Node.h"
#include "Edge.h"
#include <math.h>
#include <list>


#define PI 3.14159265
using namespace std;

void printNorme(float * arr){
  for(int i=0;i<2;i++){
    if(i!=1){
      cout<<arr[i]<<" , ";
    }
    else{
      cout<<arr[i];
    }
  }
  cout<<endl;
}
float computeDistance(Node a, Node b){
  float dx= a.getX()-b.getX();
  float dy=a.getY()-b.getY();

  return sqrt(dx*dx + dy*dy);
}

void removeSubsets(vector<Edge> v, vector<vector<Edge>>& set){

  auto index = find(begin(set), end(set), v);
   set.erase(index);
  // cout<< "size:"<<set.size()<<endl;
}

//************************************ANGLE COMPATIBILITY RELATED FUNCTIONS *****************************************************************************************************************
//************************findQuadrant returns the position of a edge on the circle ********************************
string findQuadrant(Edge a){
  float ax=a.getXComponentNorm();
  float ay=a.getYComponentNorm();

  if(ax>0 && ay>=0){
      return "first";
  }
  if (ax ==0 && ay>0){
      return "ninenty";
  }
  if(ax<0 && ay>0){
      return "second";
  }
  if(ax<0 && ay==0){
      return "pi";
  }
  if(ax<0 && ay<0){
      return "third";
  }
  if(ax ==0 && ay<0){
      return "2pi/3";
  }
  if(ax>0 && ay<0){
      return "fourth";
  }

  return " ";
}

//**************calculateAngle computes the angle of a a given edge based on its position on the circle *******************************
float calculateAngle(Edge b, string quadrant){
  float bx=b.getXComponentNorm();
  float by=b.getYComponentNorm();

  if(quadrant=="first"){
    return atan(by/bx);
  }
  if(quadrant=="ninenty"){
    return 90.0;
  }
  if(quadrant=="second"){
    return 180.0 + atan(by/bx); // This is less than 180.0 since atan() returns negative value
  }
  if(quadrant=="pi"){
    return 180.0;
  }
  if(quadrant=="third"){
    return 180.0 + atan(by/bx);
  }
  if(quadrant=="2pi/3"){
    return 270.0;
  }
  if(quadrant=="fourth"){
    return 360.0 + atan(by/bx); // This is less than 360.0
  }
  return -1.0;
}
//***********************************calculateOverallAngle returns  the angle difference between two edgeS ********************
float  calculateOverallAngle(Edge a, Edge b){
  string Q1=findQuadrant(a);
  string Q2=findQuadrant(b);

  float A1=calculateAngle(a,Q1);
  float A2=calculateAngle(b,Q2);
  float angle=abs(A1-A2);

  return  min(angle, float(360.0-angle)) ;

}
//****************calculateArcLength reurns the arc length between two points*********************************
float angleCompatibility(Edge a, Edge b){
   float r= 1.0;
   float theta=calculateOverallAngle(a,b);
   float arcLength =r*theta*PI/180.0;
   return arcLength;
}
//****************END OF ANGLE COMPATIBILITY RELATED FUNCTIONS************************************************************************************************************************************

//*******************SCALE COMPATIBILITY RELATED FUNCTIONS****************************************************************************************************************************************

float scaleCompatibility(Edge a, Edge b){
  float lavg=(a.getLength() + b.getLength())/2;
  return 2/(lavg/min(a.getLength(),b.getLength()) + max(a.getLength(),b.getLength())/lavg);
}


float postionCompatibility(Edge a, Edge b){
  float midX=a.getHead().getX() + abs(a.getTail().getX()-a.getHead().getX())/2.0;
  float midY= a.getHead().getY() + abs(a.getTail().getY()-a.getHead().getY())/2.0;
  Node Pm(a.getEdgeName(), midX,midY);

  float midXb=b.getHead().getX() + abs(b.getTail().getX()-b.getHead().getX())/2.0;
  float midYb= b.getHead().getY() + abs(b.getTail().getY()-b.getHead().getY())/2.0;
  Node Qm(b.getEdgeName(), midXb,midYb);

  float lavg=(a.getLength() + b.getLength())/2;

  return lavg/(lavg+ computeDistance(Pm,Qm));
}
Node projectNodeOnEdge(Node p, Edge a){
  float length=a.getLength();
  float l= ((a.getHead().getY()- p.getY()) * (a.getHead().getY() - a.getTail().getY()) -
          (a.getHead().getX()-p.getX()) * (a.getTail().getX()-a.getHead().getX()))/(length*length);
  float xc=a.getHead().getX() + l*(a.getTail().getX()-a.getHead().getX());
  float yc=a.getHead().getY() + l*(a.getTail().getY()-a.getHead().getY());
  Node I(a.getEdgeName(),xc,yc);
  return I;
}

float edgesVisibility(Edge a, Edge b){
  Node I0=projectNodeOnEdge(b.getHead(), a);
  Node I1=projectNodeOnEdge(b.getTail(),a);

  float midXI=I0.getX() + abs(I1.getX()-I0.getX())/2.0;
  float midYI= I0.getY() + abs(I1.getY()-I0.getY())/2.0;
  Node Im(a.getEdgeName(), midXI,midYI);

  float midXP=a.getHead().getX() + abs(a.getTail().getX()-a.getHead().getX())/2.0;
  float midYP= a.getHead().getY() + abs(a.getTail().getY()-a.getHead().getY())/2.0;
  Node midP(a.getEdgeName(), midXP,midYP);
  float value=0.0;

  float result=max(1-2*computeDistance(midP,Im)/computeDistance(I0,I1),value);

  return result;
}

float VisibilityCompatibility(Edge p, Edge q){
    return min(edgesVisibility(p,q), edgesVisibility(q,p));
}

float compatibilityscore(Edge a, Edge b){
  float result= angleCompatibility(a,b)*scaleCompatibility(a,b)*postionCompatibility(a,b)*VisibilityCompatibility(a,b);
  return result;
}

void combinedCompatibility(vector<Edge>list, vector<vector<Edge>>& listofVectors){
  for( int i=0; i<list.size(); i++){
      vector<Edge>tempvector;
      for(int j=i+1; j<list.size(); j++){
         float arclength = angleCompatibility(list[i],list[j]);
         float scale= scaleCompatibility(list[i],list[j]);
         float position= postionCompatibility(list[i],list[j]);
         float visibility=VisibilityCompatibility(list[i],list[j]);
        // cout<<list[i]<< " and "<< list[j]<<": "<<visibility<<" " <<list[i].getLength()<< " " <<list[j].getLength()<<endl;
         if(arclength<0.0008 && scale>0.75 && position > 0.71 && visibility >=.5){
           tempvector.push_back(list[j]);
         }
      }
      if(!tempvector.empty()){
        tempvector.push_back(list[i]);
        listofVectors.push_back(tempvector);
      }
  }
}

void findSubgraphs(vector<vector<Edge>>&subgraphs, vector<vector<Edge>> listofVectors){
  bool contains;
   for(int i=0; i<listofVectors.size(); i++){
     for(int j=0; j<listofVectors.size(); j++){
       if(listofVectors[i].size() <= listofVectors[j].size()  && i!=j){
         for(int k=0; k<listofVectors[i].size(); k++){
           Edge temp=listofVectors[i][k];
           if(find(listofVectors[j].begin(), listofVectors[j].end(), temp) != listofVectors[j].end()) {
             contains=true;
           }else{
             contains=false;
             break;
           }
         }
         if(contains){
           subgraphs.push_back(listofVectors[i]);
           break;
         }
       }
     }
   }
}


int main(int argc, char * argv[]){
//************************READING THE NODE FILE ****************************************************************************
     ifstream myFile("ailines_nodes.txt");
     if(!myFile.is_open()){
       cout<<"Error Opening file"<<endl;
       return 0;
     }
     string NodeIndexname, xCoordinate, yCoordinate;
     float x, y;
     string line;
     map<string, Node>nodesmap;

     while(getline(myFile, line)){
       stringstream ss(line);
       getline(ss, NodeIndexname, ' ');
       getline(ss, xCoordinate, ' ');
       getline(ss, yCoordinate, ' ');
       x=stof(xCoordinate);
       y=stof(yCoordinate);
       Node node(NodeIndexname, x,y);
       nodesmap.insert(pair<string, Node>(node.getIndex(),node));

     }
       myFile.close();

//******************************READING THE EDGE FILE****************************************************
     ifstream myfile2("ailines_edges.txt");
     if (!myfile2.is_open()){
       cout<<"Error opening Edge file"<<endl;
     }
      //map<string, Edge>edgesmap;
     vector<Edge> edgesvector;
     string edgeName, headname, tailname, textline;
     Node head,tail;
     while(getline(myfile2, textline)){
       stringstream ss(textline);
       getline(ss,edgeName, ' ');
       getline(ss,headname, ' ');
       getline(ss, tailname, ' ');
       head=nodesmap.find(headname)->second;
       tail=nodesmap.find(tailname)->second;
       Edge edge(edgeName,head,tail);
       edgesvector.push_back(edge);
     }
     myfile2.close();

//*******************************************END OF PARSER**********************************************************************

//*****************CREATINGA OUTPUT FILE WHICH CONTAINS NORMALIZED VECTOR****************************************************
    ofstream outputFile ("output_point.txt");
    if (!outputFile.is_open()){
       cout << "Unable to to create a file";
      return 0;
    }
    for(Edge x: edgesvector){
       //cout<<x<<" "<<x.getXComponentNorm()<<" "<<x.getYComponentNorm()<<endl;
       outputFile<<x.getXComponentNorm()<<" "<<x.getYComponentNorm()<<" "<<x.getLength()<<" "<<x.getmidpoint().getX()<<" "<<x.getmidpoint().getY()<<endl;
    }

    // for(Edge x: edgesvector){
    //    cout<<x<<" "<<x.getXComponentNorm()<<" " <<x.getYComponentNorm()<<endl;
    //    outputFile<<x.getXComponentNorm()<<" " <<x.getYComponentNorm()<<endl;
    // }

    // for(Edge x: edgesvector){
    //    cout<<x<<" "<<x.getLength()<<endl;
    //    outputFile<<x.getLength()<<endl;
    // }

    outputFile.close();
//************************************ COMPATIBILITIES******************************************************************************

     vector<vector<Edge>>edgesvectorList;
     combinedCompatibility(edgesvector, edgesvectorList);

  //  ******FINDING ALL SUBSETS********************************************
     vector<vector<Edge>>duplicates;
     findSubgraphs(duplicates,edgesvectorList);

//**********************NOW WE NEED TO REMOVE ALL SUBSESTS and FREE UP THE MEMORY*************************************
     for(vector<Edge>subset: duplicates){
       removeSubsets(subset, edgesvectorList);
     }
      duplicates.clear();// TO FREE UP THE MEMORY

      //********************FOR DEBUGGING PURPOSE*********************************************
     // cout<<endl;
     // for(int i=0; i<edgesvectorList.size(); i++){
     //   for(int j=0; j<edgesvectorList[i].size(); j++){
     //      cout<<edgesvectorList[i][j]<<" ";
     //    }
     //    cout<<endl;
     //  }
     //  cout<<endl;

    //******END OF CREATING A SUBGRAPHS*******************************************
 //************************************END OF aNGLE COMPATIBILITY*********************************************************************
 //cout<<"done"<<endl;
    return 0;
}
