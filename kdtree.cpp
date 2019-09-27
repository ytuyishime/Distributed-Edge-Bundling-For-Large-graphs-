#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <bits/stdc++.h>
#include <math.h>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include "Node.h"
#include "Edge.h"
#include <list>

using namespace std;
const int MAX_DIM = 5;
#define ROW 2101
#define COL 6

typedef std::vector<double> doubleVector;

template <typename Container> // we can make this generic for any container [1]
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};
//function to sort by coordinates
bool sortbyYcoordinate(const vector<double>&a, const vector<double>&b)
{
    return (a[1] < b[1]);
}
bool sortbyZcoordinate(const vector<double>&a, const vector<double>&b)
{
    return (a[2] < b[2]);
}

bool sortbyDcoordinate(const vector<double>&a, const vector<double>&b)
{
    return (a[3] < b[3]);
}
bool sortbyEcoordinate(const vector<double>&a, const vector<double>&b)
{
    return (a[4] < b[4]);
}

vector<int>getCorrectOrder(double (&data)[ROW][COL], vector<int> indices, int depth){
  int low=indices[0];
  int max=indices[0];
  for(int i=0; i<indices.size(); i++){
    if(indices[i]<low){
      low=indices[i];
    }
    if(indices[i]>max){
      max=indices[i];
    }

  }

  std::unordered_map<doubleVector, int, container_hash<doubleVector>>references;

  vector<vector<double>>temp;
  for(int i=low; i<=max;i++){
    vector<double> row;
    for(int j=0; j<COL; j++){
      //cout<<i<<" "<<j<<endl;
      row.push_back(data[i][j]);
      //cout<<data[i][j] << " ";
    }
  //  cout<<endl;
    references.insert({row,i});
    //cout<<i<<endl;
    temp.push_back(row);
  }

  unsigned cd = depth % MAX_DIM;
  if(cd==0){
    sort(temp.begin(), temp.end());
  }else if(cd==1){
    sort(temp.begin(), temp.end(), sortbyYcoordinate);
  }else if(cd==2){
    sort(temp.begin(), temp.end(), sortbyZcoordinate);
  }else if(cd==3){
    sort(temp.begin(), temp.end(), sortbyDcoordinate);
  }else{
    sort(temp.begin(), temp.end(), sortbyEcoordinate);
  }

  int k=low;
  for(int i=0; i<temp.size(); i++){
    for(int j=0; j<temp[i].size(); j++){
      data[k][j]=temp[i][j];
    }
    k++;
  }

 vector<int> orderedlist;
  for(int i=low; i<=max; i++){
    vector<double> row2;
    for(int j=0; j<COL; j++){
      row2.push_back(data[i][j]);
    }


    auto it=references.find(row2);
    //cout<<it->second<<" ";
    orderedlist.push_back(it->second);

  }

  return orderedlist;
 }

 vector<vector<int>> splitVector(vector<int> data){
   vector<vector<int>>vect;
   int size=data.size();
   int halfsize;
   if(size%2==1){
     halfsize=(size+1)/2;
   }else{
     halfsize=size/2;
   }
   vector<int>left;
   vector<int>right;
   for(int i=0; i<halfsize; i++){
     left.push_back(data[i]);
   }
   if (!left.empty()) {
       vect.push_back(left);
   }
   for(int i=halfsize; i<size; i++){
     right.push_back(data[i]);
   }
   if (!right.empty()) {
       vect.push_back(right);
   }
   return vect;
 }
 vector<vector<int>> makeTree(double (&dataset)[ROW][COL], vector<int>data, int depth){
   vector<vector<int>> temp;
   temp.push_back(data);
   vector<vector<int>> tempf;
   vector<vector<int>>vpair;
   int height=0;

   while(height<depth){
      for(int i=0; i<temp.size();i++){
        vector<int>ordered=getCorrectOrder(dataset,temp[i],height);    // I am not sure what really happened
        vpair=splitVector(ordered);
        tempf.push_back(vpair[0]);
        tempf.push_back(vpair[1]);
      }

      temp.clear();
      temp= tempf;

      tempf.clear();
      height++;
  }

return temp;
 }

 void sort(double list[][COL])
 {
      int out,in,temp,temp2;

        for(out=0;out<ROW-1;out++)
        {
         for(in=out+1;in<ROW;in++)
         {
          if(list[in][0]<list[out][0])
          {
           temp=list[out][0];
           temp2=list[out][1];

           list[out][0]=list[in][0];
           list[out][1]=list[in][1];

           list[in][0]=temp;
           list[in][1]=temp2;
          }
         }
        }
 }
int main(int argc, char* argv[]){

            if (argc != 5) {
                cout << "Usage : " << argv[0] << " input_cloud_file input_nodes_file input_edges_file partitionsNumber" << endl;
                return 1;
            }

            char * input_cloud_file = argv[1];
            //char * output_subtree_file = argv[2];
            char * input_nodes_file = argv[2];
            char * input_edges_file = argv[3];
            char * depths = argv[4];

      //************************READING THE NODE FILE ****************************************************************************
           ifstream myFile(input_nodes_file);
           if(!myFile.is_open()){
             cout<<"Error Opening file"<<endl;
             return 0;
           }
           string NodeIndexname, xCoordinate, yCoordinate;
           float xf, yf;
           string aline;
           unordered_map<string, Node>nodesmap;

           while(getline(myFile, aline)){
             stringstream ss(aline);
             getline(ss, NodeIndexname, ' ');
             getline(ss, xCoordinate, ' ');
             getline(ss, yCoordinate, ' ');
             xf=stof(xCoordinate);
             yf=stof(yCoordinate);
             Node node(NodeIndexname, xf,yf);
             nodesmap.insert(pair<string, Node>(node.getIndex(),node));

           }
             myFile.close();

      //******************************READING THE EDGE FILE****************************************************
           ifstream myfile2(input_edges_file);
           if (!myfile2.is_open()){
             cout<<"Error opening Edge file"<<endl;
           }
            //map<string, Edge>edgesmap;
           vector<Edge> edgesvector;
           unordered_map<string, Edge>edgesmap;
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
             edgesmap.insert(pair<string, Edge>(edge.getEdgeName(),edge));
             //edgesvector.push_back(edge);
           }
           myfile2.close();

           for (auto i : edgesmap)
                   cout << i.first << "   " << i.second.getHead()<<" "<<i.second.getTail()
                        << endl;

      ifstream mycloudFile(input_cloud_file);
      if(!mycloudFile.is_open()){
        cout<<"Error Opening file inputs!!"<<endl;
        return 0;
      }

      string xc,yc,zc,uc,vc,idc;
      double x,y,z,u,v,id;
      string line;
      double arr[ROW][COL];
      vector<vector<double>>buffer;
      while(getline(myFile, line)){
        stringstream ss(line);
        getline(ss, xc, ' ');
        getline(ss, yc, ' ');
        getline(ss, zc, ' ');
        getline(ss, uc, ' ');
        getline(ss, vc, ' ');
        getline(ss,idc, ' ');
        x=stod(xc);
        y=stod(yc);
        z=stod(zc);
        u=stod(uc);
        v=stod(vc);
        id=stod(idc);
        vector<double> datarow={x,y,z,u,v,id};
        buffer.push_back(datarow);
      }
      mycloudFile.close();



    sort(buffer.begin(), buffer.end());
    for(int i=0; i<ROW; i++){
      for(int j=0; j<COL; j++){
        arr[i][j]=buffer[i][j];
      }
    }

  vector<int> indexvector;
    for (int i = 0; i <ROW ; i++){
        indexvector.push_back(i);
   }
   stringstream str;
   str << depths;
   int t;
   str >> t;
  vector<vector<int>> partitions=makeTree(arr,indexvector,t);

  ofstream outputFile;
  char filename[1024];
  for(int i=0; i<pow(2,t); i++){
    sprintf(filename, "ailines_edges_%d.cvs",i);
    outputFile.open(filename);
    if (!outputFile) {
			cout << "Cannot open " << filename << endl;
			return 1;
    }else{
      cout << "Write " << filename << endl;
    }
    for (int j=0; j<partitions[i].size(); j++){
      double edgeindex= arr[partitions[i][j]][5];
      auto it=edgesmap.find(to_string(edgeindex));
      Edge tempEdge=it->second;
      outputFile << tempEdge.getEdgeName()<<tempEdge.getHead()<<tempEdge.getTail()<< endl;
    }
   		
		outputFile.close();

  }
//
// //subgraphs_cloud_file
// //  ofstream outputFile (output_subtree_file);
//   if (!outputFile.is_open()){
//      cout << "Unable to create a file";
//     return 0;
//   }else {
//      cout << "Save " << output_subtree_file << endl;
//   }
//     for(int i=0; i<partitions.size(); i++){
//       for(int j=0; j<partitions[i].size(); j++){
//          outputFile<<arr[partitions[i][j]][5]<<" ";
//       }
//       outputFile<<endl;
//     }
//   outputFile.close();

  return 0;
}
