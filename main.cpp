#include <iostream>
#include "gldraw.h"
#include "mpi.h"
#include <fstream>
#include <sstream>
#include <string>
//#include <chrono>
#include <ctime>
//typedef std::chrono::high_resolution_clock Clock;

using namespace std;
#define MY_DATA_TAG 123


int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);

      int rank = 0;
      int size = 0;

      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      char processor_name[MPI_MAX_PROCESSOR_NAME];
      int name_len;
      MPI_Get_processor_name(processor_name, &name_len);
      
      int start = 0;
      int end = 0;
      int numFiles=32;
      start = (numFiles/size)*rank;
      end = (numFiles/size)*(rank+1);
     
      
      //cout << "This is "<<processor_name<<" of " << rank <<" rank, working on    
         
       clock_t begin;
       double duration;
       begin=clock();
       for(int i = start; i < end; i++){
         
	
    
	   char node_filename[1024];
    	   sprintf(node_filename,"airlines_index_nodes.csv");

    	   char edge_filename[1024];
    	   sprintf(edge_filename, "airlines_index_edges_%d.csv", i);

    	   GLDraw my_draw;

	   my_draw.m_sDataNode = node_filename;
           my_draw.m_sDataEdge = edge_filename;

 	   my_draw.bundleEdges();
    
    	}
    	duration=(clock()-begin)/(double)CLOCKS_PER_SEC;
	cout<<"duration: "<<duration<<endl;
      MPI_Finalize();
      return 0;
}
