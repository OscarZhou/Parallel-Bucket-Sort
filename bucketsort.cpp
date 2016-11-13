#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "mpi.h"
using namespace std;

int main(int argc, char *argv[])
{


  long ndata = atol(argv[1]);

  MPI::Init(argc, argv);
  int numproc = MPI::COMM_WORLD.Get_size();
  int myid    = MPI::COMM_WORLD.Get_rank();


  ndata = ndata / numproc;
  // Fill up the array with data to send to the destination node. Note
  // that the contents of the array will
  float * sendbuf = new float[ndata * numproc];
  float * recvbuf = new float[ndata];

  int root = 0;

  const float xmin = atof(argv[2]);
  const float xmax = atof(argv[3]);

  if (myid == root) {
    std::cout << "SEND " << myid << " : ";
    for (int i = 0; i < ndata * numproc; ++i) {
      sendbuf[i] = drand48()*(xmax-xmin-1)+xmin;;
      std::cout << sendbuf[i] << " ";
    }
    std::cout << std::endl;
  }

  MPI::COMM_WORLD.Scatter(sendbuf, ndata, MPI_FLOAT, recvbuf, ndata, MPI_FLOAT, root);

  std::cout << "RECV " << myid << " : ";
  for (int i = 0; i < ndata; ++i) {
    std::cout << recvbuf[i] << " ";
  }
  std::cout << std::endl;



  MPI::Finalize();

  delete[] sendbuf;
  delete[] recvbuf;
}
