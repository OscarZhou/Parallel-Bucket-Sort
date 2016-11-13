#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "mpi.h"
using namespace std;




float* create_buckets(int nbuckets, int nitems)
{
  int i;

  int ntotal = nbuckets * nitems;

  // Pointer to an array of more pointers to each bucket
  float* bucket = calloc(ntotal, sizeof(float*));
  for (i=0; i<ntotal; ++i) bucket[i] = 0;

  // return the address of the array of pointers to float arrays
  return bucket;
}

int compare(const void* x1, const void* x2) {
  const float* f1 = x1;
  const float* f2 = x2;
  float diff = *f1 - *f2;

  return (diff < 0) ? -1 : 1;
}


void bucket_sort(float *data, int ndata, float x1, float x2, int nbuckets,
         float *bucket, int *sendflag)
{
  int i, count;

  // The range covered by one bucket
  float stepsize = (x2 - x1) / nbuckets;

  // The number of items thrown into each bucket. We would expect each
  // bucket to have a similar number of items, but they won't be
  // exactly the same. So we keep track of their numbers here.
  int* nitems = malloc(nbuckets * sizeof(int));
  for (i = 0; i < nbuckets; ++i) nitems[i] = 0;

  // Toss the data items into the correct bucket
  for (i = 0; i < ndata; ++i) {

    // What bucket does this data value belong to?
    int bktno = (int)floor((data[i] - x1) / stepsize);
    int idx = bktno * ndata + nitems[bktno];

    //printf("DATA %d %f %d %d\n", i, data[i], bktno, idx);

    // Put the data value into this bucket
    bucket[idx] = data[i];
    ++nitems[bktno];
  }

  // Sort each bucket using the standard library qsort routine. Note
  // that we need to input the correct number of items in each bucket
  count = 0;
  for (i = 0; i < nbuckets; ++i) {
    if(nitems[i]) {
      qsort(&bucket[i*ndata], nitems[i], sizeof(float), compare);
      memcpy(data,&bucket[i*ndata],nitems[i]*sizeof(float));
      data+=nitems[i];
    }
  }
  //memcpy(data,&bucket[i*ndata],nitems[i]*sizeof(float));

  sendflag = nitems;
  // Don't need the number of items anymore
  free(nitems);

}

int main(int argc,char* argv[])
{
    double t_start, t_start_p1, t_start_p2, t_end;

    const long size = atol(argv[1]);  //the size of the array of data numbers
    const float xmin = atof(argv[2]);
    const float xmax = atof(argv[3]);

    MPI::Init(argc,argv);
    int myid = MPI::COMM_WORLD.Get_rank();
    int numproc = MPI::COMM_WORLD.Get_size();
    /****************************************/
    const long per_length = size / numproc;

    float* data = new float[size];        //the array of the data
    float* per_data = new float[per_length];
    /****************************************/

    int root = 0;
    if (myid == root) {
        /****************************************/
        std::cout << "SEND " << myid << " : ";
        for(int i=0;i<size;i++)
        {
            data[i]=drand48()*(xmax-xmin-1)+xmin;
            std::cout << data[i] << " ";
        }
        std::cout << std::endl;
        /****************************************/
    }
    /****************************************/
    MPI::COMM_WORLD.Scatter(data, per_length, MPI_FLOAT, per_data, per_length, MPI_FLOAT, root);
    /****************************************/


    std::cout << "RECV " << myid << " : ";
    for (int i = 0; i < per_length; ++i) {
        std::cout << per_data[i] << " ";
    }
    std::cout << std::endl;
    /****************************************/
    int* sendflag = new int[numproc];
    int* recvflag = new int[numproc];

    float *buckets = create_buckets(numproc, per_length);
    bucket_sort(per_data, per_length ,xmin,xmax,numproc,buckets, sendflag);
    /****************************************/

    /****************************************/
    MPI::COMM_WORLD.Alltoall(sendflag, numproc, MPI_INT, recvflag, numproc, MPI_INT, MPI_COMM_WORLD);

    /****************************************/

    int recv_sdispls_length = 0;
    for(int i=0; i<numproc; i++)
    {
        recv_sdispls_length += recvflag[i];

    }
    /****************************************/
    float* per_recv_data = new float[per_length];
    per_recv_data = MPI::COMM_WORLD.Alltoallv(per_data, per_length, sendflag, MPI_FLOAT, recv_sdispls_length, recvflag, MPI_FLOAT, MPI_COMM_WORLD);
    qsort(&per_recv_data[0], recv_sdispls_length, sizeof(float), compare);
    /****************************************/



    std::cout << "AFTER ALLTOALLV " << myid << " : ";
    for (int i = 0; i < recv_sdispls_length; ++i) {
        std::cout << per_recv_data[i] << " ";
    }
    std::cout << std::endl;


    /****************************************/
    int* sendflag_forgather[1] = {recv_sdispls_length};
    int* recvflag_forgather = new int[numproc];

    MPI::COMM_WORLD.Gather(sendflag_forgather, 1, MPI_INT, recvflag_forgather, 1, MPI_INT, root);

    /****************************************/
    if (myid == root) {
        /****************************************/
        data = MPI::COMM_WORLD.Gatherv(per_recv_data, recv_sdispls_length, MPI_FLOAT, size, recvflag_forgather, MPI_FLOAT, root, MPI_COMM_WORLD);
        /****************************************/
        std::cout << "FINAL RESULT " << myid << " : ";
        for (int i = 0; i < size; ++i) {
            if(i%10 == 0)
                std::cout << std::endl;
            std::cout << data[i] << " ";

        }
    }


    MPI::Finalize();

    delete[] recvflag_forgather;
    delete[] sendflag;
    delete[] recvflag;
    delete[] data;
    delete[] per_data;
}
