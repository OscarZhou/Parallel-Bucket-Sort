#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "mpi.h"
using namespace std;




float* create_buckets(int nbuckets, int nitems)
{
  int i;

  int ntotal = nbuckets * nitems;

  // Pointer to an array of more pointers to each bucket
  float* bucket = (float*)calloc(ntotal, sizeof(float*));
  for (i=0; i<ntotal; ++i) bucket[i] = 0;

  // return the address of the array of pointers to float arrays
  return bucket;
}

int compare(const void* x1, const void* x2) {
  const float* f1 = (float*)x1;
  const float* f2 = (float*)x2;
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
  int* nitems = (int*)malloc(nbuckets * sizeof(int));
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


  for (i = 0; i < nbuckets; ++i) {
        sendflag[i] = nitems[i];
  }
  // Don't need the number of items anymore
  free(nitems);

}

int main(int argc,char* argv[])
{
    double t_start, t_start_p1, t_start_p2, t_end;

    int size = atol(argv[1]);  //the size of the array of data numbers
    const float xmin = atof(argv[2]);
    const float xmax = atof(argv[3]);

    MPI::Init(argc,argv);
    int myid = MPI::COMM_WORLD.Get_rank();
    int numproc = MPI::COMM_WORLD.Get_size();
    /**************** Create initial array and another array partitioned ************************/
    int per_length = size / numproc;
    float* data = new float[size];        //the array of the data
    float* per_data = new float[per_length];
    /****************************************/

    int root = 0;
    if (myid == root) {
        /************* Populate the initial array with random number ***************************/
        //std::cout << "SEND " << myid << " : ";
        for(int i=0;i<size;i++)
        {
            data[i]=drand48()*(xmax-xmin-1)+xmin;
            std::cout << data[i] << " ";
        }
        std::cout << std::endl;
        /****************************************/
    }
    /****************** Send each partition to the correct slave process **********************/
    MPI::COMM_WORLD.Scatter(data, per_length, MPI_FLOAT, per_data, per_length, MPI_FLOAT, root);
    /****************************************/


    //std::cout << "RECV " << myid << " : ";
    for (int i = 0; i < per_length; ++i) {
        //std::cout << per_data[i] << " ";
    }
    std::cout << std::endl;


    /******************* Each process use bucket sort in their own data collection *********************/
    int* sendcounts = new int[numproc];
    int* recvcounts = new int[numproc];

    float *buckets = create_buckets(numproc, per_length);
    bucket_sort(per_data, per_length ,xmin,xmax,numproc,buckets, sendcounts);

    std::cout << "AFTER SORT " << myid << " : ";
    for (int i = 0; i < per_length; ++i) {
        //std::cout << per_data[i] << " ";
    }
    std::cout << std::endl;

    //std::cout << "SEND FLAG " << myid << " : ";
    for (int i = 0; i < per_length; ++i) {
        //std::cout << sendflag[i] << " ";
    }
    std::cout << std::endl;
    /****************************************/



    /******************* Reverse the position number of every partition range in diagonal direction among processors *********************/
    MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
    /****************************************/

    int recv_sdispls_length = 0;
    //std::cout << "RECV FLAG " << myid << " : ";
    for(int i=0; i<numproc; i++)
    {
        std::cout << recvcounts[i] << " ";
        recv_sdispls_length += recvcounts[i];
    }
    //std::cout << "recv_sdispls_length "<<recv_sdispls_length << " ";
    std::cout << std::endl;


    /******************* Create and load the arguments to alltoallv *********************/

    int *rdispls = new int[numproc];
    int *sdispls = new int[numproc];

    for(int i=0; i<numproc; i++)
    {
        if(i == 0)
        {
            sdispls[i] = 0;
            rdispls[i] = 0;
        }
        else
        {
            sdispls[i] = sdispls[i-1] + sendcounts[i-1];
            rdispls[i] = rdispls[i-1] + recvcounts[i-1];
        }
        //std::cout << sdispls[i] << " "<<rdispls[i]<<" "<<std::endl;
    }
    /****************************************/


    /******************* Reverse the data of every partition range in diagonal direction among processors*********************/

    float* per_recv_data = new float[recv_sdispls_length];
    MPI_Alltoallv(per_data, sendcounts, sdispls, MPI_FLOAT, per_recv_data, recvcounts, rdispls, MPI_FLOAT, MPI_COMM_WORLD);

    qsort(&per_recv_data[0], recv_sdispls_length, sizeof(float), compare);
    /****************************************/


    std::cout << "AFTER ALLTOALLV " << myid << " : ";
    for (int i = 0; i < recv_sdispls_length; ++i) {
        std::cout << per_recv_data[i] << " ";
    }
    std::cout << std::endl;



    /******************* Gather the displacement number from all processor to master processor*********************/
    int sendcounts_forgather = recv_sdispls_length;
    int* recvcounts_forgather = new int[numproc];
    MPI_Gather(&sendcounts_forgather, 1, MPI_INT,recvcounts_forgather, 1, MPI_INT, root,MPI_COMM_WORLD);

    std::cout << "AFTER Gather " << myid << " : ";
    for (int i = 0; i < numproc; ++i) {
        std::cout << recvcounts_forgather[i] << " ";
    }
    std::cout << std::endl;
    /****************************************/


    /******************* Create and load the arguments to gatherv *********************/

    int *rdispls_v = new int[numproc];

    for(int i=0; i<numproc; i++)
    {
        if(i == 0)
        {
            rdispls_v[i] = 0;
        }
        else
        {
            rdispls_v[i] = rdispls_v[i-1] + recvcounts_forgather[i-1];
        }
        std::cout << rdispls_v[i] << " "<<recvcounts_forgather[i]<<" "<<std::endl;
    }
    std::cout << std::endl;
    /****************************************/
    /*************** Gather the data from each processor to master processor *************************/
    MPI_Gatherv(per_recv_data, recv_sdispls_length, MPI_FLOAT, data,  recvcounts_forgather, rdispls_v,MPI_FLOAT, root, MPI_COMM_WORLD);
    /****************************************/

    std::cout << "FINAL RESULT " << myid << " : ";
    for (int i = 0; i < size; ++i) {
        std::cout << data[i] << " ";
    }
    std::cout << std::endl;

    MPI::Finalize();


    delete[] rdispls_v;
    delete[] sendcounts;
    delete[] recvcounts;
    delete[] rdispls;
    delete[] sdispls;
    delete[] recvcounts_forgather;
    delete[] data;
    delete[] per_data;
}
