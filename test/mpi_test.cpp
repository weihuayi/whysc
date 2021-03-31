#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "time.h"

#include "mpi.h"

int main(int argc, char **argv)
{
    int nprocs, myrank;
    int t = 5;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int mydata[10] = {0};
    if(myrank != 1)
    {
      for(int i = 0; i < 10; i++)
      {
        mydata[i] = pow(i, myrank+2);
      }
    }

    if(myrank == 0)
    {
      for(int i = 1; i < 4; i++)
      {
        MPI_Send(mydata, //发送的数据地址
                   10,     //发送数据的长度
                   MPI_INT,//发送数据类型
                   i,      //发送给1进程
                   0, //发送的信息的编号
                   MPI_COMM_WORLD);
      }
    }

    if(myrank != 0)
    {
        MPI_Recv(mydata, //接收的数据存储位置 
               10,     //接收数据的长度
               MPI_INT,//接收数据类型
               0,      //接收1进程的信息
               0,      //接收信息的编号
               MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
    for(int i = 0; i < 10; i++)
    {
      if(myrank == 1)
        std::cout<< myrank << " mydata["<< i << "] = " << mydata[i] <<std::endl;
      if(myrank == 2)
        std::cout<< myrank << " mydata["<< i << "] = " << mydata[i] <<std::endl;
      if(myrank == 3)
        std::cout<< myrank << " mydata["<< i << "] = " << mydata[i] <<std::endl;
    }

     
    std::cout<< "here " << myrank << std::endl;
    MPI_Finalize();
    return 0;
}
