
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include<time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/TriangleMesh.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef WHYSC::Mesh::TriangleMesh<GK> Mesh;

int main(int argc, char **argv)
{
    clock_t start, end;
    start = clock();
    std::vector<double> nodes = {0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0};
    std::vector<int> cells = {1, 2, 0, 3, 0, 2};


    Mesh tri(nodes, cells);
    for(int i=0; i<10; i++)
    {
        tri.uniform_refine();
    }
    end = clock();
    std::cout<< "运行时间:" << double (end-start)/CLOCKS_PER_SEC << std::endl;

    //tri.print();
    return 0;
}

