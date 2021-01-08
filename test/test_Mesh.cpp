
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/Triangle.h"
#include "mesh/Mesh.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef GK::Point_2 Node;
typedef GK::Vector_2 Vector;
typedef WHYSC::Mesh::Triangle<int> Cell;
typedef WHYSC::Mesh::Mesh<GK, Vector, Node, Cell> Mesh;

int main(int argc, char **argv)
{
    clock_t start, end;
    start = clock();
    std::vector<double> nodes = {0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0};
    std::vector<int> cells = {1, 2, 0, 3, 0, 2};


    Mesh tri();

    // 插入四个点
    tri.insert(Node{0.0, 0.0});
    tri.insert(Node{1.0, 0.0});
    tri.insert(Node{1.0, 1.0});
    tri.insert(Node{0.0, 1.0});
    
    // 插入两个单元
    tri.insert(Cell{1, 2, 0});
    tri.insert(Cell{3, 0, 2});



    end = clock();
    std::cout<< "运行时间:" << double (end-start)/CLOCKS_PER_SEC << std::endl;

    //tri.print();
    return 0;
}

