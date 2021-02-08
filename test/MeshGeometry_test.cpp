#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "geometry/Geometry_kernel.h"
#include "mesh/MeshGeometry.h"

typedef WHYSC::Geometry_kernel<double, int> GK;
typedef WHYSC::Mesh::MeshGeometry<GK>  MGeo;
typedef MGeo::Node Node;

int main(int argc, char **argv)
{
    MGeo mgeo(4);
    // 插入四个点
    mgeo.push_back(0.0, 0.0);
    mgeo.push_back(1.0, 0.0);
    mgeo.push_back(1.0, 1.0);
    mgeo.push_back(0.0, 1.0);
    mgeo.print();


    std::vector<Node> nodes;
    mgeo.node(nodes);

    std::cout << "The total number of nodes is " << nodes.size() << std::endl;
    for(int i = 0; i < nodes.size(); i++)
    {
        std::cout << i << ": " << nodes[i] << std::endl;
    }

    for(auto it = mgeo.node_begin(); it != mgeo.node_end(); it++)
    {
        std::cout << *it << std::endl;
    }
    return 0;
}
