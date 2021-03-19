#include<stdio.h>
#include<iostream>
#include<list>
#include<vector>

#include"MeshToplogy.h"
using namespace std;

template<typename M>
void improve_mesh_quality(const M & mesh)
{
    typedef typename M::Toplogy Tp;
    typedef typename M::Node Node;
    typedef typename M::Container Container;

    const int GD = mesh.geo_dimension();
    const int NN = mesh.number_of_nodes();
    const int NC = mesh.number_of_cells();

    Tp cell2node;
    mesh.cell_to_node(cell2node);
    Container & nei = cell2node.neighbors();
    Container & loc = cell2node.locations();

    int wei[NN] = {0}; //向量的权重和
    int dir[NN][GD]={0}; //每个点有一个 GD 维向量

    double v[6][GD]={0};
    for(int i=0; i<NC; i++)
    {
        double S = mesh.cell_measure(i); //单元 i 的面积

        std::vector<Node> x; //单元 i 的节点位置
        for(int j = 0; j < 4; j++)
        {
            x.push_back(mesh.node(nei[loc[i]]));
        }

        for(int j = 0; j < 4; j++)
        {
            // 顺序为 i, j, k, m

            auto vij = x[j] - x[(j+1)%4];
            auto vik = x[j] - x[(j+2)%4];
            auto vim = x[j] - x[(j+3)%4];

        }










    }





};
