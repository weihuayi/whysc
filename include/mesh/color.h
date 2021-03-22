#include<stdio.h>
#include <iostream>
#include<list>
#include<vector>

using namespace std;

template<typename Mesh>
void mesh_coloring(Mesh & mesh, int * color)
{
    typedef typename Mesh::Toplogy Tp;

    int GD = mesh.geo_dimension();
    int NN = mesh.number_of_nodes();
    int NE = mesh.number_of_edges();

    int c = 0;
    int data[NN];

    list<int> edge; //没有被删除的边
    list<int> nColored; //没有被染色的点
    for(int i=0; i < NE; i++){edge.push_back(i);}
    for(int i=0; i < NN; i++){nColored.push_back(i);}

    //染色
    while(nColored.size()>0){
        c++;
        bool is_min_data[NN];
        for(int i = 0; i < NN; i++){
            is_min_data[i] = true; //初始所有未染色点都是最小点
        }

        for(list<int>::iterator i = nColored.begin(); i != nColored.end(); i++){
            data[*i] = rand();//产生随机数
        }

        for(list<int>::iterator i = edge.begin(); i != edge.end();){
            int & end0 = mesh.edge(*i)[0];
            int & end1 = mesh.edge(*i)[1];

            //判断*i是否为局部最小data值
            if(color[end0]!=0 || color[end1]!=0){
                edge.erase(i++);
            }
            else{
                if(data[end0] < data[end1]){
                    is_min_data[end1] = false;
                }
                else{
                    is_min_data[end0] = false;
                }
                i++;
            }
        }
        for(list<int>::iterator i = nColored.begin(); i != nColored.end();){
            if(is_min_data[*i])
            {
                color[*i] = c;
                nColored.erase(i++);
            }
            else{i++;}
        }

    }

    //降低颜色数量
    int cMax = c;
    Tp node2node;
    mesh.node_to_node(node2node);
    auto & nei = node2node.neighbors();
    auto & loc = node2node.locations();

    while(c>0){
        for(int i = 0; i < NN; i++){
            list<int> com;
            for(int i=0; i < cMax;i++){com.push_back(i+1);}
            for(int j = loc[i]; j < loc[i+1]; j++){
                com.remove(color[nei[j]]);
            }
            color[i] = *min_element(com.begin(), com.end());
        }

        int num = 0;
        for(int i = 0; i < NN; i++){
            if(color[i]==cMax){
                num++;
            }
        }

        if(num==0){cMax -= 1;}
        c--;
    }
};

