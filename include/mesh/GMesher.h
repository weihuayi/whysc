#ifndef GMesher_h
#define GMesher_h

#include <map>
#include <gmsh.h>
#include <memory>
#include <vector>
#include <math.h>
#include <iostream>

namespace WHYSC {
namespace Mesh {

template<class GK, typename Mesh, typename Model>
class GMesher
{
public:
  typedef typename GK::Point_3 Point;
  typedef typename GK::Vector_3 Vector;
  typedef typename GK::Float F;
  typedef typename GK::Int I;

  typedef typename std::vector<int> Line;
  typedef typename std::vector<int> Face;
  typedef typename std::vector<int> Polyhedron;
public:
  GMesher(std::shared_ptr<Model> model)
  {
    m_model = model;
    m_mesh = std::make_shared<Mesh>();
  }

  void mesher(double lc=1e-1, std::string modlename="T", bool open=false)
  {            
    auto & points = m_model->get_points();
    auto & lines = m_model->get_lines();
    auto & faces = m_model->get_faces();
    auto & volums = m_model->get_volums();

    gmsh::initialize();//初始化
    gmsh::model::add(modlename);//创建 gmsh 中的模型

    for(auto & it : points) //(0, 0, 0) 是坐标, lc 是点的大小, 1 是编号(唯一)
    {
      auto tag = it.first;
      auto p = it.second;
      gmsh::model::geo::addPoint(p[0], p[1], p[2], lc, tag);
    }

    for(auto & it : lines) //添加直线, 不同维数实体的编号独立
    {
      auto tag = it.first;
      auto line = it.second;
      gmsh::model::geo::addLine(line[0], line[1], tag);    
    }

    int NL = lines.size();
    for(auto & it : faces) //添加直线, 不同维数实体的编号独立
    {
      auto tag = it.first;
      auto face = it.second;
      gmsh::model::geo::addCurveLoop(face, tag+NL+1);//TODO 这里要求 m_lines tag 连续.
      gmsh::model::geo::addPlaneSurface({tag+NL+1}, tag);//根据loop得到曲面
    }

    int NF = faces.size();
    for(auto & it : volums)
    {
      auto tag = it.first;
      auto volum = it.second;
      gmsh::model::geo::addSurfaceLoop(volum, tag+NF+1);
      gmsh::model::geo::addVolume({tag+NF+1}, tag);
    }

    gmsh::model::geo::synchronize();//不知道干嘛
    gmsh::model::mesh::generate(3); //生成三维网格

    if(open)
      gmsh::fltk::run();//打开 gmsh

    auto & node = m_mesh->nodes();
    auto & cell = m_mesh->cells();

    std::vector<long unsigned int> nodeTag; //网格中每个节点的 tag
    std::vector<double> nodeCoord; //每个节点的坐标
    std::vector<double> nodeParaCoord; //每个节点的参数坐标, 暂时不知道什么用
    gmsh::model::mesh::getNodes(nodeTag, nodeCoord, nodeParaCoord);

    int NN = nodeTag.size();
    node.resize(NN);
    for(int i = 0; i < NN; i++)
    {
      node[i][0] = nodeCoord[3*i+0];
      node[i][1] = nodeCoord[3*i+1];
      node[i][2] = nodeCoord[3*i+2];
    }

    std::map<int, int> nTag2Nid; //node 的 tag 到 node 编号的 map
    for(int i = 0; i < NN; i++)
      nTag2Nid[nodeTag[i]] = i;

    std::vector<int> nodedim(NN);
    std::vector<int> nodetag(NN);
    for(auto it : volums)
    {
      auto id = it.first;
      std::vector<int> cellType; //单元的类型
      std::vector<std::vector<long unsigned int> > cellTag; //单元的标签
      std::vector<std::vector<long unsigned int> > c2nList; //单元的顶点 tag
      gmsh::model::mesh::getElements(cellType, cellTag, c2nList, 3, id);// 3 维, tag 为 id

      std::vector<long unsigned int> cellTag1 = cellTag[0];
      std::vector<long unsigned int> c2nList1 = c2nList[0];      
      int NC = cell.size();
      int NC0 = cellTag1.size();
      cell.resize(NC + NC0);
      for(int i = 0; i < NC0; i++)
      {
        for(int j = 0; j < 4; j++)
        {
          cell[NC+i][j] = nTag2Nid[c2nList1[4*i+j]]; //单元
          nodetag[cell[NC+i][j]] = id; //出现的点都属于几何中的 id 
          nodedim[cell[NC+i][j]] = 3; //出现的点都是 3 维的
        }
      }
    }

    for(auto it : faces)
    {
      auto id = it.first;
      std::vector<int> cellType; //单元的类型
      std::vector<std::vector<long unsigned int> > cellTag; //单元的标签
      std::vector<std::vector<long unsigned int> > c2nList; //单元的顶点 tag
      gmsh::model::mesh::getElements(cellType, cellTag, c2nList, 2, id);// 2 维, tag 为 id

      std::vector<long unsigned int> cellTag1 = cellTag[0];
      std::vector<long unsigned int> c2nList1 = c2nList[0];      
      int N = c2nList1.size();
      for(int i = 0; i < N; i++)
      {
        int A = nTag2Nid[c2nList1[i]]; //单元
        nodetag[A] = id; //出现的点都属于几何中的 id 
        nodedim[A] = 2; //出现的点都是 3 维的
      }
    }

    for(auto it : lines)
    {
      auto id = it.first;
      std::vector<int> cellType; //单元的类型
      std::vector<std::vector<long unsigned int> > cellTag; //单元的标签
      std::vector<std::vector<long unsigned int> > c2nList; //单元的顶点 tag
      gmsh::model::mesh::getElements(cellType, cellTag, c2nList, 1, id);// 2 维, tag 为 id

      std::vector<long unsigned int> cellTag1 = cellTag[0];
      std::vector<long unsigned int> c2nList1 = c2nList[0];      
      int N = c2nList1.size();
      for(int i = 0; i < N; i++)
      {
        int A = nTag2Nid[c2nList1[i]]; //单元
        nodetag[A] = id; //出现的点都属于几何中的 id 
        nodedim[A] = 1; //出现的点都是 3 维的
      }
    }

    for(auto it : points)
    {
      auto id = it.first;
      std::vector<int> cellType; //单元的类型
      std::vector<std::vector<long unsigned int> > cellTag; //单元的标签
      std::vector<std::vector<long unsigned int> > c2nList; //单元的顶点 tag
      gmsh::model::mesh::getElements(cellType, cellTag, c2nList, 0, id);// 2 维, tag 为 id

      std::vector<long unsigned int> cellTag1 = cellTag[0];
      std::vector<long unsigned int> c2nList1 = c2nList[0];      
      int N = c2nList1.size();
      for(int i = 0; i < N; i++)
      {
        int A = nTag2Nid[c2nList1[i]]; //单元
        nodetag[A] = id; //出现的点都属于几何中的 id 
        nodedim[A] = 0; //出现的点都是 3 维的
      }
    }
    m_mesh->init_top();
    auto & data = m_mesh->nodedata();
    data.gdof = nodedim;
    data.gtag = nodetag;
    gmsh::finalize();//退出 gmsh 环境
  }

  std::shared_ptr<Mesh> get_mesh()
  {
    return m_mesh;
  }
private:
  std::shared_ptr<Model> m_model;
  std::shared_ptr<Mesh> m_mesh;
};

}
}


#endif // end of GMesher_h
