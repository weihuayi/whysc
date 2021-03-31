#ifndef VTKMeshReader_h
#define VTKMeshReader_h

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h> 

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <vtkHexagonalPrism.h>
#include <vtkHexahedron.h>
#include <vtkLine.h>
#include <vtkPentagonalPrism.h>
#include <vtkPixel.h>
#include <vtkPolyLine.h>
#include <vtkPolyVertex.h>
#include <vtkPolygon.h>
#include <vtkPyramid.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkTriangleStrip.h>
#include <vtkVertex.h>
#include <vtkVoxel.h>
#include <vtkWedge.h>

#include "vtkCommonCoreModule.h"
#include "vtkObject.h"
#include "vtkVariant.h"

#include <string>
#include <vector>

namespace WHYSC {
namespace Mesh {

template<typename Mesh>
class VTKMeshReader
{
public:
  
  typedef typename Mesh::Node Node;
  typedef typename Mesh::Cell Cell;

  typedef typename Mesh::I I;
  typedef typename Mesh::F F;
  typedef typename Mesh::Toplogy Toplogy;
  typedef typename Mesh::NodeIterator NodeIterator;
  typedef typename Mesh::CellIterator CellIterator;


public:
    VTKMeshReader()
    {
        m_mesh = NULL;
        m_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        m_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    }

    VTKMeshReader(Mesh * mesh)
    {
        m_mesh = mesh;
        m_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        m_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    }
    

    void read(const std::string & fname)
    {
        m_reader->SetFileName(fname.c_str());
        m_reader->Update();
        m_ugrid = m_reader->GetOutput();
        if(m_mesh != NULL)
        {
          auto NN = m_ugrid->GetNumberOfPoints();
          auto & nodes = m_mesh->nodes();
          nodes.resize(NN);
          for(int i = 0; i < NN; i++)
          {
            m_ugrid->GetPoint(i, nodes[i].data());
          }

          auto NC = m_ugrid->GetNumberOfCells();
          auto & cells = m_mesh->cells();
          cells.resize(NC);
          for(int i = 0; i < NC; i++)
          {
            auto c = m_ugrid->GetCell(i);
            int N = c->GetNumberOfPoints();
            for(int j = 0; j < N; j++)
            {
              cells[i][j] = c->GetPointId(j);
            }
          }
        }
        m_mesh->init_top();
    }

    void get_node_data(const std::string & dname, std::vector<int> & nodedata)
    {
      auto vtkdata = m_ugrid->GetPointData()->GetArray(dname.c_str());
      int N = vtkdata->GetSize();
      nodedata.resize(N);
      for(int i = 0; i < N; i++)
      {
        nodedata[i] = vtkdata->GetComponent(i, 0);
      }
    }

    void get_cell_data(const std::string & dname, std::vector<int> & celldata)
    {
      auto vtkdata = m_ugrid->GetCellData()->GetArray(dname.c_str());
      int N = vtkdata->GetSize();
      celldata.resize(N);
      for(int i = 0; i < N; i++)
      {
        celldata[i] = vtkdata->GetComponent(i, 0);
      }

    }

private:
    Mesh *m_mesh;
    vtkSmartPointer<vtkUnstructuredGrid> m_ugrid;
    vtkSmartPointer<vtkXMLUnstructuredGridReader> m_reader;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of VTKMeshReader_h
