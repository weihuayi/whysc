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

#include <string>

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
          auto & nodes = mesh.nodes();
          nodes.resize(NN);
          for(int i = 0; i < NN; i++)
          {
            m_ugrid->GetPoint(i, nodes[i].data());
          }

          

          auto NC = m_ugrid->GetNumberOfCells();
          auto & cells = mesh.cells();
          cells.reszie(NC);
          for(int i = 0; i < NC; i++)
          {
           
          }
        }
    }

    void get_node_data(const std::string & dname)
    {
    }

    void get_cell_data(const std::string & dname)
    {
    }

private:
    Mesh *m_mesh;
    vtkSmartPointer<vtkUnstructuredGrid> m_ugrid;
    vtkSmartPointer<vtkXMLUnstructuredGridReader> m_reader;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of VTKMeshReader_h
