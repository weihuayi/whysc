#ifndef VTKMeshWriter_h
#define VTKMeshWriter_h

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

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
#include <memory>

namespace WHYSC {
namespace Mesh {

class VTKMeshWriter
{
public:


public:
    VTKMeshWriter()
    {
        m_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        m_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    }

    template<class Mesh>
    void set_mesh(Mesh & mesh)
    {
      m_ugrid->Initialize(); //把网格清空
      set_points(mesh);
      set_cells(mesh);
    }


    template<class Mesh>
    void set_points(Mesh & mesh)
    {
      auto NN = mesh.number_of_nodes();
      auto points = vtkSmartPointer<vtkPoints>::New();
      points->Allocate(NN);
      auto GD = mesh.geo_dimension();

      if(GD == 3)
      {

        int i = 0;
        for(auto it=mesh.node_begin(); it != mesh.node_end(); it++) 
        {
          std::cout << (*it)[0] << ", " << (*it)[1] << ", " << (*it)[2] << std::endl;
          points->InsertNextPoint((*it)[0], (*it)[1], (*it)[2]); // (*it) 括号是必须的
          i++;
        }
      }
      else if(GD == 2)
      {
          for(auto it=mesh.node_begin(); it != mesh.node_end(); it++) 
          {
            points->InsertNextPoint((*it)[0], (*it)[1], 0.0);
          }
      }
      m_ugrid->SetPoints(points);
    }


    template<class Mesh>
    void set_cells(Mesh & mesh)
    {
        auto NC = mesh.number_of_cells();
        auto nn = mesh.number_of_nodes_of_each_cell();

        auto cells = vtkSmartPointer<vtkCellArray>::New();
        cells->AllocateExact(NC, NC*nn); 

        int * idx = mesh.vtk_write_cell_index();
        for(auto & cell : mesh.cells())
        {
            cells->InsertNextCell(nn);
            for(int i = 0; i < nn; i++)
            {
                cells->InsertCellPoint(cell[idx[i]]);
            }
        }
        m_ugrid->SetCells(mesh.vtk_cell_type(), cells);
    }

    void set_point_data(std::vector<int> & data, int ncomponents, const std::string name)
    {
        auto n = data.size()/ncomponents;
        auto vtkdata = vtkSmartPointer<vtkIntArray>::New();
        vtkdata->SetNumberOfComponents(ncomponents);
        vtkdata->SetNumberOfTuples(n);
        vtkdata->SetName(name.c_str());
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < ncomponents; j ++)
                vtkdata->SetComponent(i, j, data[i*ncomponents + j]);
        }
        m_ugrid->GetPointData()->AddArray(vtkdata);
    }

    void set_cell_data(std::vector<int> & data, int ncomponents, const std::string name)
    {
        size_t n = data.size()/ncomponents;
        auto vtkdata = vtkSmartPointer<vtkIntArray>::New();
        vtkdata->SetNumberOfComponents(ncomponents);
        vtkdata->SetNumberOfTuples(n);
        vtkdata->SetName(name.c_str());
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < ncomponents; j ++)
                vtkdata->SetComponent(i, j, data[i*ncomponents + j]);
        }
        m_ugrid->GetCellData()->AddArray(vtkdata);
    }

    void set_point_data(std::vector<double> & data, int ncomponents, const std::string name)
    {
        auto n = data.size()/ncomponents;
        auto vtkdata = vtkSmartPointer<vtkIntArray>::New();
        vtkdata->SetNumberOfComponents(ncomponents);
        vtkdata->SetNumberOfTuples(n);
        vtkdata->SetName(name.c_str());
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < ncomponents; j ++)
                vtkdata->SetComponent(i, j, data[i*ncomponents + j]);
        }
        m_ugrid->GetPointData()->AddArray(vtkdata);
    }

    void set_cell_data(std::vector<double> & data, int ncomponents, const std::string name)
    {
        auto n = data.size()/ncomponents;
        auto vtkdata = vtkSmartPointer<vtkIntArray>::New();
        vtkdata->SetNumberOfComponents(ncomponents);
        vtkdata->SetNumberOfTuples(n);
        vtkdata->SetName(name.c_str());
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < ncomponents; j ++)
                vtkdata->SetComponent(i, j, data[i*ncomponents + j]);
        }
        m_ugrid->GetCellData()->AddArray(vtkdata);
    }


    void write(const std::string & fname)
    {
        m_writer->SetFileName(fname.c_str());
        m_writer->SetInputData(m_ugrid);
        m_writer->Write();
    }

private:
  vtkSmartPointer<vtkUnstructuredGrid> m_ugrid;
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> m_writer;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of VTKMeshWriter_h
