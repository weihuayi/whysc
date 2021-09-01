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
    VTKMeshWriter()
    {
        m_points = vtkSmartPointer<vtkPoints>::New();
        m_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        m_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    }

    /*
    VTKMeshWriter(Mesh * mesh)
    {
        m_mesh = mesh;
        m_points = vtkSmartPointer<vtkPoints>::New();
        m_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        m_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    }

    VTKMeshWriter(std::shared_ptr<Mesh> mesh)
    {
        m_mesh = &(*mesh);
        m_points = vtkSmartPointer<vtkPoints>::New();
        m_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        m_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    }
    */

    template<typename Mesh>
    void set_points(Mesh & mesh)
    {
        auto NN = mesh.number_of_nodes();
        m_points->Allocate(NN);

        auto GD = mesh.geo_dimension();

        if(GD == 3)
        {
            for(auto it=mesh.node_begin(); it != mesh.node_end(); it++) 
            {
              m_points->InsertNextPoint((*it)[0], (*it)[1], (*it)[2]); // (*it) 括号是必须的
            }
        }
        else if(GD == 2)
        {
            for(auto it=mesh.node_begin(); it != mesh.node_end(); it++) 
            {
              m_points->InsertNextPoint((*it)[0], (*it)[1], 0.0);
            }
        }
        m_ugrid->SetPoints(m_points);
    }


    template<typename Mesh>
    void set_cells(Mesh & mesh)
    {
        auto NC = mesh.number_of_cells();
        auto nn = mesh.number_of_nodes_of_each_cell();
        auto cellArray = vtkSmartPointer<vtkCellArray>::New();
        cellArray->AllocateExact(NC, NC*nn); 
        for(auto & cell : mesh.cells())
        {
            cellArray->InsertNextCell(nn);
            for(int i = 0; i < nn; i++)
            {
                cellArray->InsertCellPoint(cell[mesh.vtk_write_cell_index()[i]]);
            }
        }
        //Toplogy top;
        //mesh.cell_to_node(top);
        m_ugrid->SetCells(mesh.vtk_cell_type(), cellArray);
    }

    void set_point_data(std::vector<int> & data, int ncomponents, const std::string name)
    {
        int n = data.size()/ncomponents;
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
        int n = data.size()/ncomponents;
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
        int n = data.size()/ncomponents;
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
        int n = data.size()/ncomponents;
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
    vtkSmartPointer<vtkPoints> m_points;
    vtkSmartPointer<vtkUnstructuredGrid> m_ugrid;
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> m_writer;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of VTKMeshWriter_h
