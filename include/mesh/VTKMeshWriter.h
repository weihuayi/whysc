#ifndef VTKMeshWriter_h
#define VTKMeshWriter_h

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

namespace WHYSC {
namespace Mesh {

template<typename Mesh>
class VTKMeshWriter
{
public:
    typedef typename Mesh::I I;
    typedef typename Mesh::F F;
    typedef typename Mesh::Toplogy Toplogy;
    typedef typename Mesh::Node_iterator Node_iterator;
    typedef typename Mesh::Cell_iterator Cell_iterator;


public:
    VTKMeshWriter()
    {
        m_mesh = NULL;
        m_points = vtkSmartPointer<vtkPoints>::New();
        m_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        m_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    }

    VTKMeshWriter(Mesh * mesh)
    {
        m_mesh = mesh;
        m_points = vtkSmartPointer<vtkPoints>::New();
        m_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        m_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    }
    
    void set_points()
    {
        auto NN = m_mesh->number_of_nodes();
        m_points->Allocate(NN);

        auto GD = m_mesh->geo_dimension();

        if(GD == 3)
        {
            for(auto it=m_mesh->node_begin(); it != m_mesh->node_end(); it++) 
            {
              m_points->InsertNextPoint((*it)[0], (*it)[1], (*it)[2]); // (*it) 括号是必须的
            }
        }
        else if(GD == 2)
        {
            for(auto it=m_mesh->node_begin(); it != m_mesh->node_end(); it++) 
            {
              m_points->InsertNextPoint((*it)[0], (*it)[1], 0.0);
            }
        }
        m_ugrid->SetPoints(m_points);
    }

    void set_cells()
    {
        auto NC = m_mesh->number_of_cells();
        auto nn = m_mesh->number_of_nodes_of_each_cell();
        vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
        cellArray->AllocateExact(NC, NC*nn); 
        for(auto it=m_mesh->cell_begin(); it != m_mesh->cell_end(); it++)
        {
            cellArray->InsertNextCell(nn);
            for(auto id:*it)
            {
                cellArray->InsertCellPoint(id);
            }
        }
        //Toplogy top;
        //m_mesh->cell_to_node(top);
        m_ugrid->SetCells(m_mesh->vtk_cell_type(), cellArray);
    }

    void write(const std::sring fname = "test.vtu")
    {
        m_writer->SetFileName(fname.c_str());
        m_writer->SetInputData(unstructuredGrid);
        m_writer->Write();
    }

private:
    Mesh *m_mesh;
    vtkSmartPointer<vtkPoints> m_points;
    vtkSmartPointer<vtkUnstructuredGrid> m_ugrid;
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> m_writer;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of VTKMeshWriter_h
