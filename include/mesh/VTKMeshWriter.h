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

namespace WHYSC {
namespace Mesh {

template<typename Mesh>
class VTKMeshWriter
{
public:
    typedef typename Mesh::I I;
    typedef typename Mesh::F F;
    typedef typename Mesh::Node_iterator Node_iterator;
    typedef typename Mesh::Cell_iterator Cell_iterator;

public:
    VTKMeshWriter()
    {
        m_mesh = NULL;
    }

    VTKMeshWriter(Mesh * mesh)
    {
        m_mesh = mesh;
    }

    void write()
    {
        auto NN = m_mesh->number_of_nodes();
        auto NC = m_mesh->number_of_cells();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->Allocate(NN);

        auto GD = m_mesh->geo_dimension();

        if(GD == 3)
        {
            for(auto it=m_mesh->node_begin(); it != m_mesh->node_end(); it++) 
            {
              points->InsertNextPoint((*it)[0], (*it)[1], (*it)[2]); // (*it) 括号是必须的
            }
        }
        else if(GD == 2)
        {
            for(auto it=m_mesh->node_begin(); it != m_mesh->node_end(); it++) 
            {
              points->InsertNextPoint((*it)[0], (*it)[1], 0.0);
            }
        }
    
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

        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
        unstructuredGrid->SetPoints(points);
        unstructuredGrid->SetCells(m_mesh->vtk_cell_type(), cellArray);

        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName("test.vtu");
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }

private:
    Mesh *m_mesh;
};

} // end of namespace Mesh

} // end of namespace WHYSC
#endif // end of VTKMeshWriter_h
