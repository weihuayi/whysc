
add_executable(hello_world hello_world.cpp)
add_executable(sizeof sizeof.cpp)

add_executable(Matrix_test Matrix_test.cpp)
add_executable(dsolve_test dsolve_test.cpp)
add_executable(ceshi ceshi.cpp)

add_executable(masks masks.cpp)

add_executable(CSRMatrix_test CSRMatrix_test.cpp)
add_executable(laplace_test laplace_test.cpp)

add_executable(laplace2d_test laplace2d_test.cpp)

add_executable(random_test random_test.cpp)

add_executable(house_test house_test.cpp)

add_executable(qr_test qr_test.cpp)

add_executable(MeshGeometry_test MeshGeometry_test.cpp)

add_executable(DataArray_test DataArray_test.cpp)

add_executable(json_test json_test.cpp)

add_executable(txt_to_off_test txt_to_off_test.cpp)

# add_executable(test_TriangleMesh test_TriangleMesh.cpp)
# add_executable(test_Mesh test_Mesh.cpp)

add_executable(matplotlib_test matplotlib_test.cpp)
target_include_directories(matplotlib_test PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(matplotlib_test ${PYTHON_LIBRARIES})

add_executable(heart_plot_test heart_plot_test.cpp)
target_include_directories(heart_plot_test PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(heart_plot_test ${PYTHON_LIBRARIES})

add_executable(quiver quiver.cpp)
target_include_directories(quiver PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(quiver ${PYTHON_LIBRARIES})

if(${OPENGL_FOUND})
    add_subdirectory(opengl)
endif()

if(${VTK_FOUND})

    add_executable(vtk_linear_cell_demo vtk_linear_cell_demo.cpp )
    target_link_libraries(vtk_linear_cell_demo PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS vtk_linear_cell_demo 
        MODULES ${VTK_LIBRARIES}
    )

    #add_executable(TetrahedronMesh_test TetrahedronMesh_test.cpp)
    #target_link_libraries(TetrahedronMesh_test PRIVATE ${VTK_LIBRARIES})
    #vtk_module_autoinit(
    #    TARGETS TetrahedronMesh_test 
    #    MODULES ${VTK_LIBRARIES}
    #)

    add_executable(map_test map_test.cpp)
    target_link_libraries(map_test PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS map_test 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(hexmesh_test hexmesh_test.cpp)
    target_link_libraries(hexmesh_test PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS hexmesh_test 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(pmesh_hexmesh_generation pmesh_hexmesh_generation.cpp)
    target_link_libraries(pmesh_hexmesh_generation PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS pmesh_hexmesh_generation 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(TriRadiusRatioQuality_test TriRadiusRatioQuality_test.cpp)
    target_link_libraries(TriRadiusRatioQuality_test PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS TriRadiusRatioQuality_test 
        MODULES ${VTK_LIBRARIES}
    )

    #add_executable(TetRadiusRatioQuality_test TetRadiusRatioQuality_test.cpp)
    #target_link_libraries(TetRadiusRatioQuality_test PRIVATE ${VTK_LIBRARIES})
    #vtk_module_autoinit(
    #    TARGETS TetRadiusRatioQuality_test 
    #    MODULES ${VTK_LIBRARIES}
    #)

    add_executable(QuadJacobiQuality_test QuadJacobiQuality_test.cpp)
    target_link_libraries(QuadJacobiQuality_test PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS QuadJacobiQuality_test 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(HexJacobiQuality_test HexJacobiQuality_test.cpp)
    target_link_libraries(HexJacobiQuality_test PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS HexJacobiQuality_test 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(halfedge_mesh_test halfedge_mesh_test.cpp)
    target_link_libraries(halfedge_mesh_test PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS halfedge_mesh_test 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(IntervalMesh_test IntervalMesh_test.cpp)
    target_link_libraries(IntervalMesh_test PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS IntervalMesh_test 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(TriangleMesh_test TriangleMesh_test.cpp)
    target_link_libraries(TriangleMesh_test PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS TriangleMesh_test 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(QuadMesh_test QuadMesh_test.cpp)
    target_link_libraries(QuadMesh_test PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS QuadMesh_test 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(pmesh_generation pmesh_generation.cpp)
    target_link_libraries(pmesh_generation PRIVATE ${VTK_LIBRARIES})
    vtk_module_autoinit(
        TARGETS pmesh_generation 
        MODULES ${VTK_LIBRARIES}
    )

endif()

if(${CGAL_FOUND} AND ${VTK_FOUND})
    #add_executable(cgal_mesh_tetmesh cgal_mesh_tetmesh.cpp )
    #target_link_libraries(cgal_mesh_tetmesh PRIVATE ${CGAL_LIBRARIES}
    #    ${VTK_LIBRARIES} ${METIS_LIBRARY})

    add_executable(color_test color_test.cpp )
    message("${CGAL_LIBRARIES}")
    target_link_libraries(color_test PRIVATE ${CGAL_LIBRARIES}
        ${VTK_LIBRARIES} ${METIS_LIBRARY})

    #add_executable(cgal_surface_mesh_generation cgal_surface_mesh_generation.cpp )
    #target_link_libraries(cgal_surface_mesh_generation PRIVATE ${CGAL_LIBRARIES}
    #    ${VTK_LIBRARIES} ${METIS_LIBRARY})

    #add_executable(cgal_volume_mesh_generation cgal_volume_mesh_generation.cpp )
    #target_link_libraries(cgal_volume_mesh_generation PRIVATE ${CGAL_LIBRARIES}
    #    ${VTK_LIBRARIES} ${METIS_LIBRARY})

    add_executable(pmesh_tetmesh_generation pmesh_tetmesh_generation.cpp)
    target_link_libraries(pmesh_tetmesh_generation PRIVATE ${CGAL_LIBRARIES}
        ${VTK_LIBRARIES} ${METIS_LIBRARY})
    vtk_module_autoinit(
        TARGETS pmesh_tetmesh_generation 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(pmesh_trimesh_generation pmesh_trimesh_generation.cpp)
    target_link_libraries(pmesh_trimesh_generation PRIVATE ${CGAL_LIBRARIES}
        ${VTK_LIBRARIES} ${METIS_LIBRARY})
    vtk_module_autoinit(
        TARGETS pmesh_trimesh_generation 
        MODULES ${VTK_LIBRARIES}
    )

    #add_executable(cgal_remesh_test cgal_remesh_test.cpp )
    #message("${CGAL_LIBRARIES}")
    #target_link_libraries(cgal_remesh_test PRIVATE ${CGAL_LIBRARIES}
    #    ${VTK_LIBRARIES} ${METIS_LIBRARY})

    add_executable(PolyhedronModel_test PolyhedronModel_test.cpp )
    message("${CGAL_LIBRARIES}")
    target_link_libraries(PolyhedronModel_test PRIVATE ${CGAL_LIBRARIES}
        ${VTK_LIBRARIES} ${GMSH_LIBRARIES})

    add_executable(Cube_test Cube_test.cpp )
    message("${CGAL_LIBRARIES}")
    target_link_libraries(Cube_test PRIVATE ${CGAL_LIBRARIES}
        ${VTK_LIBRARIES} ${GMSH_LIBRARIES})

    add_executable(CubeAndSphere_test CubeAndSphere_test.cpp )
    message("${CGAL_LIBRARIES}")
    target_link_libraries(CubeAndSphere_test PRIVATE ${CGAL_LIBRARIES}
        ${VTK_LIBRARIES} ${GMSH_LIBRARIES})

    add_executable(RectangleWithHole_test RectangleWithHole_test.cpp )
    message("${CGAL_LIBRARIES}")
    target_link_libraries(RectangleWithHole_test PRIVATE ${CGAL_LIBRARIES}
        ${VTK_LIBRARIES} ${GMSH_LIBRARIES})
endif()

if(${VTK_FOUND} AND ${MPI_FOUND})
    #add_executable(parallel_coloring parallel_coloring.cpp)
    #target_link_libraries(parallel_coloring PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    #vtk_module_autoinit(
    #      TARGETS parallel_coloring 
    #      MODULES ${VTK_LIBRARIES}
    #   )

    add_executable(parallel_coloring_test parallel_coloring_test.cpp)
    target_link_libraries(parallel_coloring_test PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_coloring_test 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_tetmesh_optimization parallel_tetmesh_optimization.cpp)
    target_link_libraries(parallel_tetmesh_optimization PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_tetmesh_optimization 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_quadmesh_optimization parallel_quadmesh_optimization.cpp)
    target_link_libraries(parallel_quadmesh_optimization PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_quadmesh_optimization 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_trimesh_optimization parallel_trimesh_optimization.cpp)
    target_link_libraries(parallel_trimesh_optimization PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_trimesh_optimization 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_hexmesh_optimization parallel_hexmesh_optimization.cpp)
    target_link_libraries(parallel_hexmesh_optimization PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_hexmesh_optimization 
        MODULES ${VTK_LIBRARIES}
    )

    # 下面的是重复的测试
    add_executable(parallel_tetmesh_optimization0 parallel_tetmesh_optimization0.cpp)
    target_link_libraries(parallel_tetmesh_optimization0 PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_tetmesh_optimization0 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_quadmesh_optimization0 parallel_quadmesh_optimization0.cpp)
    target_link_libraries(parallel_quadmesh_optimization0 PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_quadmesh_optimization0 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_trimesh_optimization0 parallel_trimesh_optimization0.cpp)
    target_link_libraries(parallel_trimesh_optimization0 PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_trimesh_optimization0 
        MODULES ${VTK_LIBRARIES}
    )

    add_executable(parallel_hexmesh_optimization0 parallel_hexmesh_optimization0.cpp)
    target_link_libraries(parallel_hexmesh_optimization0 PRIVATE ${VTK_LIBRARIES} ${MPI_LIBRARIES})
    vtk_module_autoinit(
        TARGETS parallel_hexmesh_optimization0 
        MODULES ${VTK_LIBRARIES}
    )

endif()

if(${MPI_FOUND})
    add_executable(mpi_test mpi_test.cpp)
    target_link_libraries(mpi_test PRIVATE ${MPI_LIBRARIES})
endif()

#if(${GMSH_FOUND})
    add_executable(gmsh_test gmsh_test.cpp)
    target_link_libraries(gmsh_test PRIVATE ${GMSH_LIBRARIES})
#endif()

add_subdirectory(cpp)
