
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


add_executable(TriangleMesh_test TriangleMesh_test.cpp)

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

