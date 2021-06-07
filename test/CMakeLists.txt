
configure_file( TestMacro.h.in ${CMAKE_SOURCE_DIR}/test/TestMacro.h)


add_executable(test_TriangleMesh test_TriangleMesh.cpp)
add_test(TriangleMesh test_TriangleMesh)


