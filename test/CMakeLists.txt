
configure_file( TestMacro.h.in ${CMAKE_SOURCE_DIR}/test/TestMacro.h)


add_executable(test_TriangleMesh test_TriangleMesh.cpp)
add_test(TriangleMesh test_TriangleMesh)

add_executable(test_ScaledMonomialSpace2d test_ScaledMonomialSpace2d.cpp)
add_test(ScaledMonomialSpace2d test_ScaledMonomialSpace2d)

add_executable(test_InitTop test_InitTop.cpp)
add_test(InitTop test_InitTop)

