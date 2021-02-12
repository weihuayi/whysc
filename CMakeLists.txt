cmake_minimum_required(VERSION 3.0)

project(WHYSC)

# 目前已经在 Ubuntu 20.04 系统中测试过该软件包

include_directories(${PROJECT_SOURCE_DIR}/include)
message(STATUS "${PROJECT_SOURCE_DIR}")

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)

# C++ 中调要用科学计算的模块，要在系统中安装下面的 python 包
# sudo apt install python3
# sudo apt install python3-dev
# sudo apt install python3-numpy
# sudo apt install python3-scipy
# sudo apt install python3-matplotlib

find_package(PythonLibs 3)
message(STATUS "${PYTHON_INCLUDE_DIRS}")

# sudo apt install freeglut3-dev
find_package(OpenGL COMPONENTS OpenGL)
find_package(GLUT)

if(NOT ${OPENGL_FOUND})
    message("OPENGL not found")
endif()

# 
set(VTK_DIR "~/.local/vtk")
find_package(VTK)
if(NOT ${VTK_FOUND})
    message("VTK not found")
else()
    include_directories(${VTK_INCLUDE_DIRS})
    message("Found vtk! ${VTK_INCLUDE_DIRS}")
endif()

# 
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
set(CGAL_DIR "~/.local/cgal")
find_package(CGAL)
if(NOT ${CGAL_FOUND})
    message("CGAL not found")
else()
    include_directories(${CGAL_INCLUDE_DIRS})
    message("Found cgal! ${CGAL_INCLUDE_DIRS}")
endif()


add_subdirectory(test)
