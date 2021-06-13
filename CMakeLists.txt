cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(image_rasterizer)
SET(VTK_DIR <vtk build directory path>)
SET(CMAKE_BUILD_TYPE Debug)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(image_rasterizer MACOSX_BUNDLE image_rasterizer)

if(VTK_LIBRARIES)
  target_link_libraries(image_rasterizer ${VTK_LIBRARIES})
else()
  target_link_libraries(image_rasterizer vtkHybrid)
endif()

