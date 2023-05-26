cmake_minimum_required(VERSION 3.1)
project(GalaxyMaps)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

option(GUI "GUI" ON)

# OpenVolumeMesh
add_subdirectory(ext/OpenVolumeMesh)
include_directories(ext/OpenVolumeMesh/src ${CMAKE_BINARY_DIR}/ext/OpenVolumeMesh/src)

if(GUI)
    # NanoGUI
    set(NANOGUI_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(NANOGUI_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
    set(NANOGUI_INSTALL OFF CACHE BOOL "" FORCE)
    set(NANOGUI_BACKEND OpenGL CACHE STRING "")
    add_subdirectory(ext/nanogui)
    set_property(TARGET glfw glfw_objects nanogui PROPERTY FOLDER dependencies)
    add_definitions(${NANOGUI_EXTRA_DEFS})
    include_directories(${NANOGUI_EXTRA_INCS})
    find_package(OpenGL REQUIRED)
    include_directories(${OPENGL_INCLUDE_DIR})
endif()

# CGAL
find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# Eigen
find_package(Eigen3 REQUIRED CONFIG)

# GMP
find_package(GMP REQUIRED)
find_package(GMPXX REQUIRED)
include_directories(${GMPXX_INCLUDE_DIR})

add_library(GalaxyMaps)
target_sources(GalaxyMaps PRIVATE src/galaxy.cc src/foliation_para.cc src/mesh.cc src/util.cc)
target_link_libraries(GalaxyMaps
    OpenVolumeMesh::OpenVolumeMesh
    CGAL::CGAL Eigen3::Eigen ${GMP_LIBRARIES} ${GMPXX_LIBRARIES}
)
target_include_directories(GalaxyMaps PUBLIC src)

add_executable(galaxy_map)
target_sources(galaxy_map PRIVATE src/main.cc)
target_link_libraries(galaxy_map GalaxyMaps)

if(GUI)
    target_sources(galaxy_map PRIVATE src/canvas.cc)
    target_link_libraries(galaxy_map nanogui ${NANOGUI_EXTRA_LIBS} ${OPENGL_LIBRARIES})
    target_compile_definitions(galaxy_map PRIVATE GUI)
endif()
