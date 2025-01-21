cmake_minimum_required(VERSION 3.11.0)

project(SpatialContextLib
  DESCRIPTION "Spatial context library"
  LANGUAGES CXX CUDA
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CUDA REQUIRED)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE headers CONFIGURE_DEPENDS *.h)
file(GLOB_RECURSE cuSource CONFIGURE_DEPENDS *.cu)
file(GLOB_RECURSE cuHeaders CONFIGURE_DEPENDS *.cuh)

add_library(SpatialContextLib SHARED ${source} ${headers} ${cuSource} ${cuHeaders})

#Not sure what's best, add_compile_option() or set_property()
set_property(TARGET SpatialContextLib PROPERTY POSITION_INDEPENDENT_CODE ON)

if(WIN32)
    set_property(TARGET SpatialContextLib PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    set_property(TARGET SpatialContextLib PROPERTY CUDA_SEPARABLE_COMPILATION ON)
elseif(UNIX AND NOT APPLE)
    set_property(TARGET SpatialContextLib PROPERTY CUDA_USE_STATIC_CUDA_RUNTIME OFF)
endif()

target_link_libraries(SpatialContextLib PRIVATE FileUtilsLib ${CUDA_LIBRARIES})

target_include_directories(SpatialContextLib PRIVATE ${CMAKE_SOURCE_DIR}/FileUtilsLib)

target_compile_definitions(SpatialContextLib PRIVATE CUDA_CBA_EXPORTS=1)
