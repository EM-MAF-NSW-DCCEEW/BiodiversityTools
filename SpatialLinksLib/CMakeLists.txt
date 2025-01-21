cmake_minimum_required(VERSION 3.11.0)

project(SpatialLinksLib
  DESCRIPTION "Spatial links library"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE headers CONFIGURE_DEPENDS *.h)

if (MSVC)
	add_compile_options("$<$<CONFIG:RELEASE>:/O2>")
else()
	add_compile_options("-march=native" "-fPIC" "$<$<CONFIG:RELEASE>:-O3>")
endif()

add_library(SpatialLinksLib SHARED ${source} ${headers})

target_link_libraries(SpatialLinksLib PRIVATE FileUtilsLib)

target_include_directories(SpatialLinksLib PRIVATE ${CMAKE_SOURCE_DIR}/FileUtilsLib)

target_compile_definitions(SpatialLinksLib PRIVATE LINKS_EXPORTS=1)
