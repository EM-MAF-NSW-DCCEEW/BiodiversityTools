cmake_minimum_required(VERSION 3.11.0)

project(LinksRandomSampling 
  DESCRIPTION "Spatial Links with random sampling console application"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE headers CONFIGURE_DEPENDS *.h)

add_executable(LinksRandomSampling ${source} ${headers})

set_target_properties(LinksRandomSampling PROPERTIES FOLDER "ConsoleApps")

target_link_libraries(LinksRandomSampling PRIVATE FileUtilsLib SpatialLinksLib)

target_include_directories(LinksRandomSampling PRIVATE ${CMAKE_SOURCE_DIR}/FileUtilsLib ${CMAKE_SOURCE_DIR}/SpatialLinksLib)

