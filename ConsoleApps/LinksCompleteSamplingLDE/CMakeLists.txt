cmake_minimum_required(VERSION 3.11.0)

project(LinksCompleteSamplingLDE 
  DESCRIPTION "Spatial Links with complete sampling and linear decay - console application"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE headers CONFIGURE_DEPENDS *.h)

add_executable(LinksCompleteSamplingLDE ${source} ${headers})

set_target_properties(LinksCompleteSamplingLDE PROPERTIES FOLDER "ConsoleApps")

target_link_libraries(LinksCompleteSamplingLDE PRIVATE FileUtilsLib SpatialLinksLib)

target_include_directories(LinksCompleteSamplingLDE PRIVATE ${CMAKE_SOURCE_DIR}/FileUtilsLib ${CMAKE_SOURCE_DIR}/SpatialLinksLib)

