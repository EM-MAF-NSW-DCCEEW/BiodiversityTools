cmake_minimum_required(VERSION 3.11.0)

project(OccupancyCBA 
  DESCRIPTION "Occupancy CBA - console application"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE headers CONFIGURE_DEPENDS *.h)

add_executable(OccupancyCBA ${source} ${headers})

set_target_properties(OccupancyCBA PROPERTIES FOLDER "ConsoleApps")

target_link_libraries(OccupancyCBA PRIVATE FileUtilsLib SpatialContextLib)

target_include_directories(OccupancyCBA PRIVATE ${CMAKE_SOURCE_DIR}/FileUtilsLib ${CMAKE_SOURCE_DIR}/SpatialContextLib)

