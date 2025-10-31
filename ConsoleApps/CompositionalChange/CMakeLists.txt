cmake_minimum_required(VERSION 3.11.0)

project(CompositionalChange 
  DESCRIPTION "Compositional Change - console application"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE headers CONFIGURE_DEPENDS *.h)

add_executable(CompositionalChange ${source} ${headers})

set_target_properties(CompositionalChange PROPERTIES FOLDER "ConsoleApps")

target_link_libraries(CompositionalChange PRIVATE FileUtilsLib SpatialContextLib)

target_include_directories(CompositionalChange PRIVATE ${CMAKE_SOURCE_DIR}/FileUtilsLib ${CMAKE_SOURCE_DIR}/SpatialContextLib)

