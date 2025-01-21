cmake_minimum_required(VERSION 3.11.0)

project(CreateCBASegments 
  DESCRIPTION "Create CBA Segment Files - console application"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE headers CONFIGURE_DEPENDS *.h)

add_executable(CreateCBASegments ${source} ${headers})

set_target_properties(CreateCBASegments PROPERTIES FOLDER "ConsoleApps")

target_link_libraries(CreateCBASegments PRIVATE FileUtilsLib SpatialContextLib)

target_include_directories(CreateCBASegments PRIVATE ${CMAKE_SOURCE_DIR}/FileUtilsLib ${CMAKE_SOURCE_DIR}/SpatialContextLib)

