cmake_minimum_required(VERSION 3.11.0)

project(BFT 
  DESCRIPTION "Biodiversity Forecasting Tool (BFT) - console application"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE headers CONFIGURE_DEPENDS *.h)

add_executable(BFT ${source} ${headers})

set_target_properties(BFT PROPERTIES FOLDER "ConsoleApps")

target_link_libraries(BFT PRIVATE FileUtilsLib SpatialContextLib)

target_include_directories(BFT PRIVATE ${CMAKE_SOURCE_DIR}/FileUtilsLib ${CMAKE_SOURCE_DIR}/SpatialContextLib)

