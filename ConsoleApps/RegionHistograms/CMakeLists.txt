cmake_minimum_required(VERSION 3.11.0)

project(RegionHistograms 
  DESCRIPTION "Region Histograms - console application"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE headers CONFIGURE_DEPENDS *.h)

add_executable(RegionHistograms ${source} ${headers})

set_target_properties(RegionHistograms PROPERTIES FOLDER "ConsoleApps")

target_link_libraries(RegionHistograms PRIVATE FileUtilsLib)

target_include_directories(RegionHistograms PRIVATE ${CMAKE_SOURCE_DIR}/FileUtilsLib)

