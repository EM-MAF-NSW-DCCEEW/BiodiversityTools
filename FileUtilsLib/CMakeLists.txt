cmake_minimum_required(VERSION 3.11.0)

project(FileUtilsLib
  DESCRIPTION "File utils library for biodiversity tools"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GDAL REQUIRED)

file(GLOB_RECURSE source CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE headers CONFIGURE_DEPENDS *.h)

if (MSVC)
	add_compile_options("$<$<CONFIG:RELEASE>:/O2>")
else()
	add_compile_options("-march=native" "-fPIC" "$<$<CONFIG:RELEASE>:-O3>")
endif()


add_library(FileUtilsLib STATIC ${source} ${headers})

target_link_libraries(FileUtilsLib PRIVATE GDAL::GDAL)

target_include_directories(FileUtilsLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
