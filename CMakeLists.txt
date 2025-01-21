cmake_minimum_required(VERSION 3.11.0)

project(BiodiversityTools
  DESCRIPTION "Biodiversity tools"
)

# using cmake -DCMAKE_TOOLCHAIN_FILE=D:/Data/Dev/ThirdParty/vcpkg/scripts/buildsystems/vcpkg.cmake ..

# Set C++ standard per project
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Set lib and bin directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# use the FOLDER target property to organize targets into folders
# For where generators support recording folder details for targets
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Add subdirectories
add_subdirectory(FileUtilsLib)
add_subdirectory(SpatialLinksLib)
add_subdirectory(SpatialContextLib)
add_subdirectory(ConsoleApps)

