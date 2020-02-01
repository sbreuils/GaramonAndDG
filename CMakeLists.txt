cmake_minimum_required(VERSION 3.5)

project(ga_view)

# set the module path for "findLib.cmake"
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")

# add a "output" directory in the project binary directory
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# require Eigen
find_package(Eigen3 REQUIRED)
if(${EIGEN3_FOUND})
    message(STATUS "lib EIGEN3 found")
    message(STATUS "  version " ${EIGEN3_VERSION_STRING})
    message(STATUS "  include " ${EIGEN3_INCLUDE_DIR})
endif()

# require lib c2ga
find_package(C2GA)
if(${C2GA_FOUND})
    message(STATUS "lib C2GA found")
    message(STATUS "  include " ${C2GA_INCLUDE_DIRS})
    message(STATUS "  lib     " ${C2GA_LIBRARIES})
else()
    message(STATUS "lib C2GA should be set manually, please uncomment the commented lines in the CMakeLists.txt file")
    #set(C2GA_INCLUDE_DIRS "[path to]/garamon_c2ga/src/")
    #set(C2GA_LIBRARIES "[path to]/garamon_c2ga/build/libc2ga.dylib")

    # for manual setting, check if valid
    if(EXISTS "${C2GA_INCLUDE_DIRS}")
        message(STATUS "include for lib c2ga found: ${C2GA_INCLUDE_DIRS}" )
    else()
        message(FATAL_ERROR "include for lib c2ga invalid or not found: ${C2GA_INCLUDE_DIRS}" )
    endif()

    # for manual setting, check if valid
    if(EXISTS "${C2GA_LIBRARIES}")
        message(STATUS "library c2ga found: ${CC3GA_LIBRARIES}" )
    else()
        message(FATAL_ERROR "library cc3ga not found: ${CC3GA_LIBRARIES}" )
    endif()
endif()

# compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -std=c++14")
#set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")

# files to compile
file(GLOB_RECURSE source_files src/*.cpp src/*.hpp)
message(STATUS "  files " ${source_files})

# compile
add_executable(ga_view ${source_files})

#add_executable(ga_view src/Directory.cpp src/Directory.hpp src/Entry.hpp src/Entry.cpp src/main.cpp )


# includes
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${C2GA_INCLUDE_DIRS})

# link
#target_link_libraries(ga_view ${C3GA_LIBRARIES} )


