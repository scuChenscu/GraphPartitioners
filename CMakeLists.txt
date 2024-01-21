set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_OSX_ARCHITECTURES "arm64")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

cmake_minimum_required(VERSION 3.17)
project(GraphPartitioning)

set(GLOG /opt/homebrew/Cellar/glog/0.6.0)
set(GFLAGS /opt/homebrew/Cellar/gflags/2.2.2)
set(BOOST /opt/homebrew/Cellar/boost/1.83.0)

include_directories(${GLOG}/include/ ${GFLAGS}/include/ ${BOOST}/include/)
link_directories(${GLOG}/lib/ ${GFLAGS}/lib/)
add_library(SOURCE_FILES
        converter/conversions.cpp
        converter/shuffler.cpp
        baseGraph/base_graph.cpp
        partitioner/partitioner.cpp
        partitioner/edgePartitioner.cpp
        partitioner/vertexPartitioner.cpp
        utils/graph.cpp
        algorithms/ne/ne.cpp
        algorithms/model4/model4.cpp
        algorithms/model5/model5.cpp
        algorithms/model6/model6.cpp
        algorithms/model7/model7.cpp
        algorithms/dbh/dbh.cpp
        algorithms/hdrf/hdrf.cpp
        algorithms/ldg/ldg.cpp
        algorithms/fennel/fennel.cpp
        algorithms/greedy/greedy.cpp
        algorithms/rand/rand.cpp
        )
add_executable(GraphPartitioning main.cpp)
target_link_libraries(GraphPartitioning GLOG GFLAGS SOURCE_FILES)

