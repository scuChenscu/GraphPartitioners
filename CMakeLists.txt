set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_OSX_ARCHITECTURES "arm64")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

cmake_minimum_required(VERSION 3.17)
project(GraphPartitioning)

set(GLOG /opt/homebrew/Cellar/glog/0.6.0)
set(GFLAGS /opt/homebrew/Cellar/gflags/2.2.2)
set(BOOST /opt/homebrew/Cellar/boost/1.84.0_1)
set(PROFILER /opt/homebrew/Cellar/gperftools/2.15)

include_directories(${GLOG}/include/ ${GFLAGS}/include/ ${BOOST}/include/ ${PROFILER}/include/)
link_directories(${GLOG}/lib/ ${GFLAGS}/lib/ ${PROFILER}/lib/)
add_library(SOURCE_FILES
        converter/conversions.cpp
        converter/shuffler.cpp
        baseGraph/base_graph.cpp
        partitioner/partitioner.cpp
        partitioner/edgePartitioner.cpp
        partitioner/vertexPartitioner.cpp
        utils/graph.cpp
        algorithms/ne/ne.cpp
        algorithms/dne/dne.cpp
        algorithms/model4/model4.cpp
        algorithms/model5/model5.cpp
        algorithms/model6/model6.cpp
        algorithms/model7/model7.cpp
        algorithms/model8/model8.cpp
        algorithms/model9/model9.cpp
        algorithms/model10/model10.cpp
        algorithms/model11/model11.cpp
        algorithms/model12/model12.cpp
        algorithms/dbh/dbh.cpp
        algorithms/hdrf/hdrf.cpp
        algorithms/ldg/ldg.cpp
        algorithms/fennel/fennel.cpp
        algorithms/greedy/greedy.cpp
        algorithms/rand/rand.cpp
        algorithms/adwise/adwise.cpp
        algorithms/offstreamNG/offstreamNG.cpp
        algorithms/offstreamNH/offstreamNH.cpp
        algorithms/offstreamNW/offstreamNW.cpp
        algorithms/offstreamNWG/offstreamNWG.cpp
        algorithms/timene/timene.cpp
        algorithms/dcne/dcne.cpp
        algorithms/bne/bne.cpp
        algorithms/dpqne/dpqne.cpp
        algorithms/offstreamNA/offstreamNA.cpp
        )
add_executable(GraphPartitioning main.cpp)
target_link_libraries(GraphPartitioning SOURCE_FILES PROFILER GLOG GFLAGS tcmalloc)



