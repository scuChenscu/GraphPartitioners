set(CMAKE_CXX_STANDARD 20)

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
        utils/graph.cpp
        algorithms/ne/ne.cpp
        algorithms/ours/ours.cpp
        algorithms/model1/model1.cpp
        algorithms/model2/model2.cpp
        algorithms/model3/model3.cpp
        algorithms/model4/model4.cpp
        algorithms/dbh/dbh.cpp
        algorithms/hdrf/hdrf.cpp
        algorithms/ldg/ldg.cpp
        algorithms/fennel/fennel.cpp)
add_executable(GraphPartitioning main.cpp)
target_link_libraries(GraphPartitioning GLOG GFLAGS SOURCE_FILES)

