project(chemfp)

cmake_minimum_required(VERSION 2.8)

option(BUILD_SHARED "enable static build support" ON)

IF(CMAKE_COMPILER_IS_GNUCXX)
   SET(CMAKE_C_FLAGS "-O3")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

ADD_LIBRARY(chemfp SHARED bitops.c chemfp.c heapq.c searches.c fps.c
                   popcount_SSSE3.c popcount_gillies.c
                   popcount_lauradoux.c popcount_lut.c
                   popcount_popcnt.c hits.c select_popcount.c)
                   

add_executable(test_libchemfp test_libchemfp.c)
target_link_libraries(test_libchemfp chemfp)
add_test(test_libchemfp test_libchemfp)
#endforeach(test ${tests})
