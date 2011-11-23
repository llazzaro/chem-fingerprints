project(chemfp)

cmake_minimum_required(VERSION 2.8)

option(BUILD_SHARED "enable static build support" ON)


ADD_LIBRARY(chemfp SHARED bitops.c chemfp.c heapq.c searches.c fps.c)

add_executable(test_libchemfp test_libchemfp.c)
target_link_libraries(test_libchemfp chemfp)
add_test(test_libchemfp test_libchemfp)
#endforeach(test ${tests})
