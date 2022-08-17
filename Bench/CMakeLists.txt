set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "-Wall -Wfatal-errors")
set(CMAKE_CXX_FLAGS "-Wall -Wfatal-errors")
set(BENCHMARK_ENABLE_TESTING OFF)

FetchContent_Declare(
  googlebenchmark
  URL https://github.com/google/benchmark/archive/refs/tags/v1.6.1.tar.gz
)
FetchContent_MakeAvailable(googlebenchmark)

add_executable(BM_chole bench.cpp)

target_link_libraries(BM_chole blasfeo m benchmark::benchmark_main)

add_executable(bdldl bdldl.c)
target_link_libraries(bdldl blasfeo m)
target_compile_options(bdldl PRIVATE -g)

add_executable(blasfeo_test test.c)
target_link_libraries(blasfeo_test blasfeo m)
