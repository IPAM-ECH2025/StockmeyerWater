cmake_minimum_required(VERSION 3.25)
project(StockmeyerIPAM LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(tpls/kokkos)
add_subdirectory(tpls/kokkos-fft)

add_executable(main main.cpp)

target_link_libraries(main PUBLIC Kokkos::kokkos KokkosFFT::fft)
