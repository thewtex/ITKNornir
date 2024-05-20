cmake_minimum_required(VERSION 3.16)
project(itk-nornir-wasm)

set(CMAKE_CXX_STANDARD 17)

find_package(ITK REQUIRED
  COMPONENTS
    WebAssemblyInterface
    Nornir
    )
include(${ITK_USE_FILE})

add_executable(ir-refine-grid ir-refine-grid.cxx)
target_link_libraries(ir-refine-grid PUBLIC ${ITK_LIBRARIES})

enable_testing()
add_test(NAME ir-refine-grid-help COMMAND ir-refine-grid --help)