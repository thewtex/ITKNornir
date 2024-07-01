cmake_minimum_required(VERSION 3.16.3)
project(Nornir)

set(Nornir_LIBRARIES Nornir)

if(NOT ITK_SOURCE_DIR)
  find_package(ITK REQUIRED)
  list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
  include(ITKModuleExternal)
  if(WASI OR EMSCRIPTEN)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    add_subdirectory(wasm)
  endif()
else()
  set(ITK_DIR ${CMAKE_BINARY_DIR})
  itk_module_impl()
endif()

# These are not expected to pass.
# Mark the tests NornirKWStyleTest and NornirInDoxygenGroup as expected to fail
if(BUILD_TESTING)
  set_tests_properties(NornirKWStyleTest NornirInDoxygenGroup PROPERTIES WILL_FAIL TRUE)
endif()