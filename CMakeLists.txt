cmake_minimum_required(VERSION 3.18)
project(rcspp)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- CXX Standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Core C++ Static Library ---

# Sources
set(RCSPP_SOURCES
    "src/rcspp/src/vrptw_pricing/solvers/labeling/labeling.cpp"
    "src/rcspp/src/vrptw_pricing/solvers/pulsing/pulsing.cpp"
)

# Static library
add_library(rcspp_core STATIC
    ${RCSPP_SOURCES}
)

# Include directory
target_include_directories(rcspp_core
    PUBLIC "src/rcspp/include"
)

# --- Python binding ---

# Dependencies
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(nanobind REQUIRED)

# Link with core library
nanobind_add_module(_rcspp
    "src/rcspp/bindings/python/rcspp_pybind.cpp"
)

target_link_libraries(_rcspp
    PRIVATE rcspp_core
)

# Install the module to python package
install(
    TARGETS _rcspp
    DESTINATION cgvrp
)

# Generate stub
set(DLL_PATH "")

if(WIN32 AND CMAKE_COMPILER_IS_GNUCXX)
    get_filename_component(MINGW_COMPILER_BIN_DIR ${CMAKE_CXX_COMPILER} DIRECTORY)
    message(STATUS "Found MinGW compiler bin directory: ${MINGW_COMPILER_BIN_DIR}")
    set(DLL_PATH "${MINGW_COMPILER_BIN_DIR}")
endif()

add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/_rcspp.pyi"
    COMMAND
    ${Python_EXECUTABLE}
    "${CMAKE_SOURCE_DIR}/generate_stub.py"
    --module-name _rcspp
    --pyd-path "$<TARGET_FILE_DIR:_rcspp>"
    --output-file "${CMAKE_BINARY_DIR}/_rcspp.pyi"
    --dll-path "${DLL_PATH}"
    DEPENDS _rcspp
    VERBATIM
)

add_custom_target(
    rcspp_stub ALL
    DEPENDS "${CMAKE_BINARY_DIR}/_rcspp.pyi"
)

install(
    DIRECTORY "${CMAKE_BINARY_DIR}/"
    DESTINATION cgvrp
    FILES_MATCHING PATTERN "*.pyi"
)