cmake_minimum_required(VERSION 3.22)
project(Lorenz_Mie_Solver LANGUAGES CXX C)

# ---- C++ standard ----------------------------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Output Directories ----------------------------------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# ---- Sources ---------------------------------------------------------------------------------------------------------
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")
add_executable(Lorenz_Mie_Solver ${SOURCES})

target_compile_definitions(Lorenz_Mie_Solver PRIVATE
    OUTPUT_DIR="${CMAKE_SOURCE_DIR}/output"
)

# Ensure C++ 20 features
target_compile_features(Lorenz_Mie_Solver PRIVATE cxx_std_20)

# ---- Include Paths ---------------------------------------------------------------------------------------------------
target_include_directories(Lorenz_Mie_Solver PRIVATE ${CMAKE_SOURCE_DIR}/src)

# ---- Third-party Libraries -------------------------------------------------------------------------------------------

find_package(OpenMP REQUIRED)

# Link compiled libraries to executable
target_link_libraries(
    Lorenz_Mie_Solver
    PRIVATE
    OpenMP::OpenMP_CXX
    gomp
)

target_link_options(Lorenz_Mie_Solver PRIVATE -static-libgcc -static-libstdc++ -static)
