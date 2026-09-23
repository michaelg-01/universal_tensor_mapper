cmake_minimum_required(VERSION 3.26.3)
cmake_policy(VERSION 3.26.3)

# Set toolchain file before project()
set(CMAKE_TOOLCHAIN_FILE "C:/Users/Admin/vcpkg/scripts/buildsystems/vcpkg.cmake")

project(tensor_mapper)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_STANDARD 20)
set(Boost_NO_WARN_NEW_VERSIONS 1)

if(MSVC)
    add_definitions(-D_ENABLE_EXTENDED_ALIGNED_STORAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
endif()

find_package(VTK REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)

# Create executable
add_executable(tensor_mapper 
    main.cpp 
    AbaqusInpReader.cpp 
    MappingMethods.cpp
    cnpy.cpp
)

# Compiler options
target_compile_options(tensor_mapper PRIVATE 
    $<$<CXX_COMPILER_ID:MSVC>:/EHsc /MP /bigobj>
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(tensor_mapper PUBLIC OpenMP::OpenMP_CXX)
endif()

# Link libraries
target_link_libraries(tensor_mapper PRIVATE 
    ${VTK_LIBRARIES}
    ZLIB::ZLIB
)

# Set linker flags for MSVC
if(MSVC)
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
    target_link_options(tensor_mapper PUBLIC
        $<$<CONFIG:Debug>:/INCREMENTAL:NO>
    )
endif()

# Installation
install(TARGETS tensor_mapper DESTINATION bin)