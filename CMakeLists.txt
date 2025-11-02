cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

 set(CMAKE_BUILD_TYPE Debug)
# Only set default build type if none specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
cmake_policy(SET CMP0144 NEW)
cmake_policy(SET CMP0167 NEW)

# ingnore deprecation warnings
add_compile_options(-Wno-deprecated-declarations)

# Only enable AddressSanitizer for Debug builds --> this helps finding memory issues 
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

# Project name
project(Hexaly)

# Add any include file 
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(/opt/hexaly_14_0/include)

# Add all .cpp files in src folder to the executable
file(GLOB_RECURSE SRC_FILES "src/*.cpp")

# Add excutionable
add_executable(Hexaly ${SRC_FILES})

# Threads (required by Hexaly on some platforms)
find_package(Threads REQUIRED)

# Hexaly SDK linking (v14)
set(HEXALY_ROOT "$ENV{HEXALY_HOME}" CACHE PATH "Hexaly SDK root directory")
if(NOT HEXALY_ROOT)
    set(HEXALY_ROOT "/opt/hexaly_14_0")
endif()

find_path(HEXALY_INCLUDE_DIR
    NAMES optimizer/hexalyoptimizer.h
    HINTS ${HEXALY_ROOT}
    PATH_SUFFIXES include
)

find_library(HEXALY_LIBRARY
    NAMES hexaly140 hexaly
    HINTS ${HEXALY_ROOT}
    PATH_SUFFIXES lib lib64 lib/osx-arm64 lib/macos lib/darwin
)

if(HEXALY_INCLUDE_DIR)
    target_include_directories(Hexaly PRIVATE ${HEXALY_INCLUDE_DIR})
endif()

if(HEXALY_LIBRARY)
    message(STATUS "Linking Hexaly library: ${HEXALY_LIBRARY}")
    target_link_libraries(Hexaly PRIVATE ${HEXALY_LIBRARY} Threads::Threads)
    if(APPLE)
        get_filename_component(_HEXALY_LIB_DIR ${HEXALY_LIBRARY} DIRECTORY)
        set_target_properties(Hexaly PROPERTIES BUILD_RPATH ${_HEXALY_LIB_DIR} INSTALL_RPATH ${_HEXALY_LIB_DIR})
    endif()
else()
    # Fallback: try linking by name and add default lib directory
    message(WARNING "Hexaly library not found automatically. Falling back to -lhexaly140; set HEXALY_HOME if this fails.")
    link_directories(${HEXALY_ROOT}/lib)
    target_link_libraries(Hexaly PRIVATE hexaly140 Threads::Threads)
    if(APPLE)
        set_target_properties(Hexaly PROPERTIES BUILD_RPATH ${HEXALY_ROOT}/lib INSTALL_RPATH ${HEXALY_ROOT}/lib)
    endif()
endif()

# Set compiler flags 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# For inellisense config 
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



# Enable verbose linking
#set(CMAKE_VERBOSE_MAKEFILE ON)