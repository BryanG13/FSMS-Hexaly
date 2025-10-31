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

include_directories(${OpenMP_CXX_INCLUDE_DIRS})
link_libraries(${OpenMP_CXX_LIBRARIES})

# Threads 
find_package(Threads REQUIRED)

# Add any include file 
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add all .cpp files in src folder to the executable
file(GLOB_RECURSE SRC_FILES "src/*.cpp")

# Add excutionable
add_executable(Hexaly ${SRC_FILES})

# Link libraries to executable 
target_link_libraries(Hexaly ${OpenMP_CXX_LIBRARIES})
target_link_libraries(Hexaly Threads::Threads)
target_link_libraries(Hexaly ilocplex cplex concert pthread m)

# Hexaly SDK integration
# Configuration: allow overriding installation root via -DHEXALY_ROOT or HEXALY_HOME environment variable
set(_HEXALY_HINTS "$ENV{HEXALY_HOME}" "/opt/hexaly_14_0")

find_path(HEXALY_INCLUDE_DIR
    NAMES localsolver.h hexaly.h hexaly/hexaly.h
    HINTS ${_HEXALY_HINTS}
    PATH_SUFFIXES include include/hexaly
)

find_library(HEXALY_LIBRARY
    NAMES hexaly localsolver
    HINTS ${_HEXALY_HINTS}
    PATH_SUFFIXES lib lib64 lib/macos lib/osx-arm64 lib/osx-universal lib/darwin
)

if(HEXALY_INCLUDE_DIR AND HEXALY_LIBRARY)
    message(STATUS "Found Hexaly: include at ${HEXALY_INCLUDE_DIR}, lib at ${HEXALY_LIBRARY}")
    target_include_directories(Hexaly PRIVATE ${HEXALY_INCLUDE_DIR})
    target_link_libraries(Hexaly PRIVATE ${HEXALY_LIBRARY})
    # On macOS, help the runtime loader find the Hexaly dylib if it's outside system paths
    if(APPLE)
        # Add rpath pointing to the detected library directory
        get_filename_component(_HEXALY_LIB_DIR ${HEXALY_LIBRARY} DIRECTORY)
        set_target_properties(Hexaly PROPERTIES BUILD_RPATH ${_HEXALY_LIB_DIR} INSTALL_RPATH ${_HEXALY_LIB_DIR})
    endif()
else()
    message(WARNING "Hexaly SDK not found. Set HEXALY_HOME or pass -DHEXALY_ROOT to help CMake locate it. Expected root like /opt/hexaly_14_0")
endif()

# Set compiler flags 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# For inellisense config 
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



# Enable verbose linking
#set(CMAKE_VERBOSE_MAKEFILE ON)