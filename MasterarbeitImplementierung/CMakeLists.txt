cmake_minimum_required (VERSION 3.5.1)

### To use gcc/g++ on a Macintosh, you must set the Compilers
### here, not inside the project
##if(APPLE)
##       set(CMAKE_C_COMPILER "/usr/local/bin/gcc-7")
##       set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-7")
##endif()
### TODO: for now, we use CLang for Mac
###
### In order to create OpenFHE's static libraries you should enable
### the BUILD_STATIC option. For that, you run "cmake .. -DBUILD_STATIC=ON".
### After having your link completed you will find static libs
### with the suffix "_static" in ./build/libs/.
### Examples: OPENFHEpke_static.a, OPENFHEcore_static.a, etc.
### After you run "make install" in your build directory, you can build your custom application.
### If you need your application to be linked statically, then run "cmake .. -DBUILD_STATIC=ON"

project(Masterarbeit_Projekt)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug) 
option( BUILD_STATIC "Set to ON to include static versions of the library" OFF)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#create necessarcy subfolders and copy data to build folder
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/benchmark)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/demoData)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/data)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/data/keys)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/data/CiphertextOutput)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/data/CipherEvaluation)
file(COPY ${CMAKE_SOURCE_DIR}/FleetData.json DESTINATION ${CMAKE_BINARY_DIR})
#file(COPY ${CMAKE_SOURCE_DIR}/FleetDataUC2.json DESTINATION ${CMAKE_BINARY_DIR})

#find required packes
find_package(OpenFHE CONFIG REQUIRED) 
if (OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE shared lib files location: ${OpenFHE_SHARED_LIBRARIES}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif ()

find_package(jsoncpp CONFIG REQUIRED)
if(jsoncpp_FOUND)
    get_target_property(JSON_INC_PATH jsoncpp_lib INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "FOUND PACKAGE jsoncpp")
    message(STATUS "JSONCPP Include Dirs: ${JSON_INC_PATH}")
else()
    message(FATAL_ERROR "PACKAGE jsoncpp NOT FOUND")
endif ()

# find_package(benchmark REQUIRED)
# if(benchmark_FOUND)
#     #get_target_property(JSON_INC_PATH jsoncpp_lib INTERFACE_INCLUDE_DIRECTORIES)
#     message(STATUS "FOUND PACKAGE Google Benchmark")
#     #message(STATUS "JSONCPP Include Dirs: ${JSON_INC_PATH}")
# else()
#     message(FATAL_ERROR "PACKAGE Google Benchmark NOT FOUND")
# endif ()

set( CMAKE_CXX_FLAGS ${OpenFHE_CXX_FLAGS} )


include_directories( ${OPENMP_INCLUDES} ) 
include_directories( ${OpenFHE_INCLUDE} )
include_directories( ${OpenFHE_INCLUDE}/core ) 
include_directories( ${OpenFHE_INCLUDE}/pke ) 
include_directories( ${OpenFHE_INCLUDE}/binfhe )
include_directories( ${JSON_INC_PATH} )
include_directories( ${CMAKE_SOURCE_DIR}/src/include)
include_directories( ${benchmark_INCLUDE})


### add directories for other OpenFHE modules as needed for your project

link_directories( ${OpenFHE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
link_directories( ${OPENSSL_LIBRARIES})
link_directories( ${OPENSSL_CRYPTO_LIBRARY})

if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${OpenFHE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
    link_libraries( ${OpenFHE_SHARED_LIBRARIES} )
endif()

### ADD YOUR EXECUTABLE(s) HEREc
### add_executable( EXECUTABLE-NAME SOURCES )
###
### EXAMPLE:
# add_executable(FleetManager src/FleetManager.cpp)
# target_link_libraries(FleetManager PRIVATE jsoncpp_lib)

# add_executable(FleetAgent src/FleetAgent.cpp)
# target_link_libraries(FleetAgent PRIVATE jsoncpp_lib)

# add_executable(Backend src/Backend.cpp)
# target_link_libraries(Backend PRIVATE jsoncpp_lib)

add_executable(BenchmarkUC1 src/BenchmarkUC1.cpp)
target_link_libraries(BenchmarkUC1 jsoncpp_lib)


add_executable(simple-real-numbers-serial src/simple-real-numbers-serial.cpp)
target_link_libraries(simple-real-numbers-serial jsoncpp_lib)

add_executable(PearsonCorrelationHybrid src/PearsonCorrelationHybrid.cpp)
target_link_libraries(PearsonCorrelationHybrid jsoncpp_lib)

add_executable(PearsonCorrelation src/PearsonCorrelation.cpp)
target_link_libraries(PearsonCorrelation jsoncpp_lib)

add_executable(PearsonCorrelationBS src/PearsonCorrelationBS.cpp)
target_link_libraries(PearsonCorrelationBS jsoncpp_lib)


add_executable(CO2e src/CO2e.cpp)
target_link_libraries(CO2e jsoncpp_lib)

add_executable(StandardDeviation src/StandardDeviation.cpp)
target_link_libraries(StandardDeviation jsoncpp_lib)

add_executable(Variance src/Variance.cpp)
target_link_libraries(Variance jsoncpp_lib)

# add_executable(co2 src/co2.cpp)
# target_link_libraries(co2 jsoncpp_lib)

# add_executable(BenchmarkUC1hybrid src/BenchmarkUC1hybrid.cpp)
# target_link_libraries(BenchmarkUC1hybrid jsoncpp_lib)

# add_executable(BenchmarkUC2 src/BenchmarkUC2.cpp)
# target_link_libraries(BenchmarkUC2 jsoncpp_lib)

# add_executable(BenchmarkUC2hybrid src/BenchmarkUC2hybrid.cpp)
# target_link_libraries(BenchmarkUC2hybrid jsoncpp_lib)