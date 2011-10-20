cmake_minimum_required (VERSION 2.6) 

SET( CMAKE_VERBOSE_MAKEFILE on )
PROJECT( Coupling_Calibration )

link_directories (${main_BINARY_DIR}/src)

add_library( sources
    BoucWenHysteresisModel.cpp
    LevMarOptimization.cpp
    ModelOutput.cpp
    RK4Integrator.cpp)

add_executable (CalibrateCoupling CalibrateCoupling.cpp)
target_link_libraries (CalibrateCoupling sources) 

FIND_LIBRARY(LEVMAR_LIB levmar)
FIND_LIBRARY(LAPACK_LIB lapack)
FIND_LIBRARY(BLAS_LIB blas)
FIND_LIBRARY(F2C_LIB f2c)
find_package(Boost COMPONENTS system filesystem REQUIRED)

target_link_libraries(CalibrateCoupling ${LEVMAR_LIB} ${LAPACK_LIB} ${BLAS_LIB} ${F2C_LIB} ${Boost_FILESYSTEM_LIBRARY})

