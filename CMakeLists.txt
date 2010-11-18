
PROJECT(hysteresis_model)

cmake_minimum_required(VERSION 2.6)

#add definitions, compiler switches, etc.
ADD_DEFINITIONS(-Wall)

#list all source files here
#ADD_EXECUTABLE(driver main.cpp interface.cpp)
find_package(PkgConfig)

pkg_check_modules(BASE_TYPES REQUIRED "base-types")
INCLUDE_DIRECTORIES(${BASE_TYPES_INCLUDE_DIRS})
LINK_DIRECTORIES(${BASE_TYPES_LIBRARY_DIRS})

INCLUDE_DIRECTORIES(/usr/include/eigen2)

ADD_LIBRARY(hysteresis_model SHARED 
	    RK4Integrator.cpp
	    RK4Integrator.hpp
	    BoucWenBaberNooriModel.cpp
	    BoucWenBaberNooriModel.hpp)
 
CONFIGURE_FILE(Doxyfile.in Doxyfile @ONLY)
ADD_CUSTOM_TARGET(doc doxygen Doxyfile)
INSTALL(TARGETS hysteresis_model 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib)
INSTALL(FILES  
        RK4Integrator.hpp
	BoucWenBaberNooriModel.hpp
	DESTINATION include)

CONFIGURE_FILE(hysteresis_model.pc.in hysteresis_model.pc @ONLY)
INSTALL(FILES ${CMAKE_BINARY_DIR}/hysteresis_model.pc DESTINATION lib/pkgconfig)
