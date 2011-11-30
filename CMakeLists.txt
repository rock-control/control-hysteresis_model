
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

pkg_check_modules(MOTOR_CONTROLLER REQUIRED "motor_controller")
INCLUDE_DIRECTORIES(${MOTOR_CONTROLLER_INCLUDE_DIRS})
LINK_DIRECTORIES(${MOTOR_CONTROLLER_LIBRARY_DIRS})

ADD_LIBRARY(hysteresis_model SHARED 
            BoucWenHysteresisModel.cpp
            BoucWenHysteresisModel.hpp)
target_link_libraries(hysteresis_model ${MOTOR_CONTROLLER_LIBRARIES})
 
CONFIGURE_FILE(Doxyfile.in Doxyfile @ONLY)
ADD_CUSTOM_TARGET(doc doxygen Doxyfile)
INSTALL(TARGETS hysteresis_model 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib)
INSTALL(FILES  
        BoucWenHysteresisModel.hpp
	DESTINATION include)

CONFIGURE_FILE(hysteresis_model.pc.in hysteresis_model.pc @ONLY)
INSTALL(FILES ${CMAKE_BINARY_DIR}/hysteresis_model.pc DESTINATION lib/pkgconfig)
