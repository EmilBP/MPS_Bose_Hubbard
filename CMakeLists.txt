cmake_minimum_required(VERSION 3.0)
project(BoseHubbard_ITensor)

#Set ITensor directory
set(ITENSOR_DIR /home/frederik/ITensor/)

#Bring the headers into the project (full or relative path)
include_directories(
	${ITENSOR_DIR}
        ${ITENSOR_DIR}itensor/
        ${ITENSOR_DIR}itensor/detail
        ${ITENSOR_DIR}itensor/itdata
        ${ITENSOR_DIR}itensor/mps
        ${ITENSOR_DIR}itensor/mps/lattice/
        ${ITENSOR_DIR}itensor/mps/sites/
        ${ITENSOR_DIR}itensor/tensor
        ${ITENSOR_DIR}itensor/util
	include
        )

#Link the Itensor library (NOTE: Compilation now works fine without these two lines)
add_library(itensor STATIC IMPORTED)
set_property(TARGET itensor PROPERTY IMPORTED_LOCATION ${ITENSOR_DIR}lib/libitensor-g.a)

#Set a variable with all the new flags
set(ITENSOR_FLAGS "-I/home/frederik/intel/mkl/include")
set(ITENSOR_LINK_FLAGS "-L/home/frederik/intel/mkl/lib/intel64_lin -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -L${ITENSOR_DIR}lib -litensor -liomp5 -lpthread")

# Build custom library
file( GLOB LIB_SOURCES "src/*.cpp" )
file( GLOB LIB_HEADERS "include/*.h" )
add_library( MyLib ${LIB_SOURCES} ${LIB_HEADERS} )


#Append the new flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${ITENSOR_FLAGS}")

#Build executables
file( GLOB APP_SOURCES main/*.cpp )
foreach( testsourcefile ${APP_SOURCES} )
    get_filename_component(testname ${testsourcefile} NAME_WE)

    add_executable(
	${testname}
	${testsourcefile}
	)

    # Make sure YourLib is linked to each app
    target_link_libraries( "${testname}" MyLib itensor "${ITENSOR_LINK_FLAGS}")
endforeach( testsourcefile ${APP_SOURCES} )

