include_directories(${CHROMATIN_SIS_COARSE_SOURCE_DIR}/include)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(chromatin.sis.coarse chromatin.sis.coarse.cpp)
target_link_libraries(chromatin.sis.coarse XYFunctions)
##target_link_libraries(anm5 libblas.a liblapack.a)
#target_link_libraries(chromatin.sis.coarse blas lapack)
target_link_libraries(chromatin.sis.coarse
  # ${Boost_FILESYSTEM_LIBRARY}
  # ${Boost_SYSTEM_LIBRARY} 
	${Boost_LIBRARIES}
	# ${Boost_LIBRARY_DIRS}
)

set(CHROMATIN_ENSEMBLE "ensemble")
add_executable(${CHROMATIN_ENSEMBLE} ensemble.cpp)
target_link_libraries(${CHROMATIN_ENSEMBLE} XYFunctions) 
#target_link_libraries(${CHROMATIN_ENSEMBLE} blas lapack)
target_link_libraries(${CHROMATIN_ENSEMBLE} 
	${Boost_LIBRARIES}
)

#----------------------------------------------
# test octree
set(TEST_OCTREE "testoctree")
add_executable(${TEST_OCTREE} test.octree.cpp)
target_link_libraries(
    ${TEST_OCTREE} 
    XYFunctions
	${Boost_LIBRARIES})
#----------------------------------------------

#ADD_CUSTOM_COMMAND(TARGET chromatin.sis.coarse
#    POST_BUILD
#    COMMAND chromatin.sis.coarse
#    WORKING_DIRECTORY ./)
