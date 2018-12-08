add_library(sensor_range rangereading.cpp rangesensor.cpp)
target_link_libraries(sensor_range sensor_base)

#install(TARGETS sensor_range DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})
install(TARGETS sensor_range DESTINATION "${CMAKE_SOURCE_DIR}/../install/lib")