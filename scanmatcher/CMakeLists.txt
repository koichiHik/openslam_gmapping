add_library(scanmatcher eig3.cpp scanmatcher.cpp scanmatcherprocessor.cpp smmap.cpp)
target_link_libraries(scanmatcher sensor_range utils)

#install(TARGETS scanmatcher DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})
install(TARGETS scanmatcher DESTINATION "${CMAKE_SOURCE_DIR}/../install/lib")