add_library(gridfastslam
    gridslamprocessor.cpp
    gridslamprocessor_tree.cpp
    gridslamprocessor_debug_print.cpp
    motionmodel.cpp
    particle.cpp
    pose_node.cpp
)
target_link_libraries(gridfastslam scanmatcher sensor_range)

#install(TARGETS gridfastslam DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})
install(TARGETS gridfastslam DESTINATION "${CMAKE_SOURCE_DIR}/../install/lib")
