set(MQT_MISIM_TARGET_NAME mqt-misim)

# add_library( ${MQT_MISIM_TARGET_NAME} )

# set include directories target_include_directories( ${MQT_MISIM_TARGET_NAME} PUBLIC
# $<BUILD_INTERFACE:${MQT_MISIM_INCLUDE_BUILD_DIR}>)

# link to the MQT::Core and Taskflow libraries

# target_link_libraries(${MQT_MISIM_TARGET_NAME} PRIVATE MQT::ProjectWarnings MQT::ProjectOptions)

# add MQT alias add_library(MQT::Misim ALIAS ${MQT_MISIM_TARGET_NAME})

add_subdirectory(python)
