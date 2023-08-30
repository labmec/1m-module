set(MESHESDIR ${CMAKE_SOURCE_DIR}/targets/meshes)
set (CONFIG_H "${CMAKE_SOURCE_DIR}/targets/meshpath_config.h")
configure_file ("${CMAKE_SOURCE_DIR}/targets/meshpath_config.h.in" ${CONFIG_H})

add_executable(1m-module 1m-module.cpp)
target_link_libraries(1m-module PRIVATE libproject)