set(MESHESDIR ${CMAKE_SOURCE_DIR}/targets/meshes)
set (CONFIG_H "${CMAKE_SOURCE_DIR}/targets/meshpath_config.h")
configure_file ("${CMAKE_SOURCE_DIR}/targets/meshpath_config.h.in" ${CONFIG_H})

add_executable(1m-module 1m-module.cpp ProblemData.h ProblemData.cpp json.hpp)
target_link_libraries(1m-module PRIVATE libproject)

add_executable(1m-module-mixed 1m-module-mixed.cpp ProblemData.h ProblemData.cpp json.hpp)
target_link_libraries(1m-module-mixed PRIVATE libproject)

add_executable(TestSaddlePermute test-saddle-permute.cpp ProblemData.h ProblemData.cpp json.hpp)
target_link_libraries(TestSaddlePermute PRIVATE libproject)