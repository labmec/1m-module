set(MESHESDIR ${CMAKE_SOURCE_DIR}/targets/meshes)
set (CONFIG_H "${CMAKE_SOURCE_DIR}/targets/meshpath_config.h")
configure_file ("${CMAKE_SOURCE_DIR}/targets/meshpath_config.h.in" ${CONFIG_H})

add_executable(1m-module 1m-module.cpp)
target_link_libraries(1m-module PRIVATE libproject)

add_executable(1m-module-mixed 1m-module-mixed.cpp)
target_link_libraries(1m-module-mixed PRIVATE libproject)

add_executable(convergence-mixed convergence-mixed.cpp)
target_link_libraries(convergence-mixed PRIVATE libproject)

add_executable(1m-module-dynamics 1m-module-dynamics.cpp)
target_link_libraries(1m-module-dynamics PRIVATE libproject)

add_executable(Cooks-membrane-UP cooks-membrane-UP.cpp)
target_link_libraries(Cooks-membrane-UP PRIVATE libproject)

add_executable(Shear-locking-beam-UP shear-locking-beam-UP.cpp)
target_link_libraries(Shear-locking-beam-UP PRIVATE libproject)

add_executable(Bishop-beam-UP bishop-beam-UP.cpp)
target_link_libraries(Bishop-beam-UP PRIVATE libproject)

add_executable(UniformStretch uniform-stretch.cpp)
target_link_libraries(UniformStretch PRIVATE libproject)

add_executable(SimpleShear simple-shear.cpp)
target_link_libraries(SimpleShear PRIVATE libproject)

add_executable(FullyHybrid fully-hybrid.cpp)
target_link_libraries(FullyHybrid PRIVATE libproject)

add_executable(1m-module-new-mesh-hdiv 1m-module-new-mesh-hdiv.cpp)
target_link_libraries(1m-module-new-mesh-hdiv PRIVATE libproject)

add_executable(taylor-hood-test taylor-hood-test.cpp)
target_link_libraries(taylor-hood-test PRIVATE libproject)

add_executable(taylor-hood-module taylor-hood-module.cpp)
target_link_libraries(taylor-hood-module PRIVATE libproject)

