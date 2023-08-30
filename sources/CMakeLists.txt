# target_include_directories(libreorder PUBLIC
#   ${CMAKE_CURRENT_SOURCE_DIR})

set(MY_HEADERS
  TPZMixedModelProblem.h
)

set(MY_SRCS
  TPZMixedModelProblem.cpp
)

# let us create a config file with the directory containing the meshes
# create folder with meshes
# set(MESHESDIR ${PROJECT_SOURCE_DIR}/mesh)
# set (CONFIG_H "${CMAKE_CURRENT_SOURCE_DIR}/divfree_config.h")
# configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/divfree_config.h.in" ${CONFIG_H})


target_sources(libproject PRIVATE
  ${MY_HEADERS} ${MY_SRCS}
)