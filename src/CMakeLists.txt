project("OpenmeshDecimate")

include_directories(
    "C:/Program Files (x86)/OpenMesh/include/"
)

link_directories(
    "C:/Program Files (x86)/OpenMesh/lib"
)

set(LIBS_DEBUG OpenMeshCored.lib OpenMeshToolsd.lib)
set(LIBS_RELEASE OpenMeshCore.lib OpenMeshTools.lib)


add_executable(${PROJECT_NAME} main.cpp fast_quadric_mesh_simplification.h fast_quadric_mesh_simplification.cpp symetric_matrix.h)

target_link_libraries(${PROJECT_NAME} PRIVATE 
    debug ${LIBS_DEBUG} 
    optimized ${LIBS_RELEASE})