cmake_minimum_required(VERSION 2.8)
project(GridTools)

find_package(inmost REQUIRED)

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

add_executable(FixFaults test_fix_faults.cpp)
add_executable(SetLayers test_set_layers.cpp)
add_executable(MakeFaults test_make_faults.cpp)
add_executable(FixTiny fix_tiny.cpp)
add_executable(FixTinyCollapse fix_tiny_collapse.cpp)
add_executable(UniteFaces unite_faces.cpp)
add_executable(Dual dual.cpp)
add_executable(Tetra tetra.cpp)
add_executable(Slice slice.cpp)
add_executable(SliceFunc test_slice_func.cpp)
add_executable(SplitNonplanar split_nonplanar.cpp)
add_executable(CollapseDegenerate collapse_degenerate.cpp)
add_executable(Bnd2Stl bnd2stl.cpp)
add_executable(Sector sector.cpp)
add_executable(Sew sew.cpp)
add_executable(Scale scale.cpp)
add_executable(Move move.cpp)
add_executable(Convert convert.cpp)
add_executable(CubeTransform test_cube_transform.cpp)
add_executable(MeshInfo mesh_info.cpp)
add_executable(SameMeshDifference difference_same.cpp)
add_executable(MatchSameMeshDifference difference_same_match.cpp)
add_executable(MeshDifference difference_map.cpp)
add_executable(Kmeans kmeans.cpp)
add_executable(Agglomerate agglomerate.cpp)
add_executable(acute_mesh acute_grid.cpp)
add_executable(disturbed_mesh disturbed_grid.cpp)
add_executable(hex_mesh hex_grid.cpp)
add_executable(kershaw_mesh kershaw_grid.cpp)
add_executable(nonconvex_mesh nonconvex_grid.cpp)
add_executable(e_mesh e_grid.cpp)
add_executable(twist_mesh twist_grid.cpp)
add_executable(shestakov_mesh shestakov_grid.cpp)
add_executable(sinusoidal_mesh sinusoidal_grid.cpp)
add_executable(split_faces split_faces.cpp)
add_executable(glue_faces glue_faces.cpp)
add_executable(check_collapse check_collapse.cpp)
add_executable(Fracture test_fracture.cpp)
add_executable(segment_data segment_data.cpp)
add_executable(ChangeTagName change_tag_name.cpp)
add_library(FractureLib fracture.cpp fracture.h)
add_library(SliceFuncLib slice_func.cpp slice_func.h)
add_library(FixFaultsLib fix_faults.cpp fix_faults.h)
add_library(CubeTransformLib cube_transform.cpp cube_transform.h)
add_library(SetLayersLib set_layers.cpp set_layers.h)


target_link_libraries(MakeFaults ${INMOST_LIBRARIES})
target_link_libraries(MakeFaults FractureLib)
target_link_libraries(MakeFaults FixFaultsLib)
target_link_libraries(MakeFaults SliceFuncLib)
install(TARGETS MakeFaults EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(SetLayers ${INMOST_LIBRARIES})
target_link_libraries(SetLayers SetLayersLib)
install(TARGETS SetLayers EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(Fracture ${INMOST_LIBRARIES})
target_link_libraries(Fracture FractureLib)
install(TARGETS Fracture EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(FixFaults ${INMOST_LIBRARIES})
target_link_libraries(FixFaults FixFaultsLib)
install(TARGETS FixFaults EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(FixTinyCollapse ${INMOST_LIBRARIES})
install(TARGETS FixTinyCollapse EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(split_faces ${INMOST_LIBRARIES})
install(TARGETS split_faces EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(check_collapse ${INMOST_LIBRARIES})
install(TARGETS check_collapse EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(glue_faces ${INMOST_LIBRARIES})
install(TARGETS glue_faces EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(FixTiny ${INMOST_LIBRARIES})
install(TARGETS FixTiny EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(UniteFaces ${INMOST_LIBRARIES})
install(TARGETS UniteFaces EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)



target_link_libraries(Dual ${INMOST_LIBRARIES})
install(TARGETS Dual EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(Tetra ${INMOST_LIBRARIES})
install(TARGETS Tetra EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)




target_link_libraries(Slice ${INMOST_LIBRARIES})
install(TARGETS Slice EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)



target_link_libraries(SliceFunc ${INMOST_LIBRARIES})
target_link_libraries(SliceFunc SliceFuncLib)
install(TARGETS SliceFunc EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)



target_link_libraries(SplitNonplanar ${INMOST_LIBRARIES})
install(TARGETS SplitNonplanar EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)



target_link_libraries(CollapseDegenerate ${INMOST_LIBRARIES})
install(TARGETS CollapseDegenerate EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)



target_link_libraries(Bnd2Stl ${INMOST_LIBRARIES})
install(TARGETS Bnd2Stl EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(Sector ${INMOST_LIBRARIES})
install(TARGETS Sector EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(Sew ${INMOST_LIBRARIES})
install(TARGETS Sew EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(Scale ${INMOST_LIBRARIES})
install(TARGETS Scale EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)



target_link_libraries(Move ${INMOST_LIBRARIES})
install(TARGETS Move EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(Convert ${INMOST_LIBRARIES})
install(TARGETS Convert EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(CubeTransform ${INMOST_LIBRARIES})
target_link_libraries(CubeTransform CubeTransformLib)
install(TARGETS CubeTransform EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(MeshInfo ${INMOST_LIBRARIES})
install(TARGETS MeshInfo EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(SameMeshDifference ${INMOST_LIBRARIES})
install(TARGETS SameMeshDifference EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(MatchSameMeshDifference ${INMOST_LIBRARIES})
install(TARGETS MatchSameMeshDifference EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)



target_link_libraries(MeshDifference ${INMOST_LIBRARIES})
install(TARGETS MeshDifference EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(Kmeans ${INMOST_LIBRARIES})
install(TARGETS Kmeans EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(Agglomerate ${INMOST_LIBRARIES})
install(TARGETS Agglomerate EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(acute_mesh ${INMOST_LIBRARIES})
install(TARGETS acute_mesh EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(disturbed_mesh ${INMOST_LIBRARIES})
install(TARGETS disturbed_mesh EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(hex_mesh ${INMOST_LIBRARIES})
install(TARGETS hex_mesh EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(kershaw_mesh ${INMOST_LIBRARIES})
install(TARGETS kershaw_mesh EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(nonconvex_mesh ${INMOST_LIBRARIES})
install(TARGETS nonconvex_mesh EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(e_mesh ${INMOST_LIBRARIES})
install(TARGETS e_mesh EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(twist_mesh ${INMOST_LIBRARIES})
install(TARGETS twist_mesh EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(shestakov_mesh ${INMOST_LIBRARIES})
install(TARGETS shestakov_mesh EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(sinusoidal_mesh ${INMOST_LIBRARIES})
install(TARGETS sinusoidal_mesh EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

target_link_libraries(segment_data ${INMOST_LIBRARIES})
install(TARGETS segment_data EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)


target_link_libraries(ChangeTagName ${INMOST_LIBRARIES})
install(TARGETS ChangeTagName EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)

set_property(TARGET FractureLib PROPERTY PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/fracture.h")
set_property(TARGET SliceFuncLib PROPERTY PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/slice_func.h")
set_property(TARGET FixFaultsLib PROPERTY PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/fix_faults.h")
set_property(TARGET CubeTransformLib PROPERTY PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/cube_transform.h")

install(TARGETS FractureLib SliceFuncLib FixFaultsLib CubeTransformLib EXPORT inmost-targets 
        LIBRARY DESTINATION lib 
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

