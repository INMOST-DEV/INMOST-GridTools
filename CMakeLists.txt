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
add_executable(SliceWell slice_well.cpp)
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
add_executable(Anisotropy anisotropy.cpp)
add_executable(CylGrid cylinders_grid.cpp)
add_executable(XY2XZ xy2xz.cpp)



add_library(FractureLib fracture.cpp fracture.h)
add_library(SliceFuncLib slice_func.cpp slice_func.h)
add_library(FixFaultsLib fix_faults.cpp fix_faults.h)
add_library(CubeTransformLib cube_transform.cpp cube_transform.h)
add_library(SetLayersLib set_layers.cpp set_layers.h)



set_property(TARGET FractureLib PROPERTY PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/fracture.h")
set_property(TARGET SliceFuncLib PROPERTY PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/slice_func.h")
set_property(TARGET FixFaultsLib PROPERTY PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/fix_faults.h")
set_property(TARGET CubeTransformLib PROPERTY PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/cube_transform.h")

list(APPEND LIBS
	FractureLib
	SliceFuncLib
	FixFaultsLib
	CubeTransformLib)

list(APPEND EXECS
	FixFaults
	SetLayers
	MakeFaults
	FixTiny
	FixTinyCollapse
	UniteFaces
	Dual
	Tetra
	Slice
	SliceFunc
	SliceWell
	SplitNonplanar
	CollapseDegenerate
	Bnd2Stl
	Sector
	Sew
	Scale
	Move
	Convert
	CubeTransform
	MeshInfo
	SameMeshDifference
	MatchSameMeshDifference
	MeshDifference
	Kmeans
	Agglomerate
	acute_mesh
	disturbed_mesh
	hex_mesh
	kershaw_mesh
	nonconvex_mesh
	e_mesh
	twist_mesh
	shestakov_mesh
	sinusoidal_mesh
	split_faces
	glue_faces
	check_collapse
	Fracture
	segment_data
	ChangeTagName
	Anisotropy
	CylGrid
	XY2XZ)
list(APPEND LINK_Fracture
	MakeFaults
	Fracture)
list(APPEND LINK_FixFaults
	MakeFaults
	FixFaults)
list(APPEND LINK_SliceFunc
	MakeFaults
	SliceWell
	SliceFunc)
list(APPEND LINK_SetLayers
	SetLayers)
list(APPEND LINK_CubeTransform
	CubeTransform)
	
foreach(EXEC ${EXECS})
	target_link_libraries(${EXEC} ${INMOST_LIBRARIES})
	install(TARGETS ${EXEC} EXPORT inmost-targets RUNTIME DESTINATION bin/GridTools)
endforeach(EXEC ${EXECS})

foreach(EXEC ${LINK_Fracture})
	target_link_libraries(${EXEC} FractureLib)
endforeach(EXEC ${LINK_Fracture})

foreach(EXEC ${LINK_FixFaults})
	target_link_libraries(${EXEC} FixFaultsLib)
endforeach(EXEC ${LINK_FixFaults})
	
foreach(EXEC ${LINK_SliceFunc})
	target_link_libraries(${EXEC} SliceFuncLib)
endforeach(EXEC ${LINK_SliceFunc})

foreach(EXEC ${LINK_SetLayers})
	target_link_libraries(${EXEC} SetLayersLib)
endforeach(EXEC ${LINK_SetLayers})

foreach(EXEC ${LINK_CubeTransform})
	target_link_libraries(${EXEC} CubeTransformLib)
endforeach(EXEC ${LINK_CubeTransform})
	
install(TARGETS ${LIBS} EXPORT inmost-targets 
        LIBRARY DESTINATION lib 
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

