cmake_minimum_required(VERSION 3.11)
project( LungNoduleSegmentater )
find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )
find_package( VTK REQUIRED )
include( ${VTK_USE_FILE} )
if(VTK_QT_VERSION VERSION_GREATER "4")
  message("VTK was built with Qt support.")
  find_package(Qt5Widgets REQUIRED QUIET)
endif()
include_directories(../common)
add_executable( LungNoduleSegmentation
  itkLesionSegmentationCommandLineProgressReporter.cxx
  itkLesionSegmentationCommandLineProgressReporter.h
	LungNoduleSegmentation.cpp
	LesionSegmentationCLI.h
#	itkLungWallFeatureGenerator2.hxx
#	itkLungWallFeatureGenerator2.h
	../common/vtkCutPlaneWidget.h
	../common/vtkCutPlaneWidget.cxx
	../common/itkVTKViewImageAndSegmentation.cxx
	../common/itkVTKViewImageAndSegmentation.h)
target_link_libraries( LungNoduleSegmentation ${ITK_LIBRARIES} ${VTK_LIBRARIES})

include(CTest)
if (BUILD_TESTING)

  # Test 1
  file(MAKE_DIRECTORY ${LungNoduleSegmentater_BINARY_DIR}/NoduleSegmentationTest1)
  file(MAKE_DIRECTORY ${LungNoduleSegmentater_BINARY_DIR}/NoduleSegmentationTest1/nodule_image)
  add_test( NAME NoduleSegmentationTest1 COMMAND LungNoduleSegmentation
    --InputDICOMDir ${LungNoduleSegmentater_SOURCE_DIR}/../../data/E00140/ --Seeds 3 -62 66 -91.5
    --MaximumRadius 12.0 --Visualize 1 --Outline 1 --Screenshot nodule_image --Supersample 1
    --SupersampledIsotropicSpacing .2 --WriteFeatureImages 1 --OutputMesh E00140.N1.vtp
    WORKING_DIRECTORY ${LungNoduleSegmentater_BINARY_DIR}/NoduleSegmentationTest1 )

  # Test 2
  file(MAKE_DIRECTORY ${LungNoduleSegmentater_BINARY_DIR}/NoduleSegmentationTest2)
  file(MAKE_DIRECTORY ${LungNoduleSegmentater_BINARY_DIR}/NoduleSegmentationTest2/nodule_image)
  add_test( NAME NoduleSegmentationTest2 COMMAND LungNoduleSegmentation
    --InputDICOMDir ${LungNoduleSegmentater_SOURCE_DIR}/../../data/E00140/ --Seeds 3 -46.8 -56.3 -91.4
    --MaximumRadius 12.0 --Visualize 1 --Outline 1 --Screenshot nodule_image --Supersample 1
    --SupersampledIsotropicSpacing .2 --WriteFeatureImages 1 --OutputMesh E00140.N2.vtp
    WORKING_DIRECTORY ${LungNoduleSegmentater_BINARY_DIR}/NoduleSegmentationTest2 )

endif()

