project (TrefoilKnot)

cmake_minimum_required(VERSION 2.8)

#-----------------------------------------------------------------------------
# We will enforce an out of source build.
string(COMPARE EQUAL "${${project_name}_source_dir}"
  "${${PROJECT_NAME}_BINARY_DIR}" INSOURCE)
if(INSOURCE)
  message(FATAL_ERROR "${PROJECT_NAME} requires an out of source build.
Please create a separate binary directory and run CMake there.")
endif(INSOURCE)

find_package(VTK COMPONENTS
  vtkCommonColor
  vtkCommonCore
  vtkFiltersHybrid
  vtkFiltersSources
  vtkIOImage
  vtkInteractionStyle
  vtkRenderingOpenGL2
)
include(${VTK_USE_FILE})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Single output directory for all the executables."
  )

#-----------------------------------------------------------------------------
set(SRCS
  TrefoilKnot.cxx
  main.cxx
)

set(INCS
  TrefoilKnot.h
)

if (APPLE)
  add_executable(TrefoilKnot MACOSX_BUNDLE ${SRCS} ${INCS})
else()
  add_executable(TrefoilKnot ${SRCS} ${INCS})
endif()

if(VTK_LIBRARIES)
  target_link_libraries(TrefoilKnot ${VTK_LIBRARIES})
else()
  target_link_libraries(TrefoilKnot vtkHybrid )
endif()
