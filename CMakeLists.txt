cmake_minimum_required (VERSION 3.17)
project (three-point-calibration)

set (project_SRC
  src/transformation.cpp
  src/test.cpp
)
set (CMAKE_CXX_STANDARD
  20
)
set (CMAKE_CXX_STANDARD_REQUIRED
  ON
)
add_executable (
  ${PROJECT_NAME}
  ${project_SRC}
)
