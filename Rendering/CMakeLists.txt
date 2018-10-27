
set(GlRendererObjects
  "DisplayView.cpp"
  "Shader.cpp"
  "WriteFunctions.cpp"
  "HardCodedShaders.cpp"
)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
link_libraries(${PNG_LIBRARY})
find_package(OpenGL REQUIRED)
link_libraries(${OPENGL_LIBRARIES})
find_package(GLEW REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})
link_libraries(${GLEW_LIBRARIES})
find_package(glfw3 REQUIRED)
link_libraries(glfw)

add_library (GlRenderer ${GlRendererObjects})
