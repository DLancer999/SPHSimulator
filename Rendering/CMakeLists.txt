
set(GlRendererObjects
  "DisplayView.cpp"
  "Shader.cpp"
  "WriteFunctions.cpp"
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

ADD_DEPENDENCIES ( GlRenderer copy_shader_files )

ADD_CUSTOM_TARGET(
  copy_shader_files
  ${CMAKE_COMMAND}
  -D SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -D DESTINATION_DIR=${CMAKE_CURRENT_BINARY_DIR} 
  -P ${CMAKE_CURRENT_SOURCE_DIR}/copyshaders.cmake
  COMMENT "Copying Files for target: GlRenderer" VERBATIM 
)
