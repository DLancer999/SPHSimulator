cmake_minimum_required(VERSION 3.11.1)
project (SPHSimulator)

#set(CMAKE_CXX_COMPILER g++-7)
#set(CMAKE_CXX_COMPILER clang++)

set(OPTFLAGS "-O3")
set(DEBUGFLAGS "-g")
set(LINKFLAGS "-flto")
set(ERRORFLAGS "-pedantic-errors")
set(OTHERFLAGS "-std=c++17 -fopenmp")
set(BOOSTHEADERS "-isystem/home/dlancer/Programs/boost/boostGcc/include")
set(BOOSTLIBSDIR "-L/home/dlancer/Programs/boost/boostGcc/lib")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(OTHERFLAGS "${OTHERFLAGS} -stdlib=libc++")
  set(BOOSTHEADERS "-isystem/home/dlancer/Programs/boost/boostClang/include")
  set(BOOSTLIBSDIR "-L/home/dlancer/Programs/boost/boostClang/lib")
endif()

set (WARNING_FLAGS
 "-Wall -Weffc++ -pedantic \
	-Wextra -Wcast-align -Wcast-qual \
	-Wchar-subscripts  -Wcomment -Wconversion \
	-Wdisabled-optimization -Wfloat-equal  -Wformat  -Wformat=2 \
	-Wformat-nonliteral -Wformat-security -Wformat-y2k \
	-Wimport  -Winit-self  -Winvalid-pch   \
	-Wlong-long -Wmissing-braces \
	-Wmissing-field-initializers -Wmissing-format-attribute   \
	-Wmissing-include-dirs -Wmissing-noreturn \
	-Wpacked  -Wparentheses  -Wpointer-arith \
	-Wredundant-decls -Wreturn-type \
	-Wsequence-point  -Wshadow -Wsign-compare  -Wstack-protector \
	-Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch  -Wswitch-default \
	-Wswitch-enum -Wtrigraphs  -Wuninitialized \
	-Wunknown-pragmas  -Wunreachable-code -Wunused \
	-Wunused-function  -Wunused-label  -Wunused-parameter \
	-Wunused-value  -Wunused-variable  -Wvariadic-macros \
	-Wvolatile-register-var  -Wwrite-strings"
)

set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${BOOSTHEADERS} ${OPTFLAGS} ${OTHERFLAGS} ${DEBUGFLAGS} ${WARNING_FLAGS} ${ERRORFLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPTFLAGS} ${LINKFLAGS} ${OTHERFLAGS} ${DEBUGFLAGS} ${WARNING_FLAGS} ${ERRORFLAGS} ${BOOSTLIBSDIR}")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON )

include_directories(./src)
#add_subdirectory (Core)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
link_libraries(${PNG_LIBRARY})

set(SimulatorObjects
  "src/Settings.cpp"
  "src/SPHSolver.cpp"
  "src/Kernels.cpp"
  "src/WriteFunctions.cpp"
)
add_library (Simulator ${SimulatorObjects})

set(TimerUtilObjects
  "src/Statistics.cpp"
)
add_library (TimerUtil ${TimerUtilObjects})

set(GlRendererObjects
  "src/DisplayView.cpp"
  "src/Shader.cpp"
)

find_package(OpenGL REQUIRED)
link_libraries(${OPENGL_LIBRARIES})
find_package(GLEW REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})
link_libraries(${GLEW_LIBRARIES})
find_package(glfw3 REQUIRED)
link_libraries(glfw)

add_library (GlRenderer ${GlRendererObjects})

ADD_CUSTOM_TARGET(
  copy_shader_files
  ${CMAKE_COMMAND}
  -D SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -D DESTINATION_DIR=${CMAKE_CURRENT_BINARY_DIR} 
  -P ${CMAKE_CURRENT_SOURCE_DIR}/copyshaders.cmake
  COMMENT "Copying Files for target: GlRenderer" VERBATIM 
)

ADD_DEPENDENCIES ( GlRenderer copy_shader_files )

set(SimulatorComponents "Simulator" "GlRenderer" "TimerUtil")

add_executable(SPHSimulator "src/SPHSimulator.cpp")
target_link_libraries (SPHSimulator "${SimulatorComponents}")
#file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/src/shaders DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

#no tests yet
#enable_testing()
#add_subdirectory (TestDir)
