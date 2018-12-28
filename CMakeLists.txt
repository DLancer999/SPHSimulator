cmake_minimum_required(VERSION 3.11.1)
project (SPHSimulator)

#set(CMAKE_CXX_COMPILER g++-7)
#set(CMAKE_CXX_COMPILER clang++)

set(OPTFLAGS "-O3")
set(DEBUGFLAGS "-g")
set(LINKFLAGS "-flto")
set(ERRORFLAGS "-pedantic-errors")
set(OTHERFLAGS "-std=c++17 -fopenmp")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(OTHERFLAGS "${OTHERFLAGS} -stdlib=libc++")
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

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

find_package(TBB REQUIRED COMPONENTS tbb)

set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${OPTFLAGS} ${OTHERFLAGS} ${DEBUGFLAGS} ${WARNING_FLAGS} ${ERRORFLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPTFLAGS} ${LINKFLAGS} ${OTHERFLAGS} ${DEBUGFLAGS} ${WARNING_FLAGS} ${ERRORFLAGS}")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON )

include_directories(.)
add_subdirectory(./Parallel)
add_subdirectory(./Rendering)
add_subdirectory(./Simulation)
add_subdirectory(./Statistics)

set(SimulatorComponents
  "GlRenderer"
  "Parallel"
  "Simulator"
  "TimerUtil"
)

add_executable(SPHSimulator "SPHSimulator.cpp")
target_link_libraries (SPHSimulator ${SimulatorComponents} ${Boost_PROGRAM_OPTIONS_LIBRARIES} TBB::tbb)
#no tests yet
#enable_testing()
#add_subdirectory (TestDir)
