cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(SiliconTracker)

# debug
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3")

set(CMAKE_CXX_STANDARD 11)

find_package(Geant4 REQUIRED ui_all vis_all)
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include
                    ${ROOT_INCLUDE_DIRS})

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

add_executable(tracker tracker.cc ${sources} ${headers})
target_link_libraries(tracker
                      ${Geant4_LIBRARIES}
                      ${ROOT_LIBRARIES})

set(SCRIPTS
    run.mac
    draw.C)

foreach (_script ${SCRIPTS})
    configure_file(${PROJECT_SOURCE_DIR}/${_script}
                   ${PROJECT_BINARY_DIR}/${_script}
                   COPYONLY)
endforeach ()
