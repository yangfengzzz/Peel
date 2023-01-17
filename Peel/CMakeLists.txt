#  Copyright (c) 2022 Feng Yang
#
#  I am making my contributions/submissions to this project solely in my
#  personal capacity and am not conveying any rights to any intellectual
#  property of any third parties.

cmake_minimum_required(VERSION 3.12)

project(peel LANGUAGES C CXX)

file(GLOB FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

set(PHYSX_SOURCE
        ../third_party/physx/physx/source/common/src/mac
        ../third_party/physx/physx/include
        ../third_party/physx/physx/source/common/include
        ../third_party/physx/physx/source/common/src
        ../third_party/physx/physx/source/geomutils/include
        ../third_party/physx/physx/source/geomutils/src
        ../third_party/physx/physx/source/geomutils/src/contact
        ../third_party/physx/physx/source/geomutils/src/common
        ../third_party/physx/physx/source/geomutils/src/convex
        ../third_party/physx/physx/source/geomutils/src/distance
        ../third_party/physx/physx/source/geomutils/src/sweep
        ../third_party/physx/physx/source/geomutils/src/gjk
        ../third_party/physx/physx/source/geomutils/src/intersection
        ../third_party/physx/physx/source/geomutils/src/mesh
        ../third_party/physx/physx/source/geomutils/src/hf
        ../third_party/physx/physx/source/geomutils/src/pcm
        ../third_party/physx/physx/source/geomutils/src/ccd
        ../third_party/physx/physx/source/physxgpu/include
        ../third_party/physx/physx/../pxshared/include
        ../third_party/physx/physx/source/foundation/include)

set(PROJECT_FILES
        ${FILES}
        )

add_library(${PROJECT_NAME} STATIC ${PROJECT_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC ../
        ${CMAKE_CURRENT_SOURCE_DIR} ${PHYSX_SOURCE})


