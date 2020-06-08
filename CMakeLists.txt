cmake_minimum_required(VERSION 3.0.0)
project(TinyRenderer VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include(CTest)
enable_testing()

file(GLOB SOURCE_FILES "src/**")

add_executable(TinyRenderer ${SOURCE_FILES})

include_directories(include)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
