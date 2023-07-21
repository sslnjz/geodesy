# /**********************************************************************************
# *  MIT License                                                                    *
# *                                                                                 *
# *  Copyright (c) 2021 Binbin Song <ssln.jzs@gmail.com>                            *
# *                                                                                 *
# *  Geodesy tools for conversions between (historical) datums                      *
# *  (c) Chris Veness 2005-2019                                                     *
# *  www.movable-type.co.uk/scripts/latlong-convert-coords.html                     *
# *  www.movable-type.co.uk/scripts/geodesy-library.html#latlon-ellipsoidal-datum   *
# *                                                                                 *
# *  Permission is hereby granted, free of charge, to any person obtaining a copy   *
# *  of this software and associated documentation files (the "Software"), to deal  *
# *  in the Software without restriction, including without limitation the rights   *
# *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      *
# *  copies of the Software, and to permit persons to whom the Software is          *
# *  furnished to do so, subject to the following conditions:                       *
# *                                                                                 *
# *  The above copyright notice and this permission notice shall be included in all *
# *  copies or substantial portions of the Software.                                *
# *                                                                                 *
# *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
# *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       *
# *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    *
# *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         *
# *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  *
# *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  *
# *  SOFTWARE.                                                                      *
# ***********************************************************************************/
cmake_minimum_required(VERSION 3.20)

set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "build type")

project(geodesy_test
		LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 17)

find_package(geodesy REQUIRED CONFIG HINTS ..)
include_directories(${GEODESY_INCLUDE_DIR})
link_directories(${GEODESY_LIB_DIR})

################################
# GTest
################################
project(googletest-git NONE)

include(FetchContent)
FetchContent_Declare(
		googletest
		GIT_REPOSITORY https://github.com/google/googletest.git
		GIT_TAG        release-1.13.0
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

enable_testing()

find_program(CMAKE_MEMORYCHECK_COMMAND valgrind)
set(memcheck_command ${CMAKE_MEMORYCHECK_COMMAND} ${CMAKE_MEMORYCHECK_COMMAND_OPTIONS} --error-exitcode=1 --leak-check=full)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

if (WIN32)
		set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   /MP /wd4200")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /wd4200 /std:c++17")
else ()
		set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Wall -fPIC -pipe -std=gnu90 -Wunused-value")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -pipe -std=c++17 -Wno-overloaded-virtual -Wunused-value")
endif ()


set(TEST_LIST
		dms_unittest
		latlon_unittest
		latlon_spherical_unittest
		latlon_nvector_spherical_unittest
		latlon_ellipsoidal_unittest
		latlon_ellipsoidal_referenceframe_unittest
		latlon_ellipsoidal_datum_unittest
		latlon_ellipsoidal_vincenty_unittest
		vector3d_unittest
)

set(GEODESY_LIB geodesy)
set(GTEST_LIB gtest gtest_main)

foreach(src ${TEST_LIST})
	add_executable(${src} EXCLUDE_FROM_ALL ${src}.cpp)
	target_link_libraries(${src} ${GEODESY_LIB} ${GTEST_LIB})
	add_test(${src} ${src})
	add_dependencies(check ${src})
endforeach()

foreach(src ${TEST_LIST})
	add_test(${src}-memory-check ${memcheck_command} ./${src})
endforeach()