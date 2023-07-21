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

# Specifies the build type on single-configuration generators (e.g. Makefile Generators or Ninja). 
# Typical values include Debug, Release, RelWithDebInfo and MinSizeRel, but custom build types can 
# also be defined.
set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "build type")

# Do not add run time path information
set(CMAKE_SKIP_RPATH TRUE)

project(
	geodesy
	VERSION 0.0.1
	LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 17)

# compile platform check
if (CYGWIN)
	message(FATAL_ERROR "Sorry, DO NOT support Cygwin")
endif ()

if (MINGW)
	message(FATAL_ERROR "Sorry, DO NOT support MinGW")
endif ()

include(GNUInstallDirs)

set(CMAKE_CONFIG_INSTALL_FILE ${PROJECT_BINARY_DIR}/config.toinstall.cmake)
set(CMAKE_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
set(INC_DIR ${PROJECT_SOURCE_DIR}/include CACHE PATH "geodesy inc")
set(LIB_DIR ${PROJECT_SOURCE_DIR}/lib CACHE PATH "geodesy lib")

message(${INC_DIR})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIR})

add_custom_target(
	LINK_HEADERS ALL
	COMMENT "link headers..."
)

INCLUDE(CMakeLists_Headers.txt)

macro(makeLink src dest target)
	add_custom_command(
		TARGET ${target} PRE_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different ${src} ${dest}
		DEPENDS ${dest}
	)
endmacro()

add_custom_command(
	TARGET LINK_HEADERS PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory ${INC_DIR}/${PROJECT_NAME}
)

foreach(header_file ${INCLUDE_HEADERS})
	string(REPLACE "/" ";" arr ${header_file})
	list(GET arr -1 file_name)
	makeLink(${PROJECT_SOURCE_DIR}/${header_file} ${INC_DIR}/${PROJECT_NAME}/${file_name} LINK_HEADERS)
endforeach()

message("CMAKE_C_FLAGS_DEBUG is ${CMAKE_C_FLAGS_DEBUG}")
message("CMAKE_C_FLAGS_RELEASE is ${CMAKE_C_FLAGS_RELEASE}")
message("CMAKE_C_FLAGS_RELWITHDEBINFO is ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
message("CMAKE_C_FLAGS_MINSIZEREL is ${CMAKE_C_FLAGS_MINSIZEREL}")

message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
message("CMAKE_CXX_FLAGS_RELWITHDEBINFO is ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message("CMAKE_CXX_FLAGS_MINSIZEREL is ${CMAKE_CXX_FLAGS_MINSIZEREL}")

if (WIN32)
	set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   /MP /wd4200")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /wd4200 /std:c++17")
else ()
	set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Wall -fPIC -pipe -std=gnu90")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -pipe -std=c++17 -Wno-invalid-offsetof -Wno-overloaded-virtual")
endif ()

add_subdirectory(src)
#add_subdirectory(test)
add_subdirectory(ui)

####CONFIG
include(CMakePackageConfigHelpers)
set(CONFIG_INC_DIR ${INC_DIR})
set(CONFIG_LIB_DIR ${LIB_DIR})
configure_package_config_file(
	${PROJECT_NAME}-config.cmake.in
	${PROJECT_SOURCE_DIR}/${PROJECT_NAME}-config.cmake
	INSTALL_DESTINATION ${CMAKE_CONFIG_INSTALL_DIR}
	PATH_VARS CONFIG_INC_DIR CONFIG_LIB_DIR
)

set(CONFIG_INC_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(CONFIG_LIB_DIR ${CMAKE_INSTALL_LIBDIR})
configure_package_config_file(
	${PROJECT_NAME}-config.cmake.in
	${CMAKE_CONFIG_INSTALL_FILE}
	INSTALL_DESTINATION ${CMAKE_CONFIG_INSTALL_DIR}
	PATH_VARS CONFIG_INC_DIR CONFIG_LIB_DIR
)

install(
	FILES ${CMAKE_CONFIG_INSTALL_FILE}
	DESTINATION ${CMAKE_CONFIG_INSTALL_DIR}
	COMPONENT devel
	RENAME ${PROJECT_NAME}-config.cmake
)

install(
	FILES ${INCLUDE_HEADERS} ${INCLUDE_KERNEL_HEADERS}
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
	COMPONENT devel
)

install(
	FILES README.md
	DESTINATION "${CMAKE_INSTALL_DOCDIR}-${PROJECT_VERSION}"
	COMPONENT devel
)
