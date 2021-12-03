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

include_directories(${INC_DIR}/geodesy ${INC_DIR})

add_subdirectory(dms)
add_subdirectory(utils)
add_subdirectory(latlon)
add_subdirectory(vector)

add_dependencies(dms LINK_HEADERS)
add_dependencies(utils LINK_HEADERS)
add_dependencies(latlon LINK_HEADERS)
add_dependencies(vector LINK_HEADERS)

add_library(
	${PROJECT_NAME} STATIC
	$<TARGET_OBJECTS:dms>
	$<TARGET_OBJECTS:utils>
	$<TARGET_OBJECTS:latlon>
	$<TARGET_OBJECTS:vector>
)

install(
	TARGETS ${PROJECT_NAME}
	ARCHIVE
	DESTINATION ${CMAKE_INSTALL_LIBDIR}
	COMPONENT devel
)

if (APPLE)
	set(LIBSO ${LIB_DIR}/libgeodesy.a)
else ()
	set(LIBSO ${LIB_DIR}/libgeodesy.so)
	add_custom_target(
		SCRIPT_SHARED_LIB ALL
		COMMAND ${CMAKE_COMMAND} -E echo 'GROUP ( libgeodesy.a ) ' > ${LIBSO}
	)
	add_dependencies(SCRIPT_SHARED_LIB ${PROJECT_NAME})
endif ()

install(
	FILES ${LIBSO}
	DESTINATION ${CMAKE_INSTALL_LIBDIR}
	COMPONENT devel
)

install (
	TARGETS ${PROJECT_NAME} EXPORT ${CMAKE_PROJECT_NAME}-targets
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install (
	EXPORT ${CMAKE_PROJECT_NAME}-targets
	FILE ${CMAKE_PROJECT_NAME}-targets.cmake
	DESTINATION "${CMAKE_CONFIG_INSTALL_DIR}"
)