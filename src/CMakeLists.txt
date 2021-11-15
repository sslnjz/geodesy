cmake_minimum_required(VERSION 3.6)

include_directories(${OPENSSL_INCLUDE_DIR} ${INC_DIR}/geodesy)

add_subdirectory(dms)
add_subdirectory(utils)
add_subdirectory(latlon)

add_dependencies(dms LINK_HEADERS)
add_dependencies(utils LINK_HEADERS)
add_dependencies(latlon LINK_HEADERS)

add_library(
	${PROJECT_NAME} STATIC
	$<TARGET_OBJECTS:dms>
	$<TARGET_OBJECTS:utils>
	$<TARGET_OBJECTS:latlon>
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