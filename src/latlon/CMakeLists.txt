cmake_minimum_required(VERSION 3.6)
project(latlon)

set(SRC
	"latlon_spherical.cpp"
	"latlon_ellipsoidal.cpp"
)

add_library(${PROJECT_NAME} OBJECT ${SRC})