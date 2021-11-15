cmake_minimum_required(VERSION 3.6)
project(latlon)

set(SRC
	"latlon_spherical.cpp"
)

add_library(${PROJECT_NAME} OBJECT ${SRC})