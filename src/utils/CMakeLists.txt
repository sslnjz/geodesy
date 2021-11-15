cmake_minimum_required(VERSION 3.6)
project(utils)

set(SRC
	"strutil.cpp"
)

add_library(${PROJECT_NAME} OBJECT ${SRC})
