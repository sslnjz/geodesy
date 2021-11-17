cmake_minimum_required(VERSION 3.6)
project(dms)

set(SRC
    dms.cpp
)

add_library(${PROJECT_NAME} OBJECT ${SRC})