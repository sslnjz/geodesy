# Install script for directory: D:/workspace/github/geodesy

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/geodesy")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("D:/workspace/github/geodesy/build/src/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("D:/workspace/github/geodesy/build/test/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("D:/workspace/github/geodesy/build/ui/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/geodesy" TYPE FILE RENAME "geodesy-config.cmake" FILES "D:/workspace/github/geodesy/build/config.toinstall.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/geodesy" TYPE FILE FILES
    "D:/workspace/github/geodesy/src/dms/dms.h"
    "D:/workspace/github/geodesy/src/latlon/ellipsoids.h"
    "D:/workspace/github/geodesy/src/latlon/latlon.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_ellipsoidal.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_ellipsoidal_datum.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_ellipsoidal_referenceframe.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_ellipsoidal_referenceframe_txparams.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_ellipsoidal_vincenty.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_nvector_ellipsoidal.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_nvector_spherical.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_osgridref.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_spherical.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_utm.h"
    "D:/workspace/github/geodesy/src/latlon/latlon_utm_mgrs.h"
    "D:/workspace/github/geodesy/src/utils/algorithm.h"
    "D:/workspace/github/geodesy/src/utils/strutil.h"
    "D:/workspace/github/geodesy/src/utm/mgrs.h"
    "D:/workspace/github/geodesy/src/utm/osgridref.h"
    "D:/workspace/github/geodesy/src/utm/utm.h"
    "D:/workspace/github/geodesy/src/utm/utm_mgrs.h"
    "D:/workspace/github/geodesy/src/vector/cartesian.h"
    "D:/workspace/github/geodesy/src/vector/cartesian_datum.h"
    "D:/workspace/github/geodesy/src/vector/cartesian_referenceFrame.h"
    "D:/workspace/github/geodesy/src/vector/ned.h"
    "D:/workspace/github/geodesy/src/vector/nvector_cartesian.h"
    "D:/workspace/github/geodesy/src/vector/nvector_ellipsoidal.h"
    "D:/workspace/github/geodesy/src/vector/nvector_spherical.h"
    "D:/workspace/github/geodesy/src/vector/vector3d.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/geodesy-0.0.1" TYPE FILE FILES "D:/workspace/github/geodesy/README.md")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "D:/workspace/github/geodesy/build/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "D:/workspace/github/geodesy/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
