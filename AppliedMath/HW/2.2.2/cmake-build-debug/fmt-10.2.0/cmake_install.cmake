# Install script for directory: C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/2_2_2")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "C:/Program Files/JetBrains/CLion 2023.2.2/bin/mingw/bin/objdump.exe")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/cmake-build-debug/fmt-10.2.0/libfmtd.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/fmt" TYPE FILE FILES
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/args.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/chrono.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/color.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/compile.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/core.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/format.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/format-inl.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/os.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/ostream.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/printf.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/ranges.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/std.h"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/fmt-10.2.0/include/fmt/xchar.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt" TYPE FILE FILES
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/cmake-build-debug/fmt-10.2.0/fmt-config.cmake"
    "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/cmake-build-debug/fmt-10.2.0/fmt-config-version.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt/fmt-targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt/fmt-targets.cmake"
         "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/cmake-build-debug/fmt-10.2.0/CMakeFiles/Export/b834597d9b1628ff12ae4314c3a2e4b8/fmt-targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt/fmt-targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt/fmt-targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt" TYPE FILE FILES "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/cmake-build-debug/fmt-10.2.0/CMakeFiles/Export/b834597d9b1628ff12ae4314c3a2e4b8/fmt-targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/fmt" TYPE FILE FILES "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/cmake-build-debug/fmt-10.2.0/CMakeFiles/Export/b834597d9b1628ff12ae4314c3a2e4b8/fmt-targets-debug.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "C:/PHYSTECH/Books/AppliedMath/HW/2.2.2/cmake-build-debug/fmt-10.2.0/fmt.pc")
endif()
