## CMakeLists.txt - CMake file that defines the build for the app folder, works in conjunction with the main CMakeLists.txt
 # written by C. D. Degawa, last modified in 2020 - see License.htm for legal notices
 #
 # The copyright in this software is being made available under the exhale Copyright License
 # and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 # party rights, including patent rights. No such rights are granted under this License.
 #
 # Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 ##

add_executable(exhaleApp
    exhaleAppPch.h
    loudnessEstim.cpp
    basicMP4Writer.cpp
    basicMP4Writer.h
    exhaleApp.cpp
    loudnessEstim.h
    exhaleApp.ico
    exhaleApp.rc
    basicWavReader.h
    basicWavReader.cpp
    exhaleAppPch.cpp
    ${PROJECT_SOURCE_DIR}/include/exhaleDecl.h
    ${PROJECT_SOURCE_DIR}/include/version.h)

set_target_properties(exhaleApp PROPERTIES OUTPUT_NAME exhale)

if(TARGET Threads::Threads)
    target_link_libraries(exhaleApp PRIVATE Threads::Threads)
endif()
if(CMAKE_DL_LIBS)
    target_link_libraries(exhaleApp PRIVATE ${CMAKE_DL_LIBS})
endif()
target_link_libraries(exhaleApp PRIVATE exhaleLib)
target_include_directories(exhaleApp PRIVATE ${PROJECT_SOURCE_DIR}/include)
if(BUILD_SHARED_LIBS)
    target_compile_definitions(exhaleApp PRIVATE USE_EXHALELIB_DLL)
endif(BUILD_SHARED_LIBS)

# PCH requires at least 3.16
# I actually don't know if this works or not
if(CMAKE_VERSION VERSION_GREATER "3.16.0")
    target_precompile_headers(exhaleApp PUBLIC ${PROJECT_SOURCE_DIR}/src/lib/exhaleLibPch.h)
endif()

install(TARGETS exhaleApp
    RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
