## CMakeLists.txt - CMake file that defines the build for the lib folder, works in conjunction with the main CMakeLists.txt
 # written by C. D. Degawa, last modified in 2020 - see License.htm for legal notices
 #
 # The copyright in this software is being made available under the exhale Copyright License
 # and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 # party rights, including patent rights. No such rights are granted under this License.
 #
 # Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 ##

add_library(exhaleLib
    lappedTransform.cpp
    exhaleLibPch.cpp
    bitStreamWriter.cpp
    quantization.cpp
    stereoProcessing.h
    exhaleLibPch.h
    entropyCoding.cpp
    tempAnalysis.cpp
    bitAllocation.cpp
    stereoProcessing.cpp
    bitAllocation.h
    bitStreamWriter.h
    specAnalysis.h
    specAnalysis.cpp
    lappedTransform.h
    specGapFilling.cpp
    specGapFilling.h
    linearPrediction.h
    quantization.h
    entropyCoding.h
    exhaleEnc.cpp
    tempAnalysis.h
    linearPrediction.cpp
    exhaleEnc.h
    ${PROJECT_SOURCE_DIR}/include/exhaleDecl.h
    ${PROJECT_SOURCE_DIR}/include/version.h)

set_target_properties(exhaleLib PROPERTIES OUTPUT_NAME exhale)

if(TARGET Threads::Threads)
    target_link_libraries(exhaleLib PRIVATE Threads::Threads)
endif()
if(CMAKE_DL_LIBS)
    target_link_libraries(exhaleLib PRIVATE ${CMAKE_DL_LIBS})
endif()
target_include_directories(exhaleLib PRIVATE ${PROJECT_SOURCE_DIR}/include)
set_target_properties(exhaleLib PROPERTIES PUBLIC_HEADER "${PROJECT_SOURCE_DIR}/include/exhaleDecl.h")

# PCH requires at least 3.16
# I actually don't know if this works or not
if(CMAKE_VERSION VERSION_GREATER "3.16.0")
    target_precompile_headers(exhaleLib PUBLIC ${PROJECT_SOURCE_DIR}/src/lib/exhaleLibPch.h)
endif()

target_compile_definitions(exhaleLib PRIVATE EXHALE_DYN_LINK)

install(TARGETS exhaleLib
    ARCHIVE DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_FULL_INCLUDEDIR})
