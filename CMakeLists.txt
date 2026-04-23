# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ CapacitatedFacilityLocationBlock               #
#                                                                             #
#    This file allows one to build the library using CMake.                   #
#    To do so, you can use the following commands:                            #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#                                                                             #
#    The following command also installs the library in the system:           #
#                                                                             #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.15)
project(CapacitatedFacilityLocationBlock
        VERSION 0.2.0
        DESCRIPTION "SMS++ CapacitatedFacilityLocationBlock module"
        HOMEPAGE_URL https://gitlab.com/smspp/capacitatedfacilitylocationblock
        LANGUAGES C CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})
set(modNamespace "SMS++")

# Find out if it's being called by the umbrella.
get_directory_property(hasParent PARENT_DIRECTORY)

# This adds the cmake directory to the module search paths,
# allowing us to use our modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Settings ------------------------------------------------------------ #
include(CTest)
include(BuildType)

# ----- Requirements -------------------------------------------------------- #
find_package(Threads QUIET)

if (NOT hasParent)
    find_package(MCFBlock REQUIRED)
    find_package(BinaryKnapsackBlock REQUIRED)
    find_package(StochasticBlock REQUIRED)
endif ()

# ----- Configuration header ------------------------------------------------ #
configure_file(cmake/${modName}Config.h.in ${modName}Config.h)

# ----- Library ------------------------------------------------------------- #
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_17)

if (MSVC AND BUILD_SHARED_LIBS)
    target_link_options(${modName} INTERFACE "/INCLUDE:SMSpp_force_load_${modName}")
endif()

target_sources(${modName} PRIVATE src/CapacitatedFacilityLocationBlock.cpp
                                  src/ScenarioReductionSolver.cpp
                                  src/CSSCScenarioReductionSolver.cpp)

target_include_directories(
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/${modNamespace}>)

# StochasticBlock headers needed privately to compile ScenarioReductionSolver.cpp
# (CSSC algorithm uses DiscreteScenarioSet and StochasticBlock internally).
# PRIVATE: only needed to compile this library, not exported to downstream users.
if (hasParent)
    # When built from umbrella: use the sibling source directory directly
    target_include_directories(${modName} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../StochasticBlock/include)
    target_link_libraries(${modName} PRIVATE StochasticBlock)
else ()
    # When built standalone: StochasticBlock was found via find_package above
    target_link_libraries(${modName} PRIVATE ${modNamespace}::StochasticBlock)
endif ()

target_link_libraries(${modName} PUBLIC
                      ${modNamespace}::MCFBlock
                      ${modNamespace}::BinaryKnapsackBlock
                      ${CMAKE_THREAD_LIBS_INIT})

add_library(${modNamespace}::${modName} ALIAS ${modName})

# ----- Subdirectories ------------------------------------------------------ #
add_subdirectory(tools)

if (BUILD_TESTING)
    add_subdirectory(test)
endif ()

# ----- Custom commands and targets ----------------------------------------- #

set(TXT_URL "https://gitlab.com/api/v4/projects/24780184/packages/generic/txt/latest/txt.tgz")
set(TXT_ARCHIVE ${CMAKE_CURRENT_SOURCE_DIR}/data/txt.tgz)
set(TXT_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/download_txt.cmake)

file(WRITE ${TXT_SCRIPT} "
file(DOWNLOAD \"${TXT_URL}\" \"${TXT_ARCHIVE}\" SHOW_PROGRESS)
")

add_custom_command(
        OUTPUT ${TXT_ARCHIVE}
        COMMAND ${CMAKE_COMMAND} -E echo "Downloading txt.tgz from GitLab Package Registry..."
        COMMAND ${CMAKE_COMMAND} -P ${TXT_SCRIPT}
        COMMENT "Downloading txt.tgz"
        VERBATIM
)

add_custom_target(download_txt
        DEPENDS ${TXT_ARCHIVE})

if (UNIX)
    set(TAR_FLAGS --warning=no-unknown-keyword -xzf)
elseif (WIN32)
    set(TAR_FLAGS xzf)
endif ()

add_custom_command(
        OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/data/.extracted
        COMMAND tar ${TAR_FLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/data/txt.tgz
        COMMAND echo > ${CMAKE_CURRENT_SOURCE_DIR}/data/.extracted
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data
        COMMENT "Extracting txt.tgz"
        VERBATIM
        DEPENDS download_txt)

add_custom_target(extract_txt
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/data/.extracted)

add_dependencies(extract_txt download_txt)

if (UNIX)
    add_custom_command(
            OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/data/nc4/.generated
            COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tools/batch ${CMAKE_CURRENT_BINARY_DIR}/tools/txt2nc4
            COMMAND echo > ${CMAKE_CURRENT_SOURCE_DIR}/data/nc4/.generated
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools
            COMMENT "Running txt2nc4"
            VERBATIM
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/tools/txt2nc4 ${CMAKE_CURRENT_SOURCE_DIR}/data/.extracted)
elseif (WIN32)
    add_custom_command(
            OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/data/nc4/.generated
            COMMAND "C:/Program Files/Git/bin/bash.exe"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tools/batch" "$<TARGET_FILE:txt2nc4>"
            COMMAND echo > ${CMAKE_CURRENT_SOURCE_DIR}/data/nc4/.generated
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools
            COMMENT "Running txt2nc4"
            VERBATIM
            DEPENDS $<TARGET_FILE:txt2nc4> ${CMAKE_CURRENT_SOURCE_DIR}/data/.extracted)
endif ()

add_custom_target(run_txt2nc4
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/data/nc4/.generated)

add_dependencies(run_txt2nc4 txt2nc4)

# ----- Install instructions ------------------------------------------------ #
include(GNUInstallDirs)

install(TARGETS ${modName}
        EXPORT ${modName}Targets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/CapacitatedFacilityLocationBlock.h
              ${CMAKE_CURRENT_SOURCE_DIR}/include/ScenarioReductionSolver.h
              ${CMAKE_CURRENT_SOURCE_DIR}/include/CSSCScenarioReductionSolver.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

install(EXPORT ${modName}Targets
        NAMESPACE ${modNamespace}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/${modName}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/README.md
        ${CMAKE_CURRENT_LIST_DIR}/LICENSE
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${modName})

# ----- Add the build tree to RPATH ----------------------------------------- #
if (hasParent)
    set(SMSPP_BUILD_RPATH
        "${SMSPP_BUILD_RPATH};${CMAKE_CURRENT_BINARY_DIR}"
        PARENT_SCOPE)
endif ()

# ----- Packaging support --------------------------------------------------- #
if (NOT hasParent)
    set(CPACK_PACKAGE_VENDOR "SMS++ Team")
    set(CPACK_PACKAGE_CONTACT "Donato Meoli")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_LIST_DIR}/LICENSE")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_LIST_DIR}/README.md")
    set(CPACK_SOURCE_IGNORE_FILES
        \\\\.git/
        \\\\.gitignore$
        \\\\.gitattributes$)
    include(CPack)
endif ()

# --------------------------------------------------------------------------- #
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# --------------------------------------------------------------------------- #