# --------------------------------------------------------------------------- #
#    CMake file for CapacitatedFacilityLocationBlock/test_SRsolver            #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# Since we are using the block factory, objects from linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.
if (BUILD_SHARED_LIBS)
    if (UNIX AND (NOT APPLE))
        add_link_options("-Wl,--no-as-needed")
    endif ()
else ()
    if (MSVC)
        add_link_options("/WHOLEARCHIVE" "/FORCE:MULTIPLE")
    else () # Unix
        if (APPLE)
            add_link_options("-Wl,-all_load")
        else ()
            add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
        endif ()
    endif ()
endif ()

# ----- CapacitatedFacilityLocationBlock_test_SRsolver ---------------------- #
add_executable(test_SRsolver test_SRsolver.cpp)
target_link_libraries(test_SRsolver PRIVATE
                      SMS++::CapacitatedFacilityLocationBlock)

if (WIN32 AND BUILD_SHARED_LIBS AND NOT DEFINED ENV{CONDA_BUILD})
    # copy runtime DLLs next to the exe after build
    add_custom_command(TARGET test_SRsolver POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy_if_different
                       $<TARGET_RUNTIME_DLLS:test_SRsolver>
                       $<TARGET_FILE_DIR:test_SRsolver>
                       COMMAND_EXPAND_LISTS
                       COMMENT "Copying runtime DLLs next to test_SRsolver.exe")
endif()

add_test(NAME test_SRsolver
         COMMAND test_SRsolver
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# --------------------------------------------------------------------------- #
