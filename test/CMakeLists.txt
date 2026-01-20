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
    if (MSVC)
        add_link_options("/OPT:NOREF")
    else () # Unix
        if (APPLE)
            add_link_options("-Wl,-dead_strip_dylibs")
        else ()
            add_link_options("-Wl,--no-as-needed")
        endif ()
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

add_test(NAME test_SRsolver
         COMMAND test_SRsolver
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# --------------------------------------------------------------------------- #
