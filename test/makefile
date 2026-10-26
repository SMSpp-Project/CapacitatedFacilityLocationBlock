##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for CapacitatedFacilityLocationBlock)                  #
#                                                                            #
#                                Donato Meoli                                #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# module name
NAME = CapacitatedFacilityLocationBlock_test

# basic directory
DIR = .

# compiler
CC = clang++

# common flags
COMMON_SW = -std=c++20 $(SMSppSWERR) -Wno-enum-compare

# debug switches
SW_DEBUG = -g3 $(SMSppSWDBG) -fno-inline $(COMMON_SW)

# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG

# default target - build in release mode
default: SW = $(SW_RELEASE)
default: $(DIR)/$(NAME)

# debug target - build in debug mode
debug: SW = $(SW_DEBUG)
debug: $(DIR)/$(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -

# CapacitatedFacilityLocationBlock (and SMS++ with it)
CFLBkSDR = ..
include $(CFLBkSDR)/makefile-c

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(CFLBkOBJ)

# libraries
MLIB =  $(CFLBkLIB)

$(DIR)/$(NAME): $(MOBJ) $(DIR)/test.o
	$(CC) -o $(DIR)/$(NAME) $(DIR)/test.o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include directives
MINC =  $(CFLBkINC)

# includes
MH =    $(CFLBkH)

# compile command
$(DIR)/test.o: $(DIR)/test.cpp $(MH)
	$(CC) -c $(DIR)/test.cpp -o $@ $(MINC) $(SW)

# clean target
clean::
	rm -f $(DIR)/*.o $(DIR)/*~
	rm -f $(DIR)/$(NAME)

# phony targets
.PHONY: default debug clean

############################ End of makefile #################################
