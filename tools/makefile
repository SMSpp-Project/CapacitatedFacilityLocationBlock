##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of dmx2nc4                                                      #
#                                                                            #
#                                VERSION 1.00                                #
#                               25 - 10 - 2020                               #
#                                                                            #
#                              Implementation by:                            #
#                                                                            #
#                              Antonio Frangioni                             #
#                           Operations Research Group                        #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
# NAME can be main, test and dmx2nc4:, the corresponding .cpp is compiled
# and the final executable named $(NAME) is constructed 
NAME = dmx2nc4

# basic directory
DIR = ./

# switch to select the :MCFClass solver, supported values are
# -DHAVE_CSCL2      for the CS2 class
# -DHAVE_CPLEX      for the MCFCplex class
# -DHAVE_MFSMX      for the MCFSimplex class
# -DHAVE_MFZIB      for the MCFZIB class
# -DHAVE_RELAX      for the RelaxIV class
# -DHAVE_CPLEX      for the MCFCplex class
# -DHAVE_SPTRE      for the SPTree class; note that SPTree cannot solve
#                   most MCF instances, except those with SPT structure
MCFCS = -DHAVE_MFSMX

# debug switches
#SW = -g -glldb -fno-inline -std=c++17 -ferror-limit=1 $(MCFCS)
# production switches
SW = -O3 -std=c++17 -DNDEBUG $(MCFCS)

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)*.o $(DIR)*~ $(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# MCFBlock
MCFBkSDR = ../
include $(MCFBkSDR)makefile-c

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(MCFBkOBJ)

# libreries
MLIB =  $(MCFBkLIB)

$(NAME): $(MOBJ) $(DIR)$(NAME).o
	$(CC) -o $(NAME) $(DIR)$(NAME).o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

# include directives
MINC =  $(MCFBkINC)

# includes
MH =    $(MCFBkH)

# compile command

$(DIR)$(NAME).o: $(DIR)$(NAME).cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) $(SW)

############################ End of makefile #################################
