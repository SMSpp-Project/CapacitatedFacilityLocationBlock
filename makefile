##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of CapacitatedFacilityLocationBlock                             #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by CapacitatedFacilityLocationBlock, i.e., core SMS++,  #
#   BinaryKnapsackBlock and MCFBlock. These are *not* copied into            #
#   $(CFLBkINC): adding those -I directives to the compile commands will     #
#   have to done by whatever "main" makefile is using this. Analogously,     #
#   any external library and the corresponding -L< libdirs > will have to    #
#   be added to the final linking command by  whatever "main" makefile is    #
#   using this.                                                              #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by CapacitatedFacilityLocationBlock.                                #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(CFLBkSDR)    = the directory where the source is               #
#           $(BKBkOBJ)     = the final object for BinaryKnapsackBlock        #
#           $(BKBkINC)     = the -I$( BinaryKnapsackBlock source directory ) #
#           $(MCFBkOBJ)    = the final object for MCFBlock                   #
#           $(MCFBkINC)    = the -I$( MCFBlock source directory )            #
#                                                                            #
#   Output: $(CFLBkOBJ)    = the final object(s) / library                   #
#           $(CFLBkH)      = the .h files to include                         #
#           $(CFLBkINC)    = the -I$( source directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################


# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

CFLBkOBJ = $(CFLBkSDR)/CapacitatedFacilityLocationBlock.o

CFLBkINC = -I$(CFLBkSDR)

CFLBkH   = $(CFLBkSDR)/CapacitatedFacilityLocationBlock.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(CFLBkOBJ) $(CFLBkSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(CFLBkSDR)/CapacitatedFacilityLocationBlock.o: $(BKBkOBJ) $(MCFBkOBJ) \
	$(CFLBkSDR)/CapacitatedFacilityLocationBlock.cpp $(SMS++OBJ)
	$(CC) -c $(CFLBkSDR)/CapacitatedFacilityLocationBlock.cpp -o $@ \
	$(CFLBkINC) $(BKBkINC) $(MCFBkINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
