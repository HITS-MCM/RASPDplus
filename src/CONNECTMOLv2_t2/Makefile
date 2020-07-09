
INCLUDE=include
SRC=src

OBJECTS= Connect.o Pdb.o  PdbRecord.o Ring.o  

#ALL_WARNINGS= -Werror -Wall -W -Wstrict-prototypes -Wmissing-prototypes -Wpointer-arith -Wreturn-type -Wcast-qual\
	-Wwrite-strings -Wswitch -Wshadow -Wcast-align -Wuninitialized -Wbad-function-cast -Wchar-subscripts\
	-Winline -Wnested-externs -Wredundant-decls -ansi -pedantic -fbounds-check

#ALL_WARNINGS= -Wall -W -Wstrict-prototypes -Wmissing-prototypes -Wpointer-arith -Wreturn-type -Wcast-qual\
	-Wwrite-strings -Wswitch -Wshadow -Wcast-align  -Wchar-subscripts -Winline -Wredundant-decls -ansi \
	-pedantic -fbounds-check

#FLAGS=-fast  -`fpversion -foption` -mt -DNDEBUG
#FLAGS= -g -pg -I$(INCLUDE) -Wall

FLAGS= -I$(INCLUDE) -O3 -Wall 
#FLAGS= -I$(INCLUDE) -O3 -Wall -std=c99 
#FLAGS= -I$(INCLUDE) -O3 -Wall -std=c++0x 
#FLAGS= -I$(INCLUDE) -O3 -Wall -std=c++98 -ansi
#FLAGS= -I$(INCLUDE) -O3 -Wall -ansi
#FLAGS= -I$(INCLUDE) -O3 -Wall -std=gnu++98 -ansi

#FLAGS= -I$(INCLUDE) -pg
#FLAGS= $(ALL_WARNINGS)
CC=g++


#Connect.o:  $(INCLUDE)/Pdb.hpp $(INCLUDE)/constants.hpp $(INCLUDE)/Ring.h Connect.cpp
#	$(CC)  $(FLAGS) -c -I$(INCLUDE) Connect.cpp 

Connect2.0.exe:   $(OBJECTS)
	$(CC) $(FLAGS) -I$(INCLUDE) $(OBJECTS) -o Connect2.0.exe 

Connect.o: $(INCLUDE)/Pdb.hpp $(INCLUDE)/constants.hpp $(INCLUDE)/Ring.h Connect.cpp
	$(CC)  $(FLAGS) -c -I$(INCLUDE) Connect.cpp
       	
Pdb.o:  $(INCLUDE)/Pdb.hpp $(INCLUDE)/constants.hpp $(SRC)/Pdb.cpp
	$(CC)  $(FLAGS) -c -I$(INCLUDE) $(SRC)/Pdb.cpp

PdbRecord.o: $(INCLUDE)/PdbRecord.hpp $(INCLUDE)/constants.hpp $(INCLUDE)/Point_3d.hpp $(SRC)/PdbRecord.cpp
	$(CC) $(FLAGS) -c -I$(INCLUDE) $(SRC)/PdbRecord.cpp 

Ring.o:$(INCLUDE)/Ring.h $(INCLUDE)/Pdb.hpp $(SRC)/Ring.cc
	$(CC) $(FLAGS) -c -I$(INCLUDE) $(SRC)/Ring.cc 
	
#MD.o: $(INCLUDE)/MD.hpp $(INCLUDE)/Point_3d.hpp $(INCLUDE)/constants.hpp $(SRC)/MD.cpp
#	$(CC) $(FLAGS) -c -I$(INCLUDE) $(SRC)/MD.cpp  

#mymath.o: $(INCLUDE)/mymath.hpp $(SRC)/mymath.cpp
#	$(CC) $(FLAGS) -c -I$(INCLUDE) $(SRC)/mymath.cpp  

clean:
	\rm $(OBJECTS) Connect2.0.exe 

