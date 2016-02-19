# Compiler
CC = gcc
# Compilerflags
CFLAGS = -g03 -Wall -Werror -Wstrict-prototypes -Wmissing-prototypes -mtune=native
# Zu erstellende fertige Programme (mit relativem Pfad vorangestellt)
EXE = ../../bin/MorsePotential
# Objektdateien (.o bzw. .out), Librarys und Abhängigkeiten
# (wenn Abhängigkeit geändert -> make wird neu ausgeführt)
OBJ1 = MorsePotential.o
DEP1 = MorsePotential.c Makefile
LIB1 = -lm `pkg-config --cflags --libs gsl`
# Kommandoblock "all:" als Einsprungspunkt für make
all: $(EXE)
#
$(EXE): $(OBJ1) $(DEP1)
	$(CC) $(CFLAGS) $(OBJ1) $(LIB1) -o $@

clean:
	rm $(OBJ1) $(EXE)
