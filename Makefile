# Compiler and Compilerflags
CC = gcc
CFLAGS = -O2 -Wall -Wextra -Werror -march=native

EXEDIR = $(if ${MyLocalPath}, ${MyLocalPath}, bin)
EXE = $(EXEDIR)/morse-potential

LIB = -lm `pkg-config --cflags --libs gsl`

OBJ = MorsePotential.o

all: $(EXE)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIB) -o $@

clean:
	rm -f $(OBJ) $(EXE)
