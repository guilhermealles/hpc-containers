CC = gcc 
FC = gfortran -ffree-form
CFLAGS = -g -Wall -Werror
geol:	geol_reshape.c   
	$(CC) $(CFLAGS) $(COPTFLAGS) $(TIMER)  geol_reshape.c  $(CLDFLAGS) -o geol_reshape -lm
bin2ascii:   bin2ascii.c
	$(CC) $(CFLAGS) $(COPTFLAGS) $(TIMER)  bin2ascii.c  $(CLDFLAGS) -o bin2ascii -lm
parameter:   make_parameter.f90
	$(FC) $(CFLAGS) $(COPTFLAGS) $(TIMER)  make_parameter.f90  $(CLDFLAGS) -o make_parameter -lm



all: geol parameter bin2ascii
clean:
	rm -f *.o geol_reshape make_parameter bin2ascii
