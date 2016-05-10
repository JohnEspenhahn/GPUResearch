PAPER=g10
OBJECTS=odeint.o nrutil.o stiff.o $(PAPER)/derivs.o oiutil.o $(PAPER)/simulation.o jacobian.o $(PAPER)/tscalecsv.o
CFLAGS=-Wno-unused-variable -std=c99 -g -Wall

.c.o:
	gcc -c $*.c $(CFLAGS)

main.exe: main.o $(OBJECTS)
	gcc -o $@ $(**F) $(CFLAGS)
	del *.o
	del *.gch
	
clean:
	del *.exe
	del *.csv
	del *.o
	del *.gch