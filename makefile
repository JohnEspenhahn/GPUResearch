PAPER=g10
SOURCE=odeint.h odeint.c nrutil.h nrutil.c stiff.h stiff.c stifbs.h stifbs.c $(PAPER)/derivs.h $(PAPER)/derivs.c oiutil.h oiutil.c $(PAPER)/simulation.c simulation.h jacobian.h jacobian.c
FLAGS=-Wno-unused-variable -std=c99
DEBUG=-g

all: clean
	gcc $(DEBUG) $(SOURCE) main.c -o main.exe -Wall $(FLAGS)
	
test: clean
	gcc $(DEBUG) $(SOURCE) test.h test.c -o test.exe -Wall $(FLAGS)
	
clean:
	del *.exe
	del *.csv