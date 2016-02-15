SOURCE=odeint.h odeint.c nrutil.h nrutil.c stiff.h stiff.c stifbs.h stifbs.c derivs.h derivs.c oiutil.h oiutil.c simulation.c simulation.h
FLAGS=-Wno-unused-variable -std=c99
DEBUG=-g

all: clean
	gcc $(DEBUG) $(SOURCE) main.c -o main.exe -Wall $(FLAGS)
	
test: clean
	gcc $(DEBUG) $(SOURCE) test.h test.c -o test.exe -Wall $(FLAGS)
	
clean:
	del *.exe
	del *.csv