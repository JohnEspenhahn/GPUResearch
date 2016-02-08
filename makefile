SOURCE=odeint.h odeint.c nrutil.h nrutil.c bsstep.h bsstep.c stiff.h stiff.c derivs.h derivs.c 
FLAGS=-Wno-unused-variable -std=c99
DEBUG=-g

all: clean
	gcc $(DEBUG) $(SOURCE) main.c -o main.exe -Wall $(FLAGS)
	
test: clean
	gcc $(DEBUG) $(SOURCE) test.h test.c -o test.exe -Wall $(FLAGS)
	
clean:
	del *.exe
	del *.csv