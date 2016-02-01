SOURCE=odeint.h odeint.c nrutil.h nrutil.c derivs.h derivs.c
FLAGS=-Wno-unused-variable -std=c99
DEBUG=-g

all:
	gcc $(DEBUG) $(SOURCE) main.c -o main.exe -Wall $(FLAGS)
	main.exe
	
test:
	gcc $(DEBUG) $(SOURCE) test.c -o debugging.exe -Wall $(FLAGS)