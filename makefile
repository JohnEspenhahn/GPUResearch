SOURCE=odeint.h odeint.c nrutil.h nrutil.c main.c
FLAGS=-Wno-unused-variable -std=c99
DEBUG=-g

all:
	gcc $(DEBUG) $(SOURCE) -o main.exe -Wall $(FLAGS)