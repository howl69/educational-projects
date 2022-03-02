all: prog

prog: main.o file1.o file2.o
	gcc -o prog main.o file1.o file2.o -m32
main.o: main.c
	gcc -c -o main.o main.c -m32
file2.o: file2.asm
	nasm -f elf32 -o file2.o file2.asm
file1.o: file1.c
	gcc -c -o file1.o file1.c -m32
	
clean:
	rm -f *.o prog
