all: clean
	gcc main.c read.c mmio.c -o main.out

read:clean
	gcc read.c mmio.c -o read.out


example:clean
	gcc read_example.c mmio.c -o read.out
	
clean: 
	rm -f *.out