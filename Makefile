all: clean
	gcc read_example.c mmio.c -o read.out

clean: 
	rm -f *.out