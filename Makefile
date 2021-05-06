all: clean
	gcc read.c mmio.c -o read

clean: 
	rm -f read