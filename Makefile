CC = gcc  
CFLAGS += -std=gnu11 -Wall -fopenmp -g3 -DNDEBUG 
LDFLAGS += -lm -fopenmp  
 

all: clean
	$(CC) $(CFLAGS)  main.c read.c mmio.c -o main.out $(LDFLAGS)
	./main.out

omp: clean
	$(CC) $(CFLAGS)  main_omp.c read.c mmio.c queue.c -o main_omp.out $(LDFLAGS)
	./main_omp.out
read:clean
	gcc -g read.c mmio.c -o read.out


example:clean
	gcc read_example.c mmio.c -o read.out
	
clean: 
	rm -f *.out