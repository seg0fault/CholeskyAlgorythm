
#ASAN STRICT COMPILER
#FLAGS = -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format -fsanitize=address

#NO ASAN STRICT COMPILER
FLAGS = -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

#LEAK CHECK
#FLAGS = -Wall -Wextra -pedantic -Werror -lm -fsanitize=leak -g

#OPTIMISATION
#FLAGS = -O3 -Wall -Wextra -pedantic -Werror -lm -g

all: a.out
	rm -rf *.o
a.out: main.o matrix.o block.o block_calc.o matrix_calc.o
	g++ $(FLAGS) $^
main.o: main.cpp matrix.h
	g++ $(FLAGS) -c $<
matrix.o: matrix.cpp matrix.h
	g++ $(FLAGS) -c $<
block_calc.o: block_calc.cpp block.h
	g++ $(FLAGS) -c $<	
block.o: block.cpp block.h
	g++ $(FLAGS) -c $<
matrix_calc.o: matrix_calc.cpp matrix.h block.h
	g++ $(FLAGS) -c $<
clean: a.out
	rm -rf *.o
