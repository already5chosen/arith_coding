CC = g++ -O2 -Wall -DENABLE_PERF_COUNTx

all: bwt_e.exe bwt_d.exe

bwt_e.o: bwt_e.cpp bwt_e.h
	${CC} -c bwt_e.cpp

bwt_e_main.o: bwt_e_main.cpp bwt_e.h
	${CC} -c bwt_e_main.cpp

bwt_e.exe: bwt_e_main.o bwt_e.o
	${CC} bwt_e_main.o bwt_e.o -o $@

bwt_d.o: bwt_d.cpp bwt_d.h
	${CC} -c bwt_d.cpp

bwt_d_main.o: bwt_d_main.cpp bwt_d.h
	${CC} -c bwt_d_main.cpp

bwt_d.exe: bwt_d_main.o bwt_d.o
	${CC} bwt_d_main.o bwt_d.o -o $@
