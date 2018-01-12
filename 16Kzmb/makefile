CC = clang++ -std=c++11 -O2 -Wall -DENABLE_PERF_COUNTx

all: bwt_ari_encode.exe bwt_ari_decode.exe

bwt_e.o: bwt_e.cpp bwt_e.h
	${CC} -c bwt_e.cpp

arithmetic_encode.o: arithmetic_encode.cpp arithmetic_encode.h
	${CC} -c arithmetic_encode.cpp

bwt_arithmetic_encode_main.o: bwt_arithmetic_encode_main.cpp arithmetic_encode.h bwt_e.h
	${CC} -c bwt_arithmetic_encode_main.cpp

bwt_ari_encode.exe: bwt_arithmetic_encode_main.o arithmetic_encode.o bwt_e.o
	${CC} bwt_arithmetic_encode_main.o arithmetic_encode.o bwt_e.o -o $@

bwt_d.o: bwt_d.cpp bwt_d.h
	${CC} -c bwt_d.cpp

arithmetic_decode.o: arithmetic_decode.cpp arithmetic_decode.h
	${CC} -c arithmetic_decode.cpp

arithmetic_bwt_decode_main.o: arithmetic_bwt_decode_main.cpp arithmetic_decode.h bwt_d.h
	${CC} -c arithmetic_bwt_decode_main.cpp

bwt_ari_decode.exe: arithmetic_bwt_decode_main.o arithmetic_decode.o bwt_d.o
	${CC} arithmetic_bwt_decode_main.o arithmetic_decode.o bwt_d.o -o $@
