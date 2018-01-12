CC = clang++ -std=c++11 -O2 -Wall -DENABLE_PERF_COUNTx

all: ari_encode.exe ari_decode.exe

arithmetic_encode.o: arithmetic_encode.cpp arithmetic_encode.h
	${CC} -c arithmetic_encode.cpp

arithmetic_encode_main.o: arithmetic_encode_main.cpp arithmetic_encode.h
	${CC} -c arithmetic_encode_main.cpp

ari_encode.exe: arithmetic_encode_main.o arithmetic_encode.o
	${CC} arithmetic_encode_main.o arithmetic_encode.o -o $@

arithmetic_decode.o: arithmetic_decode.cpp arithmetic_decode.h
	${CC} -c arithmetic_decode.cpp

arithmetic_decode_main.o: arithmetic_decode_main.cpp arithmetic_decode.h
	${CC} -c arithmetic_decode_main.cpp

ari_decode.exe: arithmetic_decode_main.o arithmetic_decode.o
	${CC} arithmetic_decode_main.o arithmetic_decode.o -o $@
