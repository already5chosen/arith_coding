CC = g++ -O2 -Wall -DENABLE_PERF_COUNTx

all: bwt_ari_encode.exe bwt_ari_decode.exe

bwt_mtf_rle_e.o: bwt_mtf_rle_e.cpp bwt_mtf_rle_e.h
	${CC} -c $<

bwt_mtf_rle_ari_e_main.o: bwt_mtf_rle_ari_e_main.cpp bwt_mtf_rle_e.h arithmetic_encode.h
	${CC} -c $<

bwt_ari_encode.exe: bwt_mtf_rle_ari_e_main.o bwt_mtf_rle_e.o
	${CC} $+ -o $@

bwt_mtf_rle_d.o: bwt_mtf_rle_d.cpp bwt_mtf_rle_d.h
	${CC} -c $<

bwt_mtf_rle_ari_d_main.o: bwt_mtf_rle_ari_d_main.cpp bwt_mtf_rle_d.h
	${CC} -c $<

bwt_ari_decode.exe: bwt_mtf_rle_ari_d_main.o bwt_mtf_rle_d.o
	${CC} $+ -o $@
