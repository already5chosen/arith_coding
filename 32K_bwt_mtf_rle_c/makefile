CC = clang++ -std=c++11 -O2 -Wall -DENABLE_PERF_COUNTx

all: bwt_ari_encode.exe bwt_ari_decode.exe

bwt_sort.o: bwt_sort.cpp bwt_sort.h
	${CC} -c $<

bwt_mtf_rle_e.o: bwt_mtf_rle_e.cpp bwt_mtf_rle_e.h
	${CC} -c $<

bwt_mtf_rle_ari_e_main.o: bwt_mtf_rle_ari_e_main.cpp bwt_sort.h bwt_mtf_rle_e.h arithmetic_encode.h
	${CC} -c $<

arithmetic_encode.o: arithmetic_encode.cpp arithmetic_encode.h arithmetic_coder_ut.h arithmetic_coder_cfg.h
	${CC} -c $<

arithmetic_coder_ut.o: arithmetic_coder_ut.cpp arithmetic_coder_ut.h
	${CC} -c $<

bwt_ari_encode.exe: bwt_mtf_rle_ari_e_main.o bwt_sort.o bwt_mtf_rle_e.o arithmetic_encode.o arithmetic_coder_ut.o
	${CC} $+ -o $@

arithmetic_decode.o: arithmetic_decode.cpp arithmetic_decode.h arithmetic_coder_ut.h arithmetic_coder_cfg.h
	${CC} -c $<

bwt_d.o: bwt_d.cpp bwt_d.h
	${CC} -c $<

bwt_mtf_rle_ari_d_main.o: bwt_mtf_rle_ari_d_main.cpp  arithmetic_decode.h bwt_d.h
	${CC} -c $<

bwt_ari_decode.exe: bwt_mtf_rle_ari_d_main.o bwt_d.o arithmetic_decode.o arithmetic_coder_ut.o
	${CC} $+ -o $@
