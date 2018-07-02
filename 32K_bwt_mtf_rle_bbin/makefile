CC = clang++ -std=c++11 -O2 -Wall
cc = clang  -std=c99 -O2 -Wall
# CC = g++ -std=c++11 -O2 -Wall
# cc = gcc  -std=c99 -O2 -Wall
UTIL_DIR := ../32K_bwt_mtf_rle_c

all: bwt_ari_encode.exe bwt_ari_decode.exe

bwt_sort.o: $(UTIL_DIR)/bwt_sort.cpp $(UTIL_DIR)/bwt_sort.h
	${CC} -c -I $(UTIL_DIR) $<

bwt_mtf_rle_e.o: bwt_mtf_rle_e.cpp bwt_mtf_rle_e.h
	${CC} -c $<

bwt_mtf_rle_ari_e_main.o: bwt_mtf_rle_ari_e_main.cpp $(UTIL_DIR)/bwt_sort.h bwt_mtf_rle_e.h arithmetic_encode.h
	${CC} -c -I $(UTIL_DIR) $<

arithmetic_encode.o: arithmetic_encode.cpp arithmetic_encode.h arithmetic_coder_cfg.h
	${CC} -c $<

bwt_ari_encode.exe: bwt_mtf_rle_ari_e_main.o bwt_sort.o bwt_mtf_rle_e.o arithmetic_encode.o
	${CC} $+ -o $@

arithmetic_decode.o: arithmetic_decode.cpp arithmetic_decode.h arithmetic_coder_cfg.h
	${CC} -c $<

bwt_d.o: $(UTIL_DIR)/bwt_d.cpp $(UTIL_DIR)/bwt_d.h
	${CC} -c -I $(UTIL_DIR) $<

bwt_mtf_rle_ari_d_main.o: bwt_mtf_rle_ari_d_main.cpp  arithmetic_decode.h $(UTIL_DIR)/bwt_d.h
	${CC} -c -I $(UTIL_DIR) $<

bwt_ari_decode.exe: bwt_mtf_rle_ari_d_main.o bwt_d.o arithmetic_decode.o
	${CC} $+ -o $@
