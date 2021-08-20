CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -O3 -fomit-frame-pointer -I../common
RM = /bin/rm

SOURCES = kex.c kem.c indcpa.c polyvec.c poly.c ntt.c cbd.c reduce.c ../common/verify.c mkem.c
SOURCESKECCAK = $(SOURCES) ../common/fips202.c symmetric-shake.c
HEADERS = params.h kex.h kem.h indcpa.h polyvec.h poly.h ntt.h cbd.h reduce.c ../common/verify.h symmetric.h
HEADERSKECCAK = $(HEADERS) ../common/fips202.h

.PHONY: all speed shared clean

all: \
  test_kyber512 \
  test_mkem512 \
  test_kex512 \
  test_vectors512

speed: \
  test_speed512

shared: \
  libpqcrystals_kyber512_ref.so \
  libpqcrystals_fips202_ref.so

libpqcrystals_fips202_ref.so: ../common/fips202.c ../common/fips202.h
	$(CC) -shared -fPIC $(CFLAGS) ../common/fips202.c -o libpqcrystals_fips202_ref.so


libpqcrystals_kyber512_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	$(CC) -shared -fPIC $(CFLAGS) -DKYBER_K=2 $(SOURCES) symmetric-shake.c -o libpqcrystals_kyber512_ref.so

test_mkem512: $(SOURCESKECCAK) $(HEADERSKECCAK) test_mkem.c ../common/randombytes.c
	$(CC) $(CFLAGS) -DKYBER_K=2 $(SOURCESKECCAK) ../common/randombytes.c test_mkem.c -o test_mkem512

test_kyber512: $(SOURCESKECCAK) $(HEADERSKECCAK) test_kyber.c ../common/randombytes.c
	$(CC) $(CFLAGS) -DKYBER_K=2 $(SOURCESKECCAK) ../common/randombytes.c test_kyber.c -o test_kyber512

test_kex512: $(SOURCESKECCAK) $(HEADERSKECCAK) test_kex.c ../common/randombytes.c
	$(CC) $(CFLAGS) -DKYBER_K=2 $(SOURCESKECCAK) ../common/randombytes.c test_kex.c -o test_kex512

test_vectors512: $(SOURCESKECCAK) $(HEADERSKECCAK) test_vectors.c
	$(CC) $(CFLAGS) -DKYBER_K=2 $(SOURCESKECCAK) test_vectors.c -o test_vectors512

test_speed512: $(SOURCESKECCAK) $(HEADERSKECCAK) ../common/cpucycles.h ../common/cpucycles.c ../common/speed_print.h ../common/speed_print.c test_speed.c ../common/randombytes.c
	$(CC) $(CFLAGS) -DKYBER_K=2 $(SOURCESKECCAK) ../common/randombytes.c ../common/cpucycles.c ../common/speed_print.c test_speed.c -o test_speed512


clean:
	-$(RM) -rf *.gcno *.gcda *.lcov *.o *.so
	-$(RM) -rf test_kyber512
	-$(RM) -rf test_mkem512
	-$(RM) -rf test_kex512
	-$(RM) -rf test_vectors512
	-$(RM) -rf test_speed512
	-$(RM) -rf PQCgenKAT_kem512
