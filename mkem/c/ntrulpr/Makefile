LIB=libntrulpr.a
HEADERS=ntrulpr.h
OBJECTS=mkem.o ../common/verify.o ../common/randombytes.o ../common/fips202.o  ntrulpr.o ntt_1907713_1536.o

CFLAGS=-O3 -Wall -Wextra -Wpedantic -Wvla -Werror -Wmissing-prototypes -std=c99 -I../common $(EXTRAFLAGS) -DPREFIX=ntrulpr_

all:  test_mkem test_ntrulpr

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(LIB): $(OBJECTS)
	$(AR) -r $@ $(OBJECTS)

test_ntrulpr: $(LIB) test_ntrulpr.c
	$(CC) $(CFLAGS) test_ntrulpr.c $(LIB) -o test_ntrulpr

test_mkem: $(LIB) test_mkem.c
	$(CC) $(CFLAGS) test_mkem.c $(LIB) -o test_mkem

test_speed: $(LIB) test_speed.c ../common/cpucycles.c ../common/speed_print.c
	$(CC) $(CFLAGS) test_speed.c $(LIB) -o test_speed ../common/speed_print.c ../common/cpucycles.c

clean:
	$(RM) $(OBJECTS)
	$(RM) $(LIB)
