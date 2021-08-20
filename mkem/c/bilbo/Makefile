LIB=libbilbo.a
HEADERS=api.h params.h common.h
OBJECTS=kem.o matrix_shake.o noise.o util.o mkem.o ../common/verify.o ../common/randombytes.o ../common/fips202.o

CFLAGS=-O3 -Wall -Wextra -Wpedantic -Wvla -Werror -Wmissing-prototypes -std=c99 -I../common $(EXTRAFLAGS) -DPREFIX=bilbo_ -fsanitize=address

all:  test_mkem

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(LIB): $(OBJECTS)
	$(AR) -r $@ $(OBJECTS)

test_mkem: $(LIB) test_mkem.c
	$(CC) $(CFLAGS) test_mkem.c $(LIB) -o test_mkem

test_speed: $(LIB) test_speed.c ../common/cpucycles.c ../common/speed_print.c
	$(CC) $(CFLAGS) test_speed.c $(LIB) -o test_speed ../common/speed_print.c ../common/cpucycles.c

clean:
	$(RM) $(OBJECTS)
	$(RM) $(LIB)
