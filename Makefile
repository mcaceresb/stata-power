CC = gcc
SPI = 2.0
ST_C = lib/spi-$(SPI)/stplugin.c
ST_H = lib/spi-$(SPI)/stplugin.h
DEPS = src/psimci.c
OUT = lib/spi-2.0/stplugin.o lib/spi-3.0/stplugin.o src/psimci.o src/stplugin.c src/stplugin.h
CFLAGS = -fopenmp -shared -fPIC -DSYSTEM=OPUNIX
GSLFLAGS = -lgsl -lgslcblas -lm

all: clean linkc linkh $(OUT)

linkc: $(ST_C)
	ln -srf $(ST_C) src/stplugin.c

linkh: $(ST_H)
	ln -srf $(ST_H) src/stplugin.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

pluginmake: lib/spi-$(SPI)/stplugin.o src/psimci.o
	gcc $(CFLAGS) lib/spi-$(SPI)/stplugin.o src/psimci.o $(GSLFLAGS) -o psimci.plugin

.PHONY: clean
clean:
	rm -f $(OUT)
