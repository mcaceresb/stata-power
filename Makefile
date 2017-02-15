CC = gcc
SPI = 2.0
SPT = 0.1
CFLAGS = -shared -fopenmp -fPIC -DSYSTEM=OPUNIX
GSLFLAGS = -lgsl -lgslcblas -lm -static
SPT_C = lib/spt-$(SPT)/stutils.c
SPT_H = lib/spt-$(SPT)/stutils.h
ST_C = lib/spi-$(SPI)/stplugin.c
ST_H = lib/spi-$(SPI)/stplugin.h
OUT = lib/spi-$(SPI)/stplugin.o src/psimci.o src/stplugin.c src/stplugin.h
DEPS = src/psimci.

all: clean linksptc linkspth linkstc linksth $(OUT) pluginmake

linksptc: $(SPT_C)
	ln -srf $(SPT_C) src/stutils.c

linkspth: $(SPT_H)
	ln -srf $(SPT_H) src/stutils.h

linkstc: $(ST_C)
	ln -srf $(ST_C) src/stplugin.c

linksth: $(ST_H)
	ln -srf $(ST_H) src/stplugin.h

%.o: %.c $(DEPS)
	$(CC) -Wall -c -o $@ $< $(CFLAGS)

pluginmake: lib/spi-$(SPI)/stplugin.o src/psimci.o
	gcc $(CFLAGS) lib/spi-$(SPI)/stplugin.o src/psimci.o $(GSLFLAGS) -o psimci.plugin

.PHONY: clean
clean:
	rm -f $(OUT)
