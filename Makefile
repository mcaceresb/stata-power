CC = gcc
SPI = 2.0
DEPS = lib/spi-$(SPI)/stplugin.c src/psimci.c
OUT = lib/spi-2.0/stplugin.o lib/spi-3.0/stplugin.o src/psimci.o
CFLAGS = -fopenmp -shared -fPIC -DSYSTEM=OPUNIX
GSLFLAGS = -lgsl -lgslcblas -lm

link: lib/spi-$(SPI)/stplugin.h
	ln -srf lib/spi-$(SPI)/stplugin.h src/stplugin.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

pluginmake: lib/spi-$(SPI)/stplugin.o src/psimci.o
	gcc $(CFLAGS) lib/spi-$(SPI)/stplugin.o src/psimci.o $(GSLFLAGS) -o psimci.plugin

all: clean $(OUT)
.PHONY: clean
clean:
	rm -f $(OUT)
