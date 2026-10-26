# Makefile — spring network virial pressure
#
# Step 1:  make run S=0        build and run MD (one file per T)
# Step 2:  jupyter notebook pressure_analysis.ipynb

CC     = gcc
CFLAGS = -O2 -std=c99
LIBS   = -lm
S     ?= 0

pressure_vs_T: spring_networks_finite_T_dynamics.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

run: pressure_vs_T
	./pressure_vs_T $(S)

clean:
	rm -f pressure_vs_T

.PHONY: run clean
