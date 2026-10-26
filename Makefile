PYTHON ?= python3
CC ?= cc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?= -lz

.PHONY: all native-check test demo clean
all: scripts/bam_merge scripts/bam_dechimer

scripts/bam_merge: scripts/bam_merge.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

scripts/bam_dechimer: scripts/bam_dechimer.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

native-check:
	cd scripts && $(PYTHON) setup.py build_ext --inplace

test: all
	$(PYTHON) -c 'import shutil; assert all(shutil.which(t) for t in ("samtools", "pigz")), "samtools and pigz are required"'
	$(PYTHON) -m unittest discover -s tests -v

demo: all
	$(PYTHON) examples/roundtrip.py

clean:
	rm -f scripts/bam_merge scripts/bam_dechimer scripts/fastcheck*.so
	rm -rf scripts/build scripts/__pycache__ tests/__pycache__ examples/__pycache__
