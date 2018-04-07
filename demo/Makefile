DIRS = structure excitation pirr aidr ionization polariz NeL
# Note: FeL needs a lot of calculation time, so we skip this for the daily
# testing.


all:
	@for i in $(DIRS); do cd $$i; make || exit 1; cd ../; done;

openmp:
	@for i in $(DIRS); do cd $$i; make openmp || exit 1; cd ../; done;

check:
	@for i in $(DIRS); do cd $$i; make check || exit 1; cd ../; done;

clean:
	@for i in $(DIRS); do cd $$i; make clean || exit 1; cd ../; done;

ref:
	@for i in $(DIRS); do cd $$i; make ref || exit 1; cd ../; done;


.PHONY: all check clean ref
