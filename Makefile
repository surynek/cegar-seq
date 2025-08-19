SUBDIRS = libslic3r src test

all: optimized

optimized:
	for dir in $(SUBDIRS); do make -C $$dir; done

clean:
	for dir in $(SUBDIRS); do make -C $$dir clean; done
	rm -f *~ *.o *.a include/*~
