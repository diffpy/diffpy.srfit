
all: doc upload

.PHONY : doc
doc:
	sphinx-build . diffpy.srfit

.PHONY : upload
upload:
	rsync -ruvz --delete diffpy.srfit farrowch@login.cacr.caltech.edu:~/docroot/diffraction/

.PHONY : clean
clean:
	rm -rf diffpy.srfit

