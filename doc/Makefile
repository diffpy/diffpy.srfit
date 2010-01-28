
all: doc example upload

.PHONY : doc
doc:
	epydoc diffpy.srfit --html -vvv -o diffpy.srfitapi -n diffpy.srfit \
--include-log $@
	$(MAKE) -C sphinx $@

.PHONY : upload
upload:
	rsync -ruvz --delete diffpy.srfitapi farrowch@login.cacr.caltech.edu:~/docroot/diffraction/
	rsync -ruv srfit_examples.zip farrowch@login.cacr.caltech.edu:~/dev_danse_us/
	$(MAKE) -C sphinx $@

.PHONY : examples
examples:
	zip -r srfit_examples.zip ./examples -x \*.pyc -x \*svn\* -x \*threedoublepeaks\*

.PHONY : clean
clean:
	rm -rf diffpy.srfitapi
	rm -f srfit_examples.zip
	$(MAKE) -C sphinx $@

