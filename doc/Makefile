
all: doc examples upload

RELEASE = alpha5
TARGET = farrowch@login.cacr.caltech.edu
DOCROOT = ~/docroot/diffraction/
PKGROOT = ~/dev_danse_us/

.PHONY : doc
doc:
	epydoc diffpy.srfit --html -vvv -o diffpy.srfitapi -n diffpy.srfit \
--include-log $@
	$(MAKE) -C sphinx $@

.PHONY : upload
upload:
	rsync -ruvz --delete diffpy.srfitapi $(TARGET):$(DOCROOT)
	ssh $(TARGET) "rm -rf $(DOCROOT)/diffpy.srfitapi-$(RELEASE)"
	ssh $(TARGET) "cp -r $(DOCROOT)/diffpy.srfitapi $(DOCROOT)/diffpy.srfitapi-$(RELEASE)"
	rsync -ruv srfit_examples.zip $(TARGET):$(PKGROOT)
	ssh $(TARGET) "rm -rf $(PKGROOT)/srfit_examples-$(RELEASE).zip"
	ssh $(TARGET) "cp -r $(PKGROOT)/srfit_examples.zip $(PKGROOT)/srfit_examples-$(RELEASE).zip"
	$(MAKE) -C sphinx $@

.PHONY : examples
examples:
	zip -r srfit_examples.zip ./examples -x \*.pyc -x \*svn\* -x \*threedoublepeaks\*

.PHONY : clean
clean:
	rm -rf diffpy.srfitapi
	rm -f srfit_examples.zip
	$(MAKE) -C sphinx $@

