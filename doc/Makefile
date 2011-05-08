
all: doc examples upload

RELEASE = alpha9
TARGET = farrowch@login.cacr.caltech.edu
DOCROOT = ~/docroot/diffraction/
PKGROOT = ~/dev_danse_us/

.PHONY : doc
doc:
	epydoc diffpy.srfit --html -vvv -o diffpy.srfitapi -n diffpy.srfit \
--include-log --exclude diffpy.srfit.structure.cctbxstructure $@
	$(MAKE) -C devmanual $@

.PHONY : upload
upload:
	rsync -ruvz --delete diffpy.srfitapi $(TARGET):$(DOCROOT)
	ssh $(TARGET) "rm -rf $(DOCROOT)/diffpy.srfitapi-$(RELEASE)"
	ssh $(TARGET) "cp -r $(DOCROOT)/diffpy.srfitapi $(DOCROOT)/diffpy.srfitapi-$(RELEASE)"
	rsync -ruv srfit_examples.zip $(TARGET):$(PKGROOT)
	ssh $(TARGET) "rm -rf $(PKGROOT)/srfit_examples-$(RELEASE).zip"
	ssh $(TARGET) "cp -r $(PKGROOT)/srfit_examples.zip $(PKGROOT)/srfit_examples-$(RELEASE).zip"
	$(MAKE) -C devmanual $@

.PHONY : examples
examples:
	zip -r srfit_examples.zip ./examples/*.py ./examples/data -x \*svn\* -x \*threedoublepeaks\* -x \*temp\* -x \*test\*

.PHONY : clean
clean:
	rm -rf diffpy.srfitapi
	rm -f srfit_examples.zip
	$(MAKE) -C devmanual $@

