
all: doc upload

RELEASE = beta0
TARGET = farrowch@login.cacr.caltech.edu
DOCROOT = ~/docroot/diffraction/
PKGROOT = ~/dev_danse_us/

.PHONY : doc
doc:
	sphinx-build . diffpy.srfit

.PHONY : upload
upload:
	rsync -ruvz --delete diffpy.srfit $(TARGET):$(DOCROOT)
	ssh $(TARGET) "rm -rf $(DOCROOT)/diffpy.srfit-$(RELEASE)"
	ssh $(TARGET) "cp -r $(DOCROOT)/diffpy.srfit $(DOCROOT)/diffpy.srfit-$(RELEASE)"

.PHONY : clean
clean:
	rm -rf diffpy.srfit

