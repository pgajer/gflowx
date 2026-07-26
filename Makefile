.PHONY: clean attrs document test build check-fast check

VERSION := $(shell grep "^Version:" DESCRIPTION | sed 's/Version: //')
PKGNAME := gflowx
TARBALL := $(PKGNAME)_$(VERSION).tar.gz

clean:
	find src -name "*.o" -delete
	find src -name "*.so" -delete
	rm -f src/*.dll
	rm -rf $(PKGNAME).Rcheck
	rm -f $(TARBALL)

attrs:
	R -q -e "Rcpp::compileAttributes()"

document: attrs
	R -q -e "roxygen2::roxygenise(load = 'source')"

test:
	Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_dir("tests/testthat")'

build: clean document
	R CMD build .

check-fast: build
	R CMD check $(TARBALL) --as-cran --no-examples --no-tests --no-manual

check: build
	R CMD check $(TARBALL) --as-cran
