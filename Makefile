.PHONY: build doc test all install uninstall reinstall clean distclean

SETUP = ocaml setup.ml

build: setup.data
	$(SETUP) -build $(BUILDFLAGS)

all: build

doc: setup.data build
	$(SETUP) -doc $(DOCFLAGS)

test: setup.data build
	$(SETUP) -test $(TESTFLAGS)

uninstall: setup.data
	$(SETUP) -uninstall $(UNINSTALLFLAGS)

reinstall: setup.data
	ocamlfind remove biocaml
	$(SETUP) -install $(REINSTALLFLAGS)

install: reinstall

setup.ml: _oasis
	oasis setup -setup-update dynamic

setup.data: setup.ml
	$(SETUP) -configure $(CONFIGUREFLAGS)

configure: setup.data

clean:
	$(RM) -fr _build

distclean:
	$(RM) setup.data setup.log
	$(RM) configure
	$(RM) src/lib/META
	$(RM) src/lib/libbiocaml_stubs.clib
	$(RM) src/lib/doclib.odocl
	$(RM) src/lib/biocaml.mllib
	$(RM) TAGS
	$(SETUP) -distclean $(DISTCLEANFLAGS)

TAGS:
	otags -o TAGS `find src -regex ".*\.ml"`

CURR_DIR=$(shell basename $(CURDIR))
PKG=biocaml
VERSION=$(shell grep Version _oasis | cut -d' ' -f6)
PKG_VERSION=$(PKG)-$(VERSION)
DIST_FILES=Changes INSTALL LICENSE Makefile README.md TAGS _oasis _tags configure myocamlbuild.ml setup.ml doc src
.PHONY: dist
dist:
	oasis setup
	perl -pi -e 's#$(HOME)##g' myocamlbuild.ml setup.ml
	make doc
	mkdir doc
	mv _build/src/lib/doclib.docdir doc/html
	make TAGS
	cd .. ; mv $(CURR_DIR) $(PKG_VERSION); tar czf $(PKG_VERSION).tgz $(patsubst %,$(PKG_VERSION)/%,$(DIST_FILES)); mv $(PKG_VERSION) $(CURR_DIR)
	cd .. ; md5sum $(PKG_VERSION).tgz > $(PKG_VERSION).tgz.md5
	rm -rf doc
