## Copyright (C) 2015-2016 David Miguel Susano Pinto
##
## Copying and distribution of this file, with or without modification,
## are permitted in any medium without royalty provided the copyright
## notice and this notice are preserved.  This file is offered as-is,
## without any warranty.

AWK ?= awk
OCTAVE ?= octave
IMAGEJ ?= imagej

DATA := \
  Airyscan-data

ALL_VALIDATE := $(addprefix validate-, $(DATA))

.PHONY: clean help validate $(ALL_VALIDATE) centroids

IMAGES := $(shell awk '{print $$2}' data/Airyscan-data.md5sum)

help:
	@echo "Targets:"
	@echo "   centroids         - CSV files with centroid coordinates"
	@echo
	@echo "   clean             - Remove results"
	@echo
	@echo "Environment variables:"
	@echo "   AWK     - AWK program (default: awk)"
	@echo "   IMAGEJ  - ImageJ program (default: imagej)"
	@echo "   OCTAVE  - octave interpreter program (default: octave)"

$(ALL_VALIDATE):
	md5sum -c $(patsubst validate-%,data/%.md5sum, $@)

validate: $(ALL_VALIDATE)

centroids_coordinates_file = results/$(1)-centroids.tsv
centroids_log_file = results/$(1)-log.tif
centroids_results = \
  $(call centroids_coordinates_file,$(1)) \
  $(call centroids_log_file,$(1))

ALL_CENTROIDS := \
  $(foreach fname, $(patsubst data/%.czi,%, $(IMAGES)), \
      $(call centroids_results,$(fname)))

$(call centroids_results,%): data/%.czi src/weighted-centroids.m | validate results
	$(OCTAVE) src/weighted-centroids.m $< \
	  $(call centroids_log_file,$*) > $(call centroids_coordinates_file,$*)

centroids: $(ALL_CENTROIDS)


results:
	mkdir results

clean:
	$(RM) $(ALL_CENTROIDS)
