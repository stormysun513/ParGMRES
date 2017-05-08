SUBDIRS = dense_seq dense_omp sparse_omp
SUBDIRSCLEAN=$(addsuffix .clean, $(SUBDIRS))

.PHONY: $(SUBDIRS) clean

default: $(SUBDIRS)

clean: $(SUBDIRSCLEAN)

%.clean:
	@make clean -C $(basename $@)

$(SUBDIRS):
	@make -C $@
