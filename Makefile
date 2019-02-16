ifeq ($(PADDING),)
PADDING_ARG :=
else
PADDING_ARG := --padding=$(PADDING)
endif

amr_library: build_gene_files build_indices
	$(MAKE) optimize_amr_library
	$(MAKE) extract_guides ARGS="--library generated_files/library.txt \
				     --genes inputs/amr/full_guide_set/amr_staph_genes.txt \
				     --max-cuts 10 \
				     --output generated_files/amr_library.txt"

optimize_amr_library:
	$(MAKE) optimizer ARGS="--include inputs/amr/full_guide_set/all_ordered_guides.txt --output generated_files/library.txt"

extract_guides:
	python extract_guides.py $(ARGS)

build_gene_files:
	python make_genes_and_identify_all_targets.py

crispr_sites2: crispr_sites2.cpp
	g++ -O3 --std=c++11 -o crispr_sites2 crispr_sites2.cpp

generated_files/all_offtargets.txt: inputs/additional/offtargets/hg38.fa.gz crispr_sites2
	gzip -dc inputs/additional/offtargets/hg38.fa.gz | ./crispr_sites2 > generated_files/human_guides_38.txt
	gzip -dc inputs/additional/offtargets/EColi_BL21_DE3__NC_012892.2.fasta.gz | ./crispr_sites2 > generated_files/ecoli_bl21_de3_offtargets.txt
	cat generated_files/ecoli_bl21_de3_offtargets.txt generated_files/human_guides_38.txt > generated_files/all_offtargets.txt

generate_offtargets: generated_files/all_offtargets.txt

compile_offtarget_server: offtarget/matcher.go offtarget/main.go
	cd offtarget && go build

inputs/additional/offtargets/hg38.fa.gz:
	curl http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz > inputs/additional/offtargets/hg38.fa.gz.download && \
	mv inputs/additional/offtargets/hg38.fa.gz.download inputs/additional/offtargets/hg38.fa.gz && \
	touch $@

build_indices: generate_offtargets
	$(MAKE) start_offtarget_server
	python build_indices.py
	$(MAKE) stop_offtarget_server

start_offtarget_server: compile_offtarget_server
	cd offtarget && { HOST=file://`pwd`/../generated_files/all_offtargets.txt ./offtarget & echo $$! > server.PID; }

stop_offtarget_server:
	kill `cat offtarget/server.PID` && rm offtarget/server.PID

optimizer:
	python optimizer.py $(ARGS)

test:
	cd tests && pytest
	cd tests && ./regression_test.sh inputs/fasta_files/simple_guide_test.fasta \
											             outputs/expected_simple_guide_test.txt

# make library TARGETS=inputs/additional/colistin.fasta OUTPUT=library.txt
library:
	python make_genes_and_identify_all_targets.py --targets=$(TARGETS) $(PADDING_ARG) --disable-git
	$(MAKE) build_indices
	$(MAKE) optimizer ARGS="--output $(OUTPUT) $(PADDING_ARG)"

# make library_excluding TARGETS=inputs/additional/colistin.fasta OUTPUT=library.txt EXCLUDE=exclude.txt
library_excluding:
	python make_genes_and_identify_all_targets.py --targets=$(TARGETS) $(PADDING_ARG) --disable-git
	$(MAKE) build_indices
	$(MAKE) optimizer ARGS="--output $(OUTPUT) --exclude $(EXCLUDE) $(PADDING_ARG)"

library_including:
	python make_genes_and_identify_all_targets.py --targets=$(TARGETS) $(PADDING_ARG) --disable-git
	$(MAKE) build_indices
	$(MAKE) optimizer ARGS="--output $(OUTPUT) --include $(INCLUDE) $(PADDING_ARG)"
