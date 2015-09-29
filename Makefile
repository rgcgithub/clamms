all: normalize_coverage fit_models call_cnv gatk_coverage_to_bed

utils.o: utils.c
	gcc -c utils.c

hmm.o: hmm.c ltqnorm.c utils.h
	gcc -c hmm.c

normalize_coverage: normalize_coverage.o utils.o
	gcc normalize_coverage.o utils.o -o normalize_coverage -lm

normalize_coverage.o: normalize_coverage.c utils.h
	gcc -c normalize_coverage.c

fit_models: fit_models.o utils.o
	gcc fit_models.o utils.o -o fit_models -lm

fit_models.o: fit_models.c utils.h
	gcc -c fit_models.c

call_cnv: call_cnv.o utils.o hmm.o
	gcc call_cnv.o utils.o hmm.o -o call_cnv -lm

call_cnv.o: call_cnv.c utils.h hmm.h
	gcc -c call_cnv.c

gatk_coverage_to_bed: gatk_coverage_to_bed.o
	gcc gatk_coverage_to_bed.o -o gatk_coverage_to_bed

gatk_coverage_to_bed.o: gatk_coverage_to_bed.c
	gcc -c gatk_coverage_to_bed.c

clean:
	rm -rf *.o normalize_coverage fit_models call_cnv gatk_coverage_to_bed
