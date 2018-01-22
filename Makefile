
objects = object/upvc_host.o    \
	object/upvc_dpu.o    \
	object/code.o        \
	object/getgenome.o   \
	object/getread.o     \
	object/index.o       \
	object/dispatch.o    \
	object/processread.o \
	object/compare.o     \
	object/vcf.o


version = upvc_1.1

upvc : $(objects) export.txt
	rm -f $(version)/object/*                 # dl
	cp upvc.h $(version)/                     # dl
	grep dl -v Makefile | sed s/'export.txt'// >$(version)/Makefile  # dl
	tar -cf $(version).tar $(version)         # dl
	gcc -O3 -Wall -lpthread -o upvc $(objects)

export.txt :                                      # dl
	mkdir -p $(version)                       # dl
	mkdir -p $(version)/object                # dl
	date > export.txt                         # dl

object/upvc_host.o : upvc.h upvc_host.c export.txt
	cp upvc_host.c $(version)/  # dl
	gcc -O3 -Wall -lpthread -c upvc_host.c -o object/upvc_host.o

object/upvc_dpu.o : upvc.h upvc_dpu.c export.txt
	cp upvc_dpu.c $(version)/ # dl
	gcc -O3 -Wall -c upvc_dpu.c -o object/upvc_dpu.o

object/code.o : upvc.h code.c export.txt
	cp code.c $(version)/ # dl
	gcc -O3 -Wall -c code.c -o object/code.o

object/compare.o : upvc.h compare.c export.txt
	cp compare.c $(version)/ # dl
	gcc -O3 -Wall -c compare.c -o object/compare.o

object/getgenome.o : upvc.h getgenome.c export.txt
	cp getgenome.c $(version)/ # dl
	gcc -O3 -Wall -c getgenome.c -o object/getgenome.o

object/getread.o : upvc.h getread.c export.txt
	cp getread.c $(version)/ # dl
	gcc -O3 -Wall -c getread.c -o object/getread.o

object/index.o : upvc.h index.c export.txt
	cp index.c $(version)/ # dl
	gcc -O3 -Wall -c index.c -o object/index.o

object/dispatch.o : upvc.h dispatch.c export.txt
	cp dispatch.c $(version)/ # dl
	gcc -O3 -Wall -c dispatch.c -o object/dispatch.o

object/processread.o : upvc.h processread.c export.txt
	cp processread.c $(version)/ # dl
	gcc -O3 -Wall -c processread.c -o object/processread.o

object/vcf.o : upvc.h vcf.c export.txt
	cp vcf.c $(version)/ # dl
	gcc -O3 -Wall -c vcf.c -o object/vcf.o

