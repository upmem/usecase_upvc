
objects = object/upvc_host.o \
	object/upvc_dpu.o    \
	object/code.o        \
	object/getgenome.o   \
	object/getread.o     \
	object/index.o       \
	object/dispatch.o    \
	object/processread.o \
	object/compare.o     \
	object/vcf.o

version = upvc_1.0

upvc : $(objects) 
	gcc -O3 -Wall -lpthread -o upvc $(objects)


object/upvc_host.o : upvc.h upvc_host.c 
	gcc -O3 -Wall -lpthread -c upvc_host.c -o object/upvc_host.o

object/upvc_dpu.o : upvc.h upvc_dpu.c 
	gcc -O3 -Wall -c upvc_dpu.c -o object/upvc_dpu.o

object/code.o : upvc.h code.c 
	gcc -O3 -Wall -c code.c -o object/code.o

object/compare.o : upvc.h compare.c 
	gcc -O3 -Wall -c compare.c -o object/compare.o

object/getgenome.o : upvc.h getgenome.c 
	gcc -O3 -Wall -c getgenome.c -o object/getgenome.o

object/getread.o : upvc.h getread.c 
	gcc -O3 -Wall -c getread.c -o object/getread.o

object/index.o : upvc.h index.c 
	gcc -O3 -Wall -c index.c -o object/index.o

object/dispatch.o : upvc.h dispatch.c 
	gcc -O3 -Wall -c dispatch.c -o object/dispatch.o

object/processread.o : upvc.h processread.c 
	gcc -O3 -Wall -c processread.c -o object/processread.o

object/vcf.o : upvc.h vcf.c 
	gcc -O3 -Wall -c vcf.c -o object/vcf.o

