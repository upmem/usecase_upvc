objects = object/upvc_host.o    \
	object/upvc_dpu.o    \
	object/code.o        \
	object/getgenome.o   \
	object/getread.o     \
	object/index.o       \
	object/dispatch.o    \
	object/processread.o \
	object/compare.o     \
	object/normalizevar.o \
	object/vartree.o \
	object/vcf.o

.PHONY: clean $(objects)

# default optimization index is 3
OPTIMIZATION_INDEX = 3

all : upvc
debug : PREFLAGS += -g
debug : OPTIMIZATION_INDEX = 0
debug : all

#PREFLAGS = -O$(OPTIMIZATION_INDEX) -Wall -Wno-format-overflow
PREFLAGS = -O$(OPTIMIZATION_INDEX) -Wall


version = upvc_1.4

upvc : $(objects) 
	gcc $(PREFLAGS) -o upvc $(objects) -lm
#	gcc $(PREFLAGS) -pthread -o upvc $(objects) -lm


object/upvc_host.o : upvc.h upvc_host.c 
	gcc $(PREFLAGS) -pthread -c upvc_host.c -o object/upvc_host.o

object/upvc_dpu.o : upvc.h upvc_dpu.c 
	gcc $(PREFLAGS) -c upvc_dpu.c -o object/upvc_dpu.o

object/code.o : upvc.h code.c 
	gcc $(PREFLAGS) -c code.c -o object/code.o

object/compare.o : upvc.h compare.c 
	gcc $(PREFLAGS) -c compare.c -o object/compare.o

object/getgenome.o : upvc.h getgenome.c 
	gcc $(PREFLAGS) -c getgenome.c -o object/getgenome.o

object/getread.o : upvc.h getread.c 
	gcc $(PREFLAGS) -c getread.c -o object/getread.o

object/index.o : upvc.h index.c 
	gcc $(PREFLAGS) -c index.c -o object/index.o

object/dispatch.o : upvc.h dispatch.c 
	gcc $(PREFLAGS) -c dispatch.c -o object/dispatch.o

object/processread.o : upvc.h processread.c 
	gcc $(PREFLAGS) -c processread.c -o object/processread.o

object/vartree.o : upvc.h vartree.c 
	gcc $(PREFLAGS) -c vartree.c -o object/vartree.o

object/normalizevar.o :  upvc.h normalizevar.c 
	gcc $(PREFLAGS) -c normalizevar.c -o object/normalizevar.o

object/vcf.o : upvc.h vcf.c 
	gcc $(PREFLAGS) -c vcf.c -o object/vcf.o

clean:
	@(rm -f upvc $(objects))
