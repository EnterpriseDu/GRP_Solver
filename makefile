SRCS = $(shell find ./src/ -name *.c)
OBJS = $(SRCS:.c=.o)
HEADERS = $(shell find ./inc/ -name *.h)

all : GRP_Solver.a

GRP_Solver.a : $(SRCS) $(HEADERS)
	$(MAKE) --directory=./src/
	ar crv GRP_Solver.a $(OBJS)
	ranlib GRP_Solver.a


.PHONY : clean all test



clean :
	$(MAKE) --directory=./src/ clean
	rm -f *.[oa]

test :
	@echo =========================
	@echo $(SRCS)
