SRCS = $(shell find ./src/ -name *.c)
OBJS = $(SRCS:.c=.o)
HEADERS = $(shell find ./inc/ -name *.h)

all : Riemann_solver.a

Riemann_solver.a : $(SRCS) $(HEADERS)
	$(MAKE) --directory=./src/
	ar crv Riemann_solver.a $(OBJS)
	ranlib Riemann_solver.a


.PHONY : clean all test



clean :
	$(MAKE) --directory=./src/ clean
	rm -f *.[oa]

test :
	@echo =========================
	@echo $(SRCS)
