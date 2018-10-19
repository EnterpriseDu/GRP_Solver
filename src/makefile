
ARXIV = $(notdir $(CURDIR)).a
SRCS = $(wildcard *.c)
OBJS = $(SRCS:.c=.o)
CC = gcc
CFLAGS = -I../inc -c -g


all : $(OBJS)

%.o : %.c
	$(CC) $^ $(CFLAGS)










clean :
	rm -f $(OBJS) $(OBJS:.o=.d)

.PHONY : clean all test




test :
	@echo =========================
	@echo $(OBJS)
	@echo $(SRCS)
