.DEFAULT_GOAL := all

CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -D_POSIX_C_SOURCE=200809L
SRC = symnmf.c
OBJS = $(SRC:.c=.o)
EXEC = symnmf

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) -o $@ $(OBJS) $(CFLAGS) -lm 

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
