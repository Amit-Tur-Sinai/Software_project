.DEFAULT_GOAL := all

CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
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
