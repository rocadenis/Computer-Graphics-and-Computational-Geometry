NAME_STUDENT_1=Roca_Denis_E1
NAME_STUDENT_2=Fechita_Diandra_E3
NAME_HOMEWORK=cg1

#fill these in ^
#You can then `make export` to have the appropriately-named file for sending.

CC=gcc
SRC=$(NAME_HOMEWORK).cpp
BINARY=$(NAME_HOMEWORK).bin
CPPFLAGS=-O2 -g0
LDFLAGS=-lGL -lGLU -lglut

EXPORT=$(NAME_STUDENT_1)_$(NAME_STUDENT_2)_$(NAME_HOMEWORK).cpp

#Don't actually send archives, send single source files.
#ARCHIVE=$(NAME_STUDENT_1)_$(NAME_STUDENT_2)_$(NAME_HOMEWORK).zip

all: run

$(BINARY): $(SRC)
	g++ -o $@ $^ $(CPPFLAGS) $(LDFLAGS)

run: $(BINARY)
	./$<

clean:
	rm $(BINARY) *~ $(ARCHIVE) 2>/dev/null || echo -n ''

export: clean
	cp -f $(SRC) $(EXPORT)

#Again: don't send archives.
#archive: clean
#	zip $(ARCHIVE) $(SRC)
