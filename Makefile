CXXFLAGS := ${CXXFLAGS} -W -Wall -Wextra -pedantic -std=c++2a -O2
LDFLAGS := -pthread

PROG   := bin/crispy
OBJS   := build/count.o \
          build/crispy++.o \
          build/score.o \
          build/tools.o
DFILES := $(OBJS:.o=.d)


.PHONY: all clean

all: ${PROG}

clean:
	rm -vf ${PROG} ${OBJS} ${DFILES}
	rmdir bin
	rmdir build


${PROG}: ${OBJS}
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) ${LDFLAGS} $^ -o $@


# Object files
build/%.o: src/%.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -pthread -c -o $@ $<
	$(CXX) $(CXXFLAGS) -w -MM -MT $@ -MF $(@:.o=.d) $<


# Automatic header depencencies
-include $(DFILES)
