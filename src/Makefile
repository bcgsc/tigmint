all: long-to-linked-pe

CXXFLAGS+=-O3 -Wall -Werror -Wno-unknown-pragmas -std=c++11 -pthread
LDFLAGS+=-l btllib

long-to-linked-pe: long-to-linked-pe.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f long-to-linked-pe

# Check the C++ source code for errors.
lint: clang-format clang-tidy

# Check the C++ source code for errors with clang-tidy.
clang-tidy:
	clang-tidy -warnings-as-errors='*' *.cpp -- -std=c++11 -x c++ -Ibtllib/include -Iinclude $(CXXFLAGS) 

# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in *.cpp; do clang-format -style=file $$i >$$i.fixed; done
	for i in *.cpp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi
