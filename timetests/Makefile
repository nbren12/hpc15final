
CXX = c++
CXXFLAGS = -llaplace -larmadillo -O3 -L../src/ -I../src/

all: test_run


test_sep_time.exe: test_sep_time.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $^

test_run: test_sep_time.exe
	./test_sep_time.exe

clean:
	rm -rf $(EXECUTABLES)
	rm -vf *~
