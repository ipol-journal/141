COBJ	  = mt19937ar.o io_png.o
CXXOBJ	  = libauxiliar.o libPDHG.o libProximal.o addnoise_ipol.o denoisingPDHG_ipol.o imdiff_ipol.o
BIN	  = addnoise_ipol denoisingPDHG_ipol imdiff_ipol


hdrdir    = -I/opt/local/include/ -I/usr/local/include/ -I/usr/include/
libdir    = -L/opt/local/lib/ -L/usr/local/lib/ -L/usr/lib/

COPT 	  = -O3 -funroll-loops -fomit-frame-pointer -ffast-math -ftree-vectorize -fopenmp
CXXFLAGS +=  $(COPT) -Wall -Wextra -Wno-write-strings -Wno-deprecated $(hdrdir)
CFLAGS   +=  $(COPT) $(hdrdir)
LDFLAGS  +=  $(CXXFLAGS) $(libdir) -lpng -lm -lgomp

default: $(COBJ) $(CXXOBJ)  $(BIN)

$(COBJ) : %.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@


$(OBJ) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@


$(BIN) : % : %.o mt19937ar.o io_png.o libauxiliar.o libPDHG.o libProximal.o
	$(CXX)   -o $@  $^ $(LDFLAGS)

.PHONY : clean
clean:
	$(RM) $(COBJ) $(CXXOBJ) ; rm $(BIN)
