#For linux OS with ifort
#FLAGS
FC= ifort
FFLAGS = -r8 -O2
#sources files
SRC= assemble.f conv.f damt.f fem.f flux.f gauss_solver.f \
     init.f post_out.f main.f

OBJS = $(SRC:.f=.o)
LIBS =

Hbfem:$(OBJS)
    $(FC) -o execu $(FFLAGS) $(OBJS) $(LIBS)
clean:
   rm -f *.o 
