FC := gfortran
FLAGS := -fopenmp
OBJS := main.o monte_carlo.o ray_tracing.o parallel.o utils.o

all: MCRT-sphere

MCRT-sphere: $(OBJS)
    $(FC) $(FLAGS) -o $@ $(OBJS)

%.o: %.f90
    $(FC) $(FLAGS) -c $<

clean:
    rm -f *.o *.mod MCRT-sphere