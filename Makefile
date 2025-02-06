FC := gfortran
FLAGS := -fopenmp
OBJS := main.o monte_carlo.o ray_tracing.o parallel.o utils.o

all: monte_carlo_ray_tracing

monte_carlo_ray_tracing: $(OBJS)
    $(FC) $(FLAGS) -o $@ $(OBJS)

%.o: %.f90
    $(FC) $(FLAGS) -c $<

clean:
    rm -f *.o *.mod monte_carlo_ray_tracing