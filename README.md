# Monte Carlo Ray Tracing for a Single Sphere

This project implements a Monte Carlo ray tracing algorithm for light hitting a single sphere in Fortran, utilizing parallel programming techniques to enhance performance. The goal is to estimate the average phase-shift and attenuation of the light. 

## Project Structure

```
monte-carlo-ray-tracing
├── src
│   ├── main.f90          # Entry point of the application
│   ├── monte_carlo.f90   # Monte Carlo simulation implementation
│   ├── ray_tracing.f90   # Ray tracing logic
│   └── utils.f90         # Utility functions and data structures
├── Makefile               # Build instructions
└── README.md              # Project documentation
```

## Setup Instructions

1. Ensure you have a Fortran compiler installed (e.g., gfortran).
2. Clone the repository or download the project files.
3. Navigate to the project directory.

## Running the Program

To compile the project, run the following command in the terminal:

```
make
```

This will generate the executable for the Monte Carlo ray tracing application.


After building the project, you can run the application with:

```
./monte_carlo_ray_tracing
```

You may need to specify input parameters depending on the implementation. Currently, all parameters are given in the software, though the change is simple. 

## Parallel Computing 

The project utilizes OpenMP for parallel computing to speed up the Monte Carlo simulations. Ensure that your Fortran compiler supports OpenMP.

## Main Program Logic

The main program performs the following steps:

1. Shoots photons at the -z face of the sphere. 
2. Photon is transmitted into the sphere. 
4. Photon bounces around inside the sphere until it loses most of its energy. 
5. Outputs the total phase shift, attenuation and direction of all photons. 

## Future Works

- Add Brownian motion for photons in the case of a highly scattering material

## Contributions

Feel free to contribute to this project by submitting issues or pull requests. Your feedback and improvements are welcome!
