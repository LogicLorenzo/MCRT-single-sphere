# Monte Carlo Ray Tracing

This project implements a Monte Carlo ray tracing algorithm for light hitting a single sphere in Fortran, utilizing parallel programming techniques to enhance performance. The goal is to estimate the average phase-shift, and attenuation of the light. 

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

## Building the Project

To compile the project, run the following command in the terminal:

```
make
```

This will generate the executable for the Monte Carlo ray tracing application.

## Usage

After building the project, you can run the application with:

```
./monte_carlo_ray_tracing
```

You may need to specify input parameters depending on the implementation. Currently all parameters are simply coded in the software, though the change is simple. 

## Algorithms

- **Monte Carlo Simulation**: A statistical method used to estimate the expected value of a function by sampling random points.
- **Ray Tracing**: A rendering technique for generating images by tracing the path of rays through a scene.

## Contributing

Feel free to contribute to this project by submitting issues or pull requests. Your feedback and improvements are welcome!
