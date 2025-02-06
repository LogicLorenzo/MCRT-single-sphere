program main
  use monte_carlo
  use utils
  implicit none

  integer :: num_samples, radius
  real(8), allocatable :: results(:,:)
  complex(8) = m1, m2

  ! Initialize parameters
  num_samples = 1000
  radius = 20.0 !nm
  m1 = (1.0, 0.0)  ! Complex refractive index of environment
  m2 = (1.5, 1.0)  ! Complex refractive index of sphere

  ! Allocate results array
  allocate(results(num_samples, 2))  ! Columns for distance and angle

  ! Run Monte Carlo simulation
  call parallel_run(num_samples, radius, refractive_index, results)

  ! Save the results to a CSV file
  call save_results_to_csv('data/results.csv', results)

  ! Deallocate results array
  deallocate(results)

end program main