! Monte Carlo simulation of light scattering by a spherical absorbing particle
module monte_carlo
  use ray_tracing
  use parallel
  implicit none

contains

  subroutine parallel_run(num_samples, radius, m1, m2, lambda, results)
    ! Run the Monte Carlo simulation in parallel for a given number of samples. 
    ! Inputs:
    !   num_samples - the number of samples to simulate
    !   radius - the radius of the spherical absorbing particle
    !   refractive_index - the complex refractive index of the environment and the sphere
    !   lambda - the wavelength of the incident light
    integer, intent(in) :: num_samples
    real(8), intent(in) :: radius
    complex(8), intent(in) :: m1, m2
    real(8), intent(in) :: lambda
    ! Outputs:
    !   results - the output array containing the distance traveled, total phase shift, intensity of each photon, and the direction
    !             of the photon after the simulation
    real(8), allocatable, intent(out) :: results(:, :)
    ! Local variables:
    !   origin - the starting point of each photon
    !   random1, random2 - random numbers for generating the starting points
    !   r, theta - polar coordinates for generating the starting points
    !   i - loop index
    real(8) :: random1, random2, r, theta
    real(8), allocatable :: origin(:,:)
    integer :: i
    ! Constants:
    real(8) :: PI = 3.14159265358979323846

    allocate(results(num_samples, 3))
    allocate(origin(num_samples, 4))
    ! We take the center of the sphere to be the origin of the coordinate system. All we need to do is shoot photons directly
    ! at the particle at the z-direction, hitting the -z-face. It could be argued that a smaller cross-section of the 
    ! sphere could be more accurate. 
    !$omp parallel do private(i, random1, random2, r, theta)
    do i = 1, num_samples
      call random_number(random1)
      call random_number(random2)
      r = radius * random1
      theta = 2 * PI * random2
      origin(i, 1) = r * cos(theta)
      origin(i, 2) = r * sin(theta)
      origin(i, 3) = -sqrt(radius ** 2 - r ** 2)
    end do
    !$omp end parallel do

    ! We run the simulation in parallel for each photon, as in this simulation each photon's journey is independent of the others. 
    ! (note: photon-photon scattering exists and can be added if needed, though this inclusion would mean adjustments in the 
    ! parallelization architecture)
    !$omp parallel do private(i, origin)
    do i = 1, num_samples
      call simulate_photon(origin(i, :), radius, m1, m2, lambda, results(i, :))
    end do
    !$omp end parallel do

  end subroutine parallel_run

  subroutine simulate_photon(origin, radius, m1, m2, lambda, result)
    ! Simulate the entire journey of a photon through the spherical absorbing particle. End the journey of 
    ! the photon when it leaves the particle or when its intensity falls below a certain threshold.
    ! Inputs: 
    !   origin - the starting point of the photon (3D vector)
    !   radius - the radius of the sphere
    !   m1, m2 - the complex refractive index of the environment and the sphere
    !   lambda - the wavelength of the photon
    ! Additional notes:
    !   The function uses the following subroutines from the ray_tracing module: propagate_photon, refract_photon
    real(8), intent(in) :: origin(3)
    real(8), intent(in) :: radius
    complex(8), intent(in) :: m1, m2
    real(8), intent(in) :: lambda
    ! Other parameters:
    !   direction - the direction of the photon's travel (3D vector)
    !   new_direction - the new direction of the photon after refraction (3D vector)
    !   new_origin - the new origin of the photon after propagation (3D vector)
    !   distance - the distance traveled by the photon in each step
    !   phase_shift - the phase shift experienced by the photon in each step
    !   intensity - the intensity of the photon
    !   medium - the medium the photon is currently in (-1 for air, -2 for absorbed by surface, 1 for the sphere)
    real(8) :: direction(3), new_direction(3), new_origin(3), distance, phase_shift, intensity = 1.0
    integer :: medium = -1
    ! Outputs:
    !   result - the output array containing the distance traveled, total phase shift, and intensity of the photon
    real(8), intent(out) :: result(4)
    ! Constants:
    real(8) :: PI = 3.14159265358979323846
    ! Initialize photon origin and direction
    direction = [0.0, 0.0, 1.0]  ! Direction of propagation

    ! First step: The photon propagates to the sphere and reaches its surface, some are reflected away from the surface
    call refract_photon(origin, direction, radius, m1, m2, lambda, medium, intensity, new_direction, phase_shift)

    ! Bouncing in the sphere: The photon moves in a linear fashion directly towards the next boundary. It can either reach 
    ! that boundary or be completely absorbed by the spherical volume. If it does reach the boundary, it can either be reflected 
    ! back with total internal reflection or be transmitted into the surrounding medium. We are interested in how much distance 
    ! each photon travels, how much phase shift it accumulates, and how much intensity it loses. These results will be stored in 
    ! the result array.   
    do while (intensity > 10 ** -3 .and. medium > 0)
      
      ! Protons move to second impact point, losing intensity due to volume absorption effects modeled by Beer-Lambert's law. 
      call propagate_photon(origin, direction, radius, distance)
      intensity = intensity * exp(-2 * PI * aimag(m2) * distance / lambda)
      ! Proton comes to contact with the sphere. Some of its intensity will be lost to surface interactions,
      ! modeled by Fresnel equations.
      call refract_photon(origin, direction, radius, m2, m1, lambda, medium, intensity, new_direction, phase_shift)
      direction = new_direction ! The photon continues in the new direction after the interaction
      result(1) = result(1) + distance
      result(2) = result(2) + phase_shift
    end do
    result(3) = intensity
    result(4) = direction
  
  end subroutine simulate_photon

end module monte_carlo
