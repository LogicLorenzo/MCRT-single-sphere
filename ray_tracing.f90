module ray_tracing
  implicit none

contains

  subroutine propagate_photon(origin, direction, radius, distance)
    ! The photon starts from one point on the sphere and starts moving away from the surface. 
    ! In/out parameters:
    !   origin - the impact point of the photon
    real(8), intent(inout) :: origin(3)
    ! Inputs:
    !   direction - the direction of the photon
    !   radius - the radius of the sphere
    real(8), intent(in) :: direction(3)
    real(8), intent(in) :: radius
    ! Outputs:
    !   distance - the distance the photon travels before hitting the sphere
    real(8), intent(out) :: distance
    ! Local variables:
    real(8) :: a, b  ! Coefficients for the quadratic equation
    
    ! Geometrically, the photon travels in a straight line. We can calculate
    ! the intersection point of the photon with the sphere by solving a 
    ! a simple equation.
    ! (x + da) ** 2 + (y + db) ** 2 + (z + dc) ** 2 = r ** 2
    ! (a^2 + b^2 + c^2) d^2 + 2(a*x + b*y + c*z) d + (x^2 + y^2 + z^2 - r^2) = 0
    ! d = -2(a*x + b*y + c*z) / (a^2 + b^2 + c^2)
    a = dot_product(direction, direction) 
    b = 2 * dot_product(origin, direction)
    distance = -b / a
    origin = origin + distance * direction

  end subroutine propagate_photon

  subroutine refract_photon(origin, direction, radius, m1, m2, lambda, medium, intensity, new_direction, phase_shift)
    ! The photon hits the sphere and is transmitted, reflected, and absorbed based on the Fresnel equations.
    ! In/out parameters:
    !   medium - the medium the photon is currently in (-1 for air, 1 for the sphere)
    !   intensity - the intensity of the photon
    integer, intent(inout) :: medium
    real(8), intent(inout) :: intensity
    ! Inputs:
    !   direction - the direction of the photon (normalized)
    !   origin - the impact point of the photon
    !   radius - the radius of the sphere
    !   m1, m2 - the complex refractive index of medium 1 and medium 2
    !   lambda - the wavelength of the photon
    real(8), intent(in) :: direction(3)
    real(8), intent(in) :: origin(3)
    real(8), intent(in) :: radius
    complex(8), intent(in) :: m1, m2
    real(8), intent(in) :: lambda
    ! Outputs:
    !   new_direction - the new direction of the photon after refraction
    !   phase_shift - the phase shift experienced by the photon
    real(8), intent(out) :: new_direction(3)
    real(8), intent(out) :: phase_shift
    ! Local variables:
    !  normal - the normal vector to the sphere at the point of intersection
    !  theta1 - the angle between the photon direction and the normal vector
    !  theta2 - the angle of refraction
    real(8) :: normal(3)
    real(8) :: theta1, sintheta2, costheta2
    real(8) :: k_par(3), k_per(3), k_t(3)
    complex(8) :: r_s, t_s, a_s, r_p, t_p, a_p
    real(8) :: R_s, T_s, A_s, R_p, T_p, A_p
    real(8) :: chance

    ! Find the angle of incidence, which is the angle between the normal vector and the photon direction. 
    normal = origin / radius
    theta1 = acos(dot_product(direction, normal))
    
    ! Calculate the angle of refraction using Snell's law. Here, cos(theta2) is calculated based on the attenuation
    ! of the electric field. The wavevector of light is m cos(theta) which should be positive for decay to occur.
    sintheta2 = m1 / m2 * sin(theta1)
    costheta2 = sqrt(1 - sintheta2**2)
    if imag(m2 * costheta2) < 0 then
      costheta2 = -costheta2
    end if
    ! Calculate the Fresnel coefficients
    r_s = (m1 * cos(theta1) - m2 * costheta2) / (m1 * cos(theta1) + m2 * costheta2)
    r_p = (m2 * cos(theta1) - m1 * costheta2) / (m2 * cos(theta1) + m1 * costheta2)
    t_s = (2 * cos(theta1)) / (m1 * cos(theta1) + m2 * costheta2)
    t_s = (2 * cos(theta1)) / (m2 * cos(theta1) + m1 * costheta2)
    ! Calculate the reflection, transmission, and absorption coefficients
    R_s = abs(r_s)**2
    T_s = m2 * cos(theta2) / (m1 * cos(theta1)) * abs(t_s)**2
    R_p = abs(r_p)**2
    T_p = m2 * cos(theta2) / (m1 * cos(theta1)) * abs(t_p)**2
    A_s = 1 - R_s - T_s
    A_p = 1 - R_p - T_p

    ! Decompose the wave into parallel and perpendicular components and calculate the 
    ! reflection and transmission directions
    k_par = dot_product(direction, normal) * normal
    k_per = direction - k_par
    k_t = 2 * PI / lambda * m2 * (k_par + costheta2 * normal)
    
    call random_number(chance)
    if (chance < R_s) then
      ! Reflect the photon
      new_direction = direction - 2 * k_per
      phase_shift = atan(aimag(r_s) / real(r_s))
    else if (chance < R_s + T_s) then
      ! Transmit the photon
      new_direction = real(k_t) / norm2(real(k_t))
      phase_shift = atan(aimag(r_s) / real(r_s))
      medium = -medium
    else
      ! Absorb the photon
      new_direction = [0.0, 0.0, 0.0]
      phase_shift = 0.0
      intensity = 0.0
    end if

  end subroutine refract_photon
end module ray_tracing