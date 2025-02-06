module utils
  implicit none

contains

  subroutine save_results_to_csv(filename, results)
    character(len=*), intent(in) :: filename
    real(8), intent(in) :: results(:, :)
    integer :: i, unit

    open(newunit=unit, file=filename, status='replace', action='write')
    do i = 1, size(results, 1)
      write(unit, '(F8.4, ",", F8.4)') results(i, 1), results(i, 2)
    end do
    close(unit)

  end subroutine save_results_to_csv

end module utils
