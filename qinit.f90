! qinit routine for initializing lake into solution, q
! polygon METHOD (with 3 points)
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use geoclaw_module, only: grav

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Lake Height (min=2265.0, max=2280.0)
    real(kind=8), parameter :: h = 2267.d0	

    ! Northwest point (Upper and left bound)
    real(kind=8), parameter :: nw_x = 94.72d0
    real(kind=8), parameter :: nw_y = 30.315d0

    !New Dam
    ! Southwest point on right bank at outlet
    real(kind=8), parameter :: sw_x = 94.9324d0
    real(kind=8), parameter :: sw_y = 30.1752d0
    ! Southeast point on left bank at outlet  (make sure this is high enough in lat.)
    real(kind=8), parameter :: se_x = 94.9651d0
    real(kind=8), parameter :: se_y = 30.1970

    
    ! Other storage
    integer :: i,j
    real(kind=8) :: x,y
    real(kind=8) :: m,b
 
    !Calculate Slope/Intercept
    m = (se_y - sw_y)/(se_x - sw_x)    
    b = sw_y - m*sw_x

    do i=1-mbc,mx+mbc
        x = xlower + (i - 0.5d0)*dx
        do j=1-mbc,my+mbc
            y = ylower + (j - 0.5d0)*dy
	    !Fill x<sw_x, sw_y<y<nw_y
	    if ((x>nw_x).and.(x<sw_x).and.(y>sw_y).and.(y<nw_y)) then
                q(1,i,j) = max(0.d0,(h - aux(1,i,j))) 
	    else if ((x>sw_x).and.(x<se_x).and.(y>(m*x+b)).and.(y<se_y)) then
	        q(1,i,j) = max(0.d0,(h - aux(1,i,j)))
            else
		q(1,i,j) = 0.d0
	    endif
        enddo
    enddo
    
end subroutine qinit
