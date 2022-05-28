!    ! subroutine SOR   
!         subroutine sor(aT,aR,aC,aL,aB,s,phi,po,n,g)
!   implicit none
!   integer, parameter :: dp = selected_real_kind(15, 307)
!        integer ::  i, j,it,jj 
!        real(dp)::tol,dif,sumj,w
!        integer, intent(in) :: n,g
!      real(dp), dimension(:,:), intent(in) :: aT,aR,aC,aL,aB,s
!  !    real(dp), dimension(:,:), intent(in) :: fac
!      real(dp),pointer, dimension(:,:), intent(out) :: phio
!      real(dp), dimension(:,:) :: po
!  !    real(dp), dimension(n,n) :: phi1
! !     allocate (phi1(n,n))
!!       if (.not. allocated(po)) allocate(po(n))
!      tol=1d-03
!       dif=999                                                 
!      w=1.7
!       it=0 
!        do while (dif .GT. tol)
!           it=it+1
!            po=phio
!               
!        do i=1,n
!            do j=1,n 
!                phio(i,j)=((s(i,j)-aL(i,j)*phio(i-1,j)-aB(i,j)*phio(i,j-1)-aR(i,j)*po(i+1,j)-aT(i,j)*po(i,j+1))*w/aC(i,j) + (1-w)*po(i,j))
!             enddo   
!          enddo
!                 phio(0,:)=0.0d0
!                 phio(:,n+1)=0.0d0
!                 phio(:,0)=phio(:,1)
!                 phio(n+1,:)=phio(n,:)
!    
!                   
!           !  dif= maxval((phi))-maxval((po))
!          !  dif= abs(sum(abs(phi))-sum(abs(po)))
!          ! dif=maxval(abs(phi/po-1))
!                dif= maxval(abs((phio-po)/phio))   
!                    
!                  
!      enddo 
! !   print*,'SOR it', it
!    end subroutine sor