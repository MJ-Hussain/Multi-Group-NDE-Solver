    ! subroutine SOR   
         subroutine sar(aT,aR,aC,aL,aB,s,phi,po,n)
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 307)
        integer ::  i, j,it,jj 
        real(dp)::tol,dif,sumj,w
        integer, intent(in) :: n
      real(dp), dimension(:,:), intent(in) :: aT,aR,aC,aL,aB,s
!      real(dp), dimension(:,:), intent(in) :: fac
      real(dp), dimension(:,:), intent(out) :: phi
      real(dp), dimension(:,:) :: po
!       if (.not. allocated(po)) allocate(po(n))
      tol=1d-03
       dif=999                                                 
      w=1.7
       it=0 
        do while (dif .GT. tol)
           it=it+1
           po=phi
               phi(1,1)=((s(1,1)-aR(1,1)*po(2,1)-aT(1,1)*po(1,2))*w/aC(1,1) + (1-w)*po(1,1))!*fac(1,1)
               
                do j=2,n-1
                    phi(1,j)=((s(1,j)-aB(1,j)*phi(1,j-1)-aR(1,j)*po(2,j)-aT(1,j)*po(1,j+1))*w/aC(1,j) + (1-w)*po(1,j))!*fac(1,j)
                    enddo
               phi(1,n)=((s(1,n)-aB(1,n)*phi(1,n-1)-aR(1,n)*po(2,n))*w/aC(1,n) + (1-w)*po(1,n))!*fac(1,n)
                   do i=2,n-1
                    phi(i,1)=((s(i,1)-aL(i,1)*phi(i-1,1)-aR(i,1)*po(i+1,1)-aT(i,1)*po(i,2))*w/aC(i,1) + (1-w)*po(i,1))!*fac(i,1)
                    enddo  
               phi(n,1)=((s(n,1)-aL(n,1)*phi(n-1,1)-aT(n,1)*po(n,2))*w/aC(n,1) + (1-w)*po(n,1))!*fac(n,1)
               
                      do i=2,n-1
            do j=2,n-1 
                phi(i,j)=((s(i,j)-aL(i,j)*phi(i-1,j)-aB(i,j)*phi(i,j-1)-aR(i,j)*po(i+1,j)-aT(i,j)*po(i,j+1))*w/aC(i,j) + (1-w)*po(i,j))!*fac(i,j)
             enddo   
          enddo
              
               do j=2,n-1 
                    phi(n,j)=((s(n,j)-aL(n,j)*phi(n-1,j)-aB(n,j)*phi(n,j-1)-aT(n,j)*po(n,j+1))*w/aC(n,j) + (1-w)*po(n,j) )!*fac(n,j)
               enddo
               do i=2,n-1
                    phi(i,n)=((s(i,n)-aL(i,n)*phi(i-1,n)-aB(i,n)*phi(i,n-1)-aR(i,n)*po(i+1,n))*w/aC(i,n) + (1-w)*po(i,n))!*fac(i,n)
               enddo
                  phi(n,n)=((s(n,n)-aL(n,n)*phi(n-1,n)-aB(n,n)*phi(n,n-1))*w/aC(n,n) + (1-w)*po(n,n) )!*fac(n,n)
                   
           !  dif= maxval((phi))-maxval((po))
          !  dif= abs(sum(abs(phi))-sum(abs(po)))
          ! dif=maxval(abs(phi/po-1))
                  dif= maxval(abs((phi-po)/phi))
      enddo 
 !   print*,'SOR it', it
    end subroutine sar