  !subroutine CMFD
  !  use dataset, only : ngroup,nx,ny,D,nu_sgf,sig_abs,Phi
  !  implicit none
  !  real(dp), allocatable, dimension(:,:) :: D_bar,nu_sgf_bar,sig_abs_bar,phi_bar
  !  real(dp) :: g,i,j
  !  
  !  allocate (D_bar(ny,nx),nu_sgf_bar(ny,nx),sig_abs_bar(ny,nx),phi_bar(ny,nx))
  !  
  !   D_bar=0.0d0
  !   nu_sgf_bar=0.0d0
  !   sig_abs_bar=0.0d0
  !   phi_bar=0.0d0
  !   
  !   
  !       do i=1,ny
  !           do j=1,nx
  !               do g=1,ngroup
  !                 phi_bar(i,j)=phi_bar(i,j)+Phi(g,i,j)
  !                 D_bar(i,j)=D_bar(i,j)+D(g,i,j)*Phi(g,i,j)
  !                 nu_sgf_bar(i,j)=nu_sgf_bar(i,j)+nu_sgf(g,i,j)*Phi(g,i,j)
  !                 sig_abs_bar(i,j)=sig_abs_bar(i,j)+sig_abs(g,i,j)*Phi(g,i,j)
  !               enddo
  !           enddo
  !       enddo
  !       
         !D_bar=D_bar/phi_bar
         !nu_sgf_bar=nu_sgf_bar/phi_bar
         !sig_abs_bar=sig_abs_bar/phi_bar
         !
         