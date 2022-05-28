   module dataset
    implicit none
        ! main program variables
              integer, parameter :: dp = selected_real_kind(15, 307)
              integer :: ngroup,nx,ny,nz,it
              real(dp), allocatable, dimension(:) :: chi
              real(dp), allocatable, dimension(:,:,:) :: Phi
               real(dp), allocatable, dimension(:,:,:) ::  D,sig_R,nu_sgf,sig_abs
              real :: start, finish
              real(dp) :: K_eps, Phi_eps, K ,BK
    end module dataset	
    
      module typedata
    implicit none
     integer, parameter :: dp = selected_real_kind(15, 307)
     
          TYPE geometry
			   integer :: nmat,group,gmtry
               integer, allocatable, dimension(:) :: nreg
               integer, allocatable, dimension(:,:,:) :: core
               real(dp), allocatable, dimension(:,:) :: msize,reglen
          end TYPE geometry
          
          TYPE (geometry) :: geometry_inp
          
          TYPE xsection
              ! real(dp), allocatable, dimension(:) :: chi
			   real(dp), allocatable, dimension(:,:) ::  D1,sig_a1,nu_sgf1
			   real(dp), allocatable, dimension(:,:,:) :: sig_scatt1
          end TYPE xsection
          
          TYPE (xsection) :: xsec_inp
          
          TYPE scattering
              real(dp), allocatable, dimension(:) :: scatt
          end TYPE scattering
          
          TYPE coeffs
              real(dp), allocatable, dimension(:,:) :: aC,aL,aR,aT,aB,F
              TYPE(Scattering), allocatable, dimension(:,:) :: sig_down
              real(dp), allocatable, dimension(:,:,:) :: aC3d,aL3d,aR3d,aT3d,aB3d,aF3d,aBa3d,F3d
              real(dp), allocatable, dimension(:,:,:) :: D3d,sig_R3d
              TYPE(Scattering), allocatable, dimension(:,:,:) :: sig_down3d
          End TYPE coeffs
          
          TYPE (coeffs),allocatable,dimension(:) :: coeff
          
          TYPE boundary
              integer :: Left1D,Right1D
              integer :: Left,Right,Top,Bottom,Front,Back
              real(dp) :: ext_len
          End TYPE boundary
          
          TYPE (boundary) :: B_C 
          
          !interface
       !   subroutine BiCGSTAB       
       !real(dp), allocatable,dimension(:,:),intent(in) :: Sinp
       !real(dp), allocatable,dimension(:,:,:),intent(inout) :: phiprev
       !integer , intent(in) :: g
       !   end subroutine
       !   end interface
    end module typedata
    
                 