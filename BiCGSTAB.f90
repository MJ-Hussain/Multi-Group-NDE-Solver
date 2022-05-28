!   subroutine BiCGSTAB (g,Sinp,phiprev)
!     use dataset, only : ngroup,nx,ny,nz,K_eps,Phi_eps ,chi,it,K ,Phi
!                 use typedata
!                 implicit none
!       real(dp), allocatable,dimension(:,:),intent(in) :: Sinp
!       integer , intent(in) :: g
!       real(dp), allocatable,dimension(:,:,:),intent(inout) :: phiprev
!        real    (dp), allocatable,dimension(:,:)   :: r, rs, v, p, s, t,AX
!        real    (dp), parameter                     :: e = 1d-33
!        real    (dp)                                :: rho      , rho_prev
!        real    (dp)                                :: alpha    , omega   , beta
!        real    (dp)                                :: norm_r   , norm_b       
!        real    (dp)                                :: summesion, temp
!
!        integer                                 :: iit=0,err,i,j,l 
!        
!        
!    allocate(r(ny,nx), rs(ny,nx), v(ny,nx), p(ny,nx), s(ny,nx), t(ny,nx),AX(ny,nx))
!        
!    !-------------------------------------------------------!
!        Phi(g,:,:)  = phiprev(g,:,:)                                     !-------> INITIAL GUESS
!!-------------------------------------------------------!
!        do i=1,ny
!            do j=1,nx
!        AX=  coeff(g)%aL(i,j)*phi(g,i-1,j)+coeff(g)%aB(i,j)*phi(g,i,j-1)+coeff(g)%aR(i,j)*Phi(g,i+1,j)+coeff(g)%aT(i,j)*Phi(g,i,j+1)
!            enddo
!        enddo
!        
!        r  = Sinp - AX                            !-------> LINE 1
!        rs = r                                          !
!!-------------------------------------------------------!
!        rho   = 1.0d0; alpha = 1.0d0; omega = 1.0d0  !-------> LINE 2
!!-------------------------------------------------------!
!        v  = 0.0d0; p  = 0.0d0                        !-------> LINE 3
!!                                                       !
!        norm_r = sqrt(sum(r*r))                 !
!        norm_b = sqrt(sum(Sinp*Sinp))                 !
!!-------------------------------------------------------!
!          do while(norm_r .GT. e*norm_b)                          !-------> START OF LOOP
!
!        !-------------------------------------------------------!
!            rho_prev = rho                                      !-------> LINE 5
!            rho      = sum(rs*r)                        !
!        !-------------------------------------------------------!
!            beta     = (rho/rho_prev) * (alpha/omega)           !-------> LINE 6
!        !-------------------------------------------------------!
!            p        = r + beta * (p - omega*v)                 !-------> LINE 7
!        !-------------------------------------------------------!
!           
!            v        =   coeff(g)%aL(i,j)*p(i-1,j)+coeff(g)%aB(i,j)*p(i,j-1)+coeff(g)%aR(i,j)*p(i+1,j)+coeff(g)%aT(i,j)*P(i,j+1)                            !-------> LINE 8
!        !-------------------------------------------------------!
!            alpha    = rho/sum(rs*v)                    !-------> LINE 9
!        !-------------------------------------------------------!
!            s        = r - alpha*v                              !-------> LINE 10
!        !-------------------------------------------------------!
!            t        = coeff(g)%aL(i,j)*s(i-1,j)+coeff(g)%aB(i,j)*s(i,j-1)+coeff(g)%aR(i,j)*s(i+1,j)+coeff(g)%aT(i,j)*s(i,j+1)                               !-------> LINE 11
!        !-------------------------------------------------------!
!            omega    = sum(t*s)/sum(t*t)        !-------> LINE 12
!        !-------------------------------------------------------!
!            phi(g,:,:)        = phi(g,:,:) + alpha*p + omega*s                    !-------> LINE 13
!        !-------------------------------------------------------!
!            r        = s - omega*t                              !-------> LINE 17
!        !-------------------------------------------------------!
!            norm_r   = sqrt(sum(r*r))                   !
!            norm_b   = sqrt(sum(Sinp*Sinp))                   !
!        !-------------------------------------------------------!
!            iit = iit + 1                                         !
!        !-------------------------------------------------------!
!
!          end do 
!          phiprev(g,:,:)=phi(g,:,:)
!end subroutine