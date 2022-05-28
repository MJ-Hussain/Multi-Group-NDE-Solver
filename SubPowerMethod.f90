         subroutine power
    !(coeff,ngroup,nx,ny,chi,K_eps,Phi_eps)
                 use dataset, only : ngroup,nx,ny,nz,K_eps,Phi_eps ,chi,it,K ,Phi
                 use typedata
                 implicit none
                 
                      !TYPE (coeffs),allocatable,dimension(:),intent(in) :: coeff
                      !integer, intent(in) :: ngroup, nx, ny
                      !real(dp), allocatable, dimension(:), intent(in) :: chi
                      !real(dp), intent(in) :: K_eps, Phi_eps
                      integer           :: g,i,j,l,ll
                    !  real(dp), allocatable, dimension(:,:,:) :: Phi 
                      real(dp), allocatable, dimension(:) :: wg ,ut,zt
                      real(dp), allocatable, dimension(:,:)   :: Source,S_prev,S_pp,S_inp,dScatt_sum,uScatt_sum
                      real(dp), allocatable, dimension(:,:,:) :: Source3d,S_prev3d,Phi_prev3d,S_inp3d,dScatt_sum3d,uScatt_sum3d,Phi_prev ,oold_phi
                      real(dp), allocatable,dimension(:,:,:,:) :: Phi3d ,old_phi
                      real(dp)                               :: K_prev,Phi_tol,phi_tol_prev,K_tol,w,muo,alpha_prev,beta_prev,dom,theta,ap,bp
                    !-------------BiCGSTAB ----------------------!
                      real    (dp), allocatable,dimension(:,:)   :: r, rs, v, p, s, t,AX
                      real    (dp), parameter :: e = 1d-33
                      real    (dp) :: rho , rho_prev
                      real    (dp) :: alpha    , omega   , beta
                      real    (dp) :: norm_r   , norm_b       
                      real    (dp) :: summesion, temp

                      integer :: iit=0,err 
                      
    
     !*****************************************************************************************************!
     !----------------------- 1D Calculation --------------------------------------------------------------!
     !_____________________________________________________________________________________________________!
         if (geometry_inp%gmtry==1) then
             !print*, geometry_inp%gmtry;pause
           allocate(Source(1,nx),Phi(ngroup,1,0:nx+1),S_prev(1,nx),Phi_prev(ngroup,1,0:nx+1),dScatt_sum(1,nx),uScatt_sum(1,nx),S_inp(1,nx),wg(ngroup))
                  Phi_prev=0.0
                    Phi=1.0
                   ! print*, Phi(1,0,:);stop
                     Phi(:,1,0)=0.0d0
                     Phi(:,1,nx+1)=0.0d0
                     
                    ! print*, Phi(1,1,:);stop
                    Source=1.0d0
                    K_tol=999
                !    K_eps=1d-07
                 !   Phi_eps=1d-03
                    K=1.0d0
       ! Relaxation parameter calcualtion
                    do g=1,ngroup
                   muo= sum(abs(coeff(g)%aL+coeff(g)%aR))/sum(coeff(g)%aC)
                   wg(g)=2/(1+sqrt(1-muo))
                   enddo
                   ! w=2/(1+sqrt(1-muo))  
                   ! print*, muo,wg; pause
                    it=0    
        do while (K_tol>K_eps)
                it=it+1
                S_prev(1,1:nx)=Source(1,1:nx)
                K_prev=K
            do g=1,ngroup
                     dScatt_sum=0.0d0
                     uScatt_sum=0.0d0
                do l=1,g-1
                   do j=1,nx    
                        dScatt_sum(1,j)=dScatt_sum(1,j)+coeff(g)%sig_down(1,j)%scatt(l)*Phi(l,1,j)
                   enddo
                enddo
                do l=g+1,ngroup
                   do j=1,nx    
                        uScatt_sum(1,j)=uScatt_sum(1,j)+coeff(g)%sig_down(1,j)%scatt(l)*Phi(l,1,j)
                   enddo
                enddo
                !     do l=1,ngroup
                !       do j=1,nx    
                !         dScatt_sum(1,j)=dScatt_sum(1,j)+coeff(g)%sig_down(1,j)%scatt(l)*Phi(l,1,j)
                !   enddo
                !enddo    
                     
              ! print*, coeff(2)%sig_down(1,1)%scatt(:);stop
               !print*, uScatt_sum(1,:)   ;pause
              ! print*, dScatt_sum(1,:)   ;pause                                                         
                      
                   S_inp(1,1:nx)=Source(1,1:nx)*chi(g)/k+dScatt_sum(1,1:nx) +uScatt_sum(1,1:nx)
                    ! print*, S_inp(1,:);pause
                       Phi_tol=999
                       do while (Phi_tol>Phi_eps)
                         Phi_prev(g,1,0:nx+1)=Phi(g,1,0:nx+1)
                     
                           do j=1,nx
                        
                            phi(g,1,j)=((S_inp(1,j)-coeff(g)%aL(1,j)*phi(g,1,j-1)-coeff(g)%aR(1,j)*Phi_prev(g,1,j+1))*wg(g)/coeff(g)%aC(1,j) + (1-wg(g))*Phi_prev(g,1,j))
                       
                           enddo
                       !  phi_tol_prev=phi_tol
                         Phi_tol= maxval(abs((phi(g,1,0:nx+1)-phi_prev(g,1,0:nx+1))/phi(g,1,0:nx+1)))
                       enddo
            enddo
                     !      call CMFD
               !print*, phi(1,1,1:nx);pause
               ! print*, phi(2,1,1:nx);pause
                Source=0.0d0
                
                do i=1,ngroup
              Source(1,1:nx)=Source(1,1:nx)+coeff(i)%F(1,1:nx)*Phi(i,1,1:nx)
             ! print*,coeff(i)%F(1,1:nx)
                enddo
               ! print*, Source(1,:);pause
            ! print*, Source(61,:);pause   
           K=K_prev*sum(Source(1,1:nx))/sum(S_prev(1,1:nx))
            K_tol=abs((K-K_prev)/K_prev)
            print*,'k ',k ,'      tol',K_tol
     enddo
          open (unit=22, file='Flux.txt')
          do g=1,ngroup
              write(22,'(A,I3)') 'Flux Group ',g
              do i=1,nx
                  write(22,'(ES13.6)') Phi(g,1,i)
              enddo
          enddo
          close(22)
                     
                      
      !*****************************************************************************************************!
      !----------------------- 2D Calculation --------------------------------------------------------------!
      !_____________________________________________________________________________________________________!
                 elseif (geometry_inp%gmtry==2) then
       allocate(Source(ny,nx),Phi(ngroup,0:ny+1,0:nx+1),S_prev(ny,nx),Phi_prev(ngroup,0:ny+1,0:nx+1),dScatt_sum(ny,nx),uScatt_sum(ny,nx),S_inp(ny,nx),wg(ngroup))
       allocate(old_phi(ngroup,1,ny,nx),ut(ngroup),zt(ngroup),S_pp(ny,nx),oold_phi(ngroup,ny,nx))
        allocate(r(ny,nx), rs(ny,nx), v(ny,nx), p(0:ny+1,0:nx+1), s(0:ny+1,0:nx+1), t(ny,nx),AX(ny,nx)) 
                   Phi_prev=1.0
                    Phi=1.0
                   ! print*, Phi(1,0,:);stop
                     Phi(:,0,:)=0.0d0
                     Phi(:,:,nx+1)=0.0d0
                     Phi(:,:,0)=0.0d0
                     Phi(:,ny+1,:)=0.0d0
                     Phi_prev=Phi
                    Source=1.0d0
                    S_prev=0.0d0
                    old_phi=0.0d0
                      alpha=1.0d0
                      beta=0.5d0
                      theta=1.0d0
                    K_tol=999
                !    K_eps=1d-07
                 !   Phi_eps=1d-03
                    K=1.0d0
                    ! Relxation parameter calclation
                    do g=1,ngroup
                    muo= sum(abs(coeff(g)%aL+coeff(g)%aR+coeff(g)%aT+coeff(g)%aB))/sum(coeff(g)%aC)
                    wg(g)=2/(1+sqrt(1-muo))
                    enddo
                    !  print*,muo, wg ;pause
                      ap=1.0d0; bp=0.0d0  
                    !do g=1,ngroup
                    !do i=1,ny
                    !    do j=1,nx
                    !      muo=abs(coeff(g)%aL(i,j))/coeff(g)%aC(i,j)+abs(coeff(g)%aR(i,j))/coeff(g)%aC(i,j)+abs(coeff(g)%aT(i,j)/coeff(g)%aC(i,j))+abs(coeff(g)%aB(i,j))/coeff(g)%aC(i,j)  
                    !enddo
                    !enddo
                    !muo=muo/nx
                    ! wg(g)=2/(1+sqrt(1-muo))
                    ! print*,muo, wg ;pause
                    !enddo
                    !  stop  
                    
                    
                    
                    it=0
                !do i=1,ngroup
                ! Source=Source+coeff(i)%F*Phi(i,1:ny,1:nx)
                !enddo   
                open(unit=23, file='comparison.txt')
                write(23,'(A)')  ' it     K_tol'
             do while (K_tol>K_eps)
                it=it+1
                S_pp=S_prev
                S_prev=Source
                K_prev=K

      
                do g=1,ngroup
                     dScatt_sum=0.0d0
                     uScatt_sum=0.0d0
                     
                do l=1,g-1
                    do i=1,ny
                        do j=1,nx
                            
                                dScatt_sum(i,j)=dScatt_sum(i,j)+coeff(g)%sig_down(i,j)%scatt(l)*Phi(l,i,j)
                            enddo
                        enddo
                enddo
                do l=g+1,ngroup
                    do i=1,ny
                        do j=1,nx
                            
                                uScatt_sum(i,j)=uScatt_sum(i,j)+coeff(g)%sig_down(i,j)%scatt(l)*Phi(l,i,j)
                            enddo
                        enddo
                enddo
                    
                    S_inp=Source*chi(g)/K+dScatt_sum+uScatt_sum
                    phi_tol=1.0d0
!                    !call BiCGSTAB (g,S_inp,phi_prev)
!     !--------------------------------------------------------!
!     !                       BiCGSTAB                         !
!     !--------------------------------------------------------!                     
!                   
!    !-------------------------------------------------------!
!        Phi(g,:,:)  = phi_prev(g,:,:)                                     !-------> INITIAL GUESS
!    !-------------------------------------------------------!
!        do i=1,ny
!            do j=1,nx
!        AX(i,j)=  coeff(g)%aL(i,j)*phi(g,i-1,j)+coeff(g)%aB(i,j)*phi(g,i,j-1)+coeff(g)%aR(i,j)*Phi(g,i+1,j)+coeff(g)%aT(i,j)*Phi(g,i,j+1)
!            enddo
!        enddo
!        
!        r  = abs(S_inp - AX) !-------> LINE 1
!        rs = r 
!        print*, rs(1,1:10) ;pause
!!-------------------------------------------------------!
!        rho   = 1.0d0; alpha = 1.0d0; omega = 1.0d0  !-------> LINE 2
!!-------------------------------------------------------!
!        v  = 0.0d0; p  = 0.0d0                        !-------> LINE 3
!!                                                       !
!        !norm_r = sqrt(sum(r*r))                 !
!        !norm_b = sqrt(sum(S_inp*S_inp))                 !
!!-------------------------------------------------------!
!          do while(Phi_tol>Phi_eps)                         !-------> START OF LOOP
!
!        !-------------------------------------------------------!
!            rho_prev = rho                                      !-------> LINE 5
!            rho      = sum(rs*r)                        !
!        !-------------------------------------------------------!
!            beta     = (rho/rho_prev) * (alpha/omega)           !-------> LINE 6
!        !-------------------------------------------------------!
!            p(1:ny,1:nx)        = r + beta * (p(1:ny,1:nx) - omega*v)                 !-------> LINE 7
!        !-------------------------------------------------------!
!            do i=1,ny
!                do j=1,nx
!            v(i,j)        =   coeff(g)%aL(i,j)*p(i-1,j)+coeff(g)%aB(i,j)*p(i,j-1)+coeff(g)%aR(i,j)*p(i+1,j)+coeff(g)%aT(i,j)*P(i,j+1)  !-------> LINE 8
!                enddo
!            enddo
!            
!        !-------------------------------------------------------!
!            alpha    = rho/sum(rs*v)                    !-------> LINE 9
!        !-------------------------------------------------------!
!            s(1:ny,1:nx)        = r - alpha*v                              !-------> LINE 10
!        !-------------------------------------------------------!
!          do i=1,ny
!              do j=1,nx
!            t (i,j)       = coeff(g)%aL(i,j)*s(i-1,j)+coeff(g)%aB(i,j)*s(i,j-1)+coeff(g)%aR(i,j)*s(i+1,j)+coeff(g)%aT(i,j)*s(i,j+1)                               !-------> LINE 11
!              enddo
!          enddo
!          
!            !-------------------------------------------------------!
!            omega    = sum(t*s(1:ny,1:nx))/sum(t*t)        !-------> LINE 12
!        !-------------------------------------------------------!
!            do i=1,ny
!                do j=1,nx
!            phi(g,i,j)        =abs( phi(g,i,j) + alpha*p(i,j) + omega*s(i,j) )
!                enddo
!            enddo
!            
!        !-------------------------------------------------------!
!            r        = abs(s(1:ny,1:nx) - omega*t )                             !-------> LINE 17
!        !-------------------------------------------------------!
!            norm_r   = sqrt(sum(r*r))                   !
!            norm_b   = sqrt(sum(S_inp*S_inp))                   !
!        !-------------------------------------------------------!
!            iit = iit + 1                                         !
!        !-------------------------------------------------------!
!             Phi_tol= maxval(abs((phi(g,:,:)-phi_prev(g,:,:))/phi(g,:,:)))
!             phi_prev(g,:,:)=phi(g,:,:) 
!           !  print*, Phi_tol, norm_r, norm_b;pause
!          end do 
!                   
!                    
!           print*,  phi(g,1,1:10); pause         
!                    
                    
                    
                    
     !--------------------------------------------------------!
     !                        SOR                             !
     !--------------------------------------------------------!  
                       
                   do while (Phi_tol>Phi_eps)
                          Phi_prev(g,:,:)=Phi(g,:,:)
                     do i=1,ny
                       do j=1,nx
                        
                        phi(g,i,j)=((S_inp(i,j)-coeff(g)%aL(i,j)*phi(g,i-1,j)-coeff(g)%aB(i,j)*phi(g,i,j-1)-coeff(g)%aR(i,j)*Phi_prev(g,i+1,j)-coeff(g)%aT(i,j)*Phi_prev(g,i,j+1))*wg(g)/coeff(g)%aC(i,j) + (1-wg(g))*Phi_prev(g,i,j))
                       enddo
                     enddo
                      phi_tol_prev=phi_tol 
                      Phi_tol= maxval(abs((phi(g,:,:)-phi_prev(g,:,:))/phi(g,:,:)))
                         
                   enddo
     !---------------------------------END SOR-------------------------!              
                   ! ut(g)=  phi_tol/phi_tol_prev
                   ! zt(g)=ut(g)/(1-ut(g))
                enddo
                   
                   
                 !if (it<20)then
                 !    zt=1.5d0
                 !endif
                Source=0.0d0
                
              
                do i=1,ngroup
                   ! Phi(i,1:ny,1:nx)=old_Phi(i,1,:,:)+ap*( Phi(i,1:ny,1:nx)-old_Phi(i,1,:,:))+bp*(old_phi(i,1,:,:)-oold_phi(i,:,:))
                    !Phi(i,1:ny,1:nx)=old_Phi(i,1,:,:)+1.65*( Phi(i,1:ny,1:nx)-old_Phi(i,1,:,:))
                     !Phi(i,1:ny,1:nx)=Phi(i,1:ny,1:nx)+0.5*( Phi(i,1:ny,1:nx)-old_Phi(i,1,:,:))
              Source=Source+coeff(i)%F*Phi(i,1:ny,1:nx)
                  !oold_phi(i,:,:)= old_phi(i,1,:,:)
                  !old_phi(i,1,:,:)=Phi(i,1:ny,1:nx)
              
                enddo
               !   alpha_prev=alpha
               !beta_prev=beta
               !alpha=maxval(K_prev*Source/S_prev)
               !beta=minval(K_prev*Source/S_prev)
               !dom=(alpha-beta)/(alpha_prev-beta_prev)
               !ap=1/((1-0.5d0*dom)-(dom*dom*ap/16.0))
               !bp=(1-0.5d0*dom)*ap-1
               ! if (it>25) then
               ! Source=S_prev+ap*(Source-S_prev)+bp*(S_prev-S_pp)
               ! endif
                
           K=K_prev*sum(Source)/sum(S_prev)
           
               
               !if (it>10) then
               !    theta=2/(2-dom)
               !endif
              ! print*,dom,ap,bp 
            K_tol=abs((K-K_prev)/K_prev)
            print*,'k ',k ,'      tol',K_tol
             
            write(23,'(I3,ES12.4)')  it, K_tol
             enddo
              close(23)
             
      !*****************************************************************************************************!
      !----------------------- 3D Calculation --------------------------------------------------------------!
      !_____________________________________________________________________________________________________!
                 elseif (geometry_inp%gmtry==3) then 
                 allocate(Source3d(nz,ny,nx),S_prev3d(nz,ny,nx),dScatt_sum3d(nz,ny,nx),uScatt_sum3d(nz,ny,nx),S_inp3d(nz,ny,nx),wg(ngroup))
                 allocate(Phi3d(ngroup,0:nz+1,0:ny+1,0:nx+1),Phi_prev3d(0:nz+1,0:ny+1,0:nx+1))   
                ! Phi_prev3d=10.0
                 !do g=1,ngroup
                 !   Phi3d(g,1:nz,1:ny,1:nx)= coeff(g)%F3d(:,:,:)
                 !enddo
                  Phi3d=1.0d0
                   ! print*, Phi(1,0,:);stop
                     Phi3d(:,:,0,:)=0.0d0
                     Phi3d(:,:,:,nx+1)=0.0d0
                     Phi3d(:,:,:,0)=0.0d0
                     Phi3d(:,:,ny+1,:)=0.0d0
                     Phi3d(:,0,:,:)=0.0d0
                     Phi3d(:,nz+1,:,:)=0.0d0
                     !print*, Phi3d(1,0,:,:);stop
                    Source3d=1.0d0
                    K_tol=999
                !    K_eps=1d-07
                 !   Phi_eps=1d-03
                    K=1.0d0
                    
                       ! Relaxation parameter calcualtion 
                    do g=1,ngroup
                   muo= sum(abs(coeff(g)%aL3d+coeff(g)%aR3d+coeff(g)%aT3d+coeff(g)%aB3d+coeff(g)%aF3d+coeff(g)%aBa3d))/sum(coeff(g)%aC3d)
                    wg(g)=2/(1+sqrt(1-muo))
                   enddo
                   
                    !  print*, muo,w ;pause
                   ! w=1.5
                    it=0     
              do while (K_tol>K_eps)
                  it=it+1
                  S_prev3d=Source3d
                  K_prev=K

          !print*,phi(1,1,1:30);pause      
                do g=1,ngroup
                     dScatt_sum3d=0.0d0
                     uScatt_sum3d=0.0d0
                     
                do ll=1,g-1
                    do l=1,nz
                      do i=1,ny
                        do j=1,nx
                            
                          dScatt_sum3d(l,i,j)=dScatt_sum3d(l,i,j)+coeff(g)%sig_down3d(l,i,j)%scatt(ll)*Phi3d(ll,l,i,j)
                        enddo
                      enddo
                    enddo
                enddo
                 do ll=g+1,ngroup
                    do l=1,nz
                      do i=1,ny
                        do j=1,nx
                            
                         uScatt_sum3d(l,i,j)=uScatt_sum3d(l,i,j)+coeff(g)%sig_down3d(l,i,j)%scatt(ll)*Phi3d(ll,l,i,j)
                        enddo
                      enddo
                    enddo 
                 enddo
                    !print*, coeff(1)%sig_down(1,1)%scatt(:)  ;stop
                    !print*,dScatt_sum(1,1:10)
                    ! print*,uScatt_sum(1,1:10);pause
                    S_inp3d=Source3d*chi(g)/K+dScatt_sum3d+uScatt_sum3d  
                    
                         Phi_tol=999
                   do while (Phi_tol>Phi_eps)
                       Phi_prev3d=Phi3d(g,:,:,:)
                     do l=1,nz
                      do i=1,ny
                       do j=1,nx
                        
                        Phi3d(g,l,i,j)=((S_inp3d(l,i,j)-coeff(g)%aL3d(l,i,j)*Phi3d(g,l,i-1,j)-coeff(g)%aB3d(l,i,j)*Phi3d(g,l,i,j-1)-coeff(g)%aBa3d(l,i,j)*Phi3d(g,l-1,i,j)-coeff(g)%aR3d(l,i,j)*Phi_prev3d(l,i+1,j)-coeff(g)%aT3d(l,i,j)*Phi_prev3d(l,i,j+1)-coeff(g)%aF3d(l,i,j)*Phi_prev3d(l+1,i,j))*wg(g)/coeff(g)%aC3d(l,i,j) + (1-wg(g))*Phi_prev3d(l,i,j))
                       enddo
                      enddo
                    enddo
                       
                      Phi_tol= maxval(abs((Phi3d(g,:,:,:)-phi_prev3d)/phi3d(g,:,:,:)))
                  enddo
                enddo
                Source3d=0.0d0
                
                do i=1,ngroup
              Source3d=Source3d+coeff(i)%F3d*Phi3d(i,1:nz,1:ny,1:nx)
              !print*, coeff(i)%F(11,:);pause
                enddo
                !Source3d=Source3d+0.5*(Source3d-S_prev3d)
             !print*, Source(61,:);pause   
           K=K_prev*sum(Source3d)/sum(S_prev3d)
            K_tol=abs((K-K_prev)/K_prev)
            print*,'k ',k ,'      tol',K_tol
             enddo
            
                 endif 
           !open (unit=2, file='flux.txt')
           !do g=1,ngroup
           !    write(2,'(a,I3)') 'Group',g
           !    do l=1,nz
           !        write(2,'(a,I3)') 'Plane',l
           !        do i=1,ny
           !            write(2,'(200F10.6)') (Phi3d(g,l,i,j), j=1,nx)
           !        enddo
           !    enddo
           !enddo
         !   open (unit=2, file='flux.txt')
         !  do g=1,ngroup
         !      write(2,'(a,I3)') 'Group',g
         !    !  do l=1,nz
         !       !   write(2,'(a,I3)') 'Plane',l
         !          do i=1,ny
         !              write(2,'(200ES12.5)') (Phi3d(g,20,ny-i+1,j), j=1,nx)
         !          enddo
         !    !  enddo
         !  enddo       
         !  
         !close(2)          
         ! open(unit=3, file='AxialFlux.txt')
         ! write(3,'(a)') 'Axial Flux CODE'
         ! do g=1,ngroup
         !     write(3,'(a,I3)') 'Group',g
         ! do l=1,nz
         !     write(3,'(ES12.5)') Phi3d(g,l,ny,1)
         ! enddo
         ! enddo
         ! close(3)
         ! 
         ! open(unit=3, file='DiagonalFlux.txt')
         ! write(3,'(a)') 'Diagonal Flux CODE'
         ! do g=1,ngroup
         !     write(3,'(a,I3)') 'Group',g
         ! do i=1,ny
         !     write(3,'(ES12.5)') Phi3d(g,20,i,i)
         ! enddo
         ! enddo
         ! close(3)
             print*, 'k-------------'
            print*, k
            print*,'Number of iterations'
            print*, it
             
          
           
                     
       end subroutine power