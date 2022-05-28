   subroutine read_data
    !(geometry_inp,xsec_inp,chi,B_C,K_eps,Phi_eps)
              use dataset, only : K_eps,Phi_eps,chi,BK
		      use typedata
              implicit none
                 !TYPE(geometry), intent(out) :: geometry_inp
                 !TYPE(xsection), intent(out) :: xsec_inp
                 !TYPE (boundary), intent(out) :: B_C
                 !real(dp), allocatable, dimension(:), intent(out) :: chi
                 !real(dp), intent(out) :: K_eps, Phi_eps
			 integer :: ios,nm,g,i,j,k,l,max_nreg
             integer, allocatable, dimension(:) :: mat
             character(len=150) :: line,title
             character(len=20) :: word,nword
		     
			 ios=0
		    open (unit=1, STATUS='OLD', file='inp.txt')
		    read(1,'(A)') title
			 read(1,'(A)',iostat=ios) line
            read(line,*) word
                    if (word == 'K_tolerence') then
                      read(line,*) word,K_eps
                    else
                        print*,'No K_tolerence input'
                        print*,'using default value 1d-06'
                        K_eps=1d-06
                    endif
            read(1,'(A)',iostat=ios) line
            read(line,*) word
                    if (word == 'flux_tolerence') then
                      read(line,*) word,Phi_eps
                    else
                        print*,'No flux_tolerence'
                        print*,'using default value 1d-03'
                        Phi_eps=1d-03
                    endif
			    ios=0
        do while (ios==0)
            
            read(1,'(A)',iostat=ios) line        
         if (ios==0) then
              select case (line)  
              case ('[GEOMETRY]') 
                read(1,'(A)',iostat=ios) line
                read(line,'(9x,a20)') nword
         !______________________________________________________________________________!       
		 !-------------------  FOR THE CASE OF 1D PROBLEM ------------------------------!
         !______________________________________________________________________________!
                
              if (nword == '1D_Slab') then
                  !print*, nword ;stop
                  geometry_inp%gmtry=1
                  allocate(geometry_inp%nreg(geometry_inp%gmtry))
                  
                  read(1,'(a)',iostat=ios) line
                  read(line,*) word
                    if (word == 'n-regions') then
                      read(line,*) word,geometry_inp%nreg
                      !print*, geometry_inp%nreg;stop
                       max_nreg=maxval(geometry_inp%nreg)
                      ! print*,max_nreg;stop 
                      allocate (geometry_inp%reglen(1,max_nreg),geometry_inp%msize(1,max_nreg))
                      allocate (geometry_inp%core(1,1,max_nreg))
                     
                    else
                        print*,'no regions'
                        stop
                    endif  
                  
                       read(1,'(a)',iostat=ios) line  
                      read(line,*) word
                    if (word == 'reg-length') then
                      read(line,*) word,geometry_inp%reglen(1,:)
					 
                     !   print*,'reg-length',geometry_inp%reglen(2,:) ;stop
                    else
                        print*,'no reg-length'
                        stop
                    endif
                      
                      read(1,'(a)',iostat=ios) line  
                      read(line,*) word
                    if (word == 'mesh-size') then
                      read(line,*) word,geometry_inp%msize(1,:)
					  
                  !     print*,'mesh-size',geometry_inp%msize(2,:); stop
                    else
                        print*,'no mesh-size'
                        stop
                    endif
                    
                       read(1,'(a)',iostat=ios) line
                       read(line,*) word
                    if (word == 'core') then
                    
                         read(1,'(a)',iostat=ios) line
                          read(line,*) (geometry_inp%core(1,1,i), i=1,max_nreg)
                    
!			nmat=maxval(core)
                    else
                        print*,'no core configuration'
                        stop
                    endif
                    
                     read(1,'(a)',iostat=ios) line
                       read(line,*) word
                    if (word == 'boundary') then
                       read(line,*) word,B_C%Left,B_C%Right
                    else
                        print*, 'No boundary condition'
                        stop
                    endif
                    
                    
                    
                    
                    
         !____________________________________________________________________________!         
         !  FOR THE CASE OF 2D PROBLEM
         !____________________________________________________________________________!           
              elseif (nword == '2D_Slab') then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			      
                  geometry_inp%gmtry=2
                 
                  allocate(geometry_inp%nreg(geometry_inp%gmtry))
                  
                  read(1,'(a)',iostat=ios) line
                  read(line,*) word
                    if (word == 'n-regions') then
                      read(line,*) word,geometry_inp%nreg
                       max_nreg=maxval(geometry_inp%nreg)
                      allocate (geometry_inp%reglen(2,max_nreg),geometry_inp%msize(2,max_nreg))
                      allocate (geometry_inp%core(1,geometry_inp%nreg(2),geometry_inp%nreg(1)))
!                       print*,'no-regions',nreg
                    else
                        print*,'no regions'
                        stop
                    endif
                    
                      read(1,'(a)',iostat=ios) line  
                      read(line,*) word
                    if (word == 'reg-length') then
                      read(line,*) word,geometry_inp%reglen(1,:)
					  read(1,'(a)',iostat=ios) line  
					  read(line,*) geometry_inp%reglen(2,:)
                     !   print*,'reg-length',geometry_inp%reglen(2,:) ;stop
                    else
                        print*,'no reg-length'
                        stop
                    endif
                    
                    
                      read(1,'(a)',iostat=ios) line  
                      read(line,*) word
                    if (word == 'mesh-size') then
                      read(line,*) word,geometry_inp%msize(1,:)
					  read(1,'(a)',iostat=ios) line 
                      read(line,*) geometry_inp%msize(2,:)
                  !     print*,'mesh-size',geometry_inp%msize(2,:); stop
                    else
                        print*,'no mesh-size'
                        stop
                    endif

                    
                       read(1,'(a)',iostat=ios) line
                       read(line,*) word
                    if (word == 'core') then
                    do i=1,geometry_inp%nreg(2)
                         read(1,'(a)',iostat=ios) line
                          read(line,*) (geometry_inp%core(1,i,j),j=1,geometry_inp%nreg(1))
                    enddo
!			nmat=maxval(core)
                    else
                        print*,'no core configuration'
                        stop
                    endif
                    
                     read(1,'(a)',iostat=ios) line
                       read(line,*) word
                    if (word == 'buckling') then
                       read(line,*) word,BK
                    else
                        print*, 'No buckling, value set to zero'
                        BK=0.0d0
                    endif
                    
                     read(1,'(a)',iostat=ios) line
                       read(line,*) word
                    if (word == 'boundary') then
                       read(line,*) word,B_C%Left,B_C%Right,B_C%Top,B_C%Bottom
                    else
                        print*, 'No boundary condition'
                        stop
                    endif
                   
                   if ((B_C%Left==2).OR.(B_C%Right==2).OR.(B_C%Top==2).OR.(B_C%Bottom==2)) then
                       read(1,'(a)',iostat=ios) line
                       read(line,*) word
                    if (word == 'extp-length') then
                       read(line,*) word,B_C%ext_len
                    else
                        print*, 'Using default value 0.4692'
                           B_C%ext_len=0.4692d0
                           
                    endif 
                  endif
         !____________________________________________________________________________!         
         !  FOR THE CASE OF 3D PROBLEM
         !____________________________________________________________________________!
              elseif (nword == '3D_Slab') then 
			      
                  geometry_inp%gmtry=3
                 
                  allocate(geometry_inp%nreg(geometry_inp%gmtry))
                  
                  read(1,'(a)',iostat=ios) line
                  read(line,*) word
                    if (word == 'n-regions') then
                      read(line,*) word,geometry_inp%nreg
                       max_nreg=maxval(geometry_inp%nreg)
                      allocate (geometry_inp%reglen(3,max_nreg),geometry_inp%msize(3,max_nreg))
                      allocate (geometry_inp%core(geometry_inp%nreg(3),geometry_inp%nreg(2),geometry_inp%nreg(1)))
!                       print*,'no-regions',nreg
                    else
                        print*,'no regions'
                        stop
                    endif
                    
                      read(1,'(a)',iostat=ios) line  
                      read(line,*) word
                    if (word == 'reg-length') then
                      read(line,*) word,geometry_inp%reglen(1,:)
					  read(1,'(a)',iostat=ios) line  
					  read(line,*) geometry_inp%reglen(2,:)
                      read(1,'(a)',iostat=ios) line  
					  read(line,*) geometry_inp%reglen(3,1:geometry_inp%nreg(3))
                        !print*,'reg-length',geometry_inp%reglen(3,:) ;stop
                    else
                        print*,'no reg-length'
                        stop
                    endif
                    
                    
                      read(1,'(a)',iostat=ios) line  
                      read(line,*) word
                    if (word == 'mesh-size') then
                      read(line,*) word,geometry_inp%msize(1,:)
					  read(1,'(a)',iostat=ios) line 
                      read(line,*) geometry_inp%msize(2,:)
                      read(1,'(a)',iostat=ios) line 
                      read(line,*) geometry_inp%msize(3,1:geometry_inp%nreg(3))
                  !     print*,'mesh-size',geometry_inp%msize(2,:); stop
                    else
                        print*,'no mesh-size'
                        stop
                    endif

                    
                       read(1,'(a)',iostat=ios) line
                       read(line,*) word
                    if (word == 'core') then
                  do l=1,geometry_inp%nreg(3)      
                    do i=1,geometry_inp%nreg(2)
                         read(1,'(a)',iostat=ios) line
                          read(line,*) (geometry_inp%core(l,i,j),j=1,geometry_inp%nreg(1))
                    enddo
                    read(1,'(a)',iostat=ios) line
                  enddo
                  
!			nmat=maxval(core)
                    else
                        print*,'no core configuration'
                        stop
                    endif
                     read(1,'(a)',iostat=ios) line
                       read(line,*) word
                    if (word == 'boundary') then
                       read(line,*) word,B_C%Left,B_C%Right,B_C%Top,B_C%Bottom,B_C%Front,B_C%Back
                    else
                        print*, 'No boundary condition'
                        stop
                    endif
                    
              
              
              
              
              endif	
       !*************************************************************************************!
       !----------------END GEOMETRY---------------------------------------------------------!
       !*************************************************************************************! 
              
				case ('[MATERIALS]') 
				  read(1,'(a)',iostat=ios) line  
                      read(line,*) word
                    if (word == 'n-materials') then
					  read(line,*) word,geometry_inp%nmat
					else
					  print*, 'No number of materials'
					  stop
					endif
					
					
				  read(1,'(a)',iostat=ios) line  
                      read(line,*) word
                    if (word == 'n-groups') then
                      read(line,*) word,geometry_inp%group
					else
					  print*, 'No number of groups'
					  stop
                    endif
					g=geometry_inp%group
                    nm=geometry_inp%nmat
                   ! print*,g, nm; stop
					allocate(xsec_inp%D1(g,nm),xsec_inp%sig_a1(g,nm))
                    allocate(xsec_inp%nu_sgf1(g,nm),xsec_inp%sig_scatt1(nm,g,g))
                    allocate(mat(nm),chi(g))
					
		          read(1,'(a)',iostat=ios) line  
                      read(line,*) word
                    if (word == 'cross-sections') then
					 do i=1,nm
					   read(1,'(a)',iostat=ios) line
					   read(line,*) mat(i),xsec_inp%D1(1,i),xsec_inp%sig_a1(1,i),xsec_inp%nu_sgf1(1,i),(xsec_inp%sig_scatt1(i,1,k), k=1,g)
					   do j=2,g
					     read(1,'(a)',iostat=ios) line
     					     read(line,*) xsec_inp%D1(j,i),xsec_inp%sig_a1(j,i),xsec_inp%nu_sgf1(j,i),(xsec_inp%sig_scatt1(i,j,k), k=1,g)
					   enddo
                     enddo
                    else
                        print*, 'No Cross-sections Found'
                    endif
                   read(1,'(a)',iostat=ios) line  
                      read(line,*) word
                    if (word == 'spectrum') then 
                    read(line,*) word,chi
                    else
                        print*,'No Spectrum Found'
                    endif
                    
			end select
        endif
			
        enddo
            
            
		   open (unit=11, file='out.txt')
        write(11,*),'---------------------------------------------------------'
        write(11,*),'|                        OUTPUT                         |'
        write(11,*),'---------------------------------------------------------'
        write(11,*),title
        write(11,'(A,3I3)'),'Number of regions',geometry_inp%nreg
        write(11,'(A)'),'Region Length (cm)'
		do i=1,geometry_inp%gmtry
		  write(11,'(10(2X,F6.2))') geometry_inp%reglen(i,:)
        enddo		  
        write(11,'(A)'),'Mesh size (cm)'
		do i=1,geometry_inp%gmtry
		  write(11,'(10(2X,F6.2))') geometry_inp%msize(i,:)
        enddo
        write(11,'(A,I3)'),'Number of Materials',nm
		write(11,'(A,I3)'),'Number of Groups',g
        write(11,'(A,ES8.1)'),'Eigen value Tolerence ', K_eps
        write(11,'(A,ES8.1)'),'Flux Tolerence ', Phi_eps
        write(11,*),'Cross Sections'
        do i=1,nm
		   write(11,'(I3,10(2X,F10.7))'), mat(i),xsec_inp%D1(1,i),xsec_inp%sig_a1(1,i),xsec_inp%nu_sgf1(1,i),(xsec_inp%sig_scatt1(i,1,k), k=1,g)
			do j=2,g
			   write(11,'(3X,10(2X,F10.7))'), xsec_inp%D1(j,i),xsec_inp%sig_a1(j,i),xsec_inp%nu_sgf1(j,i),(xsec_inp%sig_scatt1(i,j,k), k=1,g)
			enddo
        enddo
         write(11,*),'Spectrum'
         write(11,'(10F6.2)') chi 
         
 !       close(11)
		end subroutine read_data