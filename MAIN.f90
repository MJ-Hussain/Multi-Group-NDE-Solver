
    
    program Console21
      use dataset
      use typedata
      implicit none
           
         
      
       
        call read_data
        
        call coeff_calc
        
     call CPU_TIME(start)
        
        call power
        
     call CPU_TIME(finish)
             print*,'CPU time ', finish - start
             
             write(11,'(A,I6,A)') 'Solution converged in ',it,' iterations'
             write(11,'(A,F8.6)') 'Final K-effective ', K
             write(11,'(A,F8.2,A)') 'CPU Time ', finish-start, ' Seconds'
    end program Console21

    

    
     
     