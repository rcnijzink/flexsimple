MODULE mo_multi_obj

IMPLICIT NONE

CONTAINS

	subroutine multi_obj(par_ini, par_max,par_min,incon, prec_data,temp_data,etp_data,&
                             dem, cellsize,Qobs_data, dates_data, optim, sumax_data)
		
USE mo_model
USE mo_readdata
USE mo_objectives
USE mo_eval_signatures
USE mo_init_random_seed

	real*8,dimension(4), intent(in)		                :: incon	 ! initial conditions states
	real*8,dimension(:), allocatable, intent(in)	    :: etp_data	 ! evaporation data
	real*8,dimension(:), allocatable, intent(in)	    :: prec_data ! precipitation data
	real*8,dimension(:), allocatable, intent(in)	    :: temp_data ! temperature data
        real*8, dimension(:,:), allocatable, intent(in) :: dem
    real*8,intent(in)                                   :: cellsize
	real*8,dimension(:), allocatable, intent(in)	    :: Qobs_data	! observed discharge
	character*10,dimension(:), allocatable, intent(in)  :: dates_data	! observed discharge
	real*8,dimension(17), intent(in)		            :: par_ini	    ! maximum parameter values
	real*8,dimension(17), intent(in)		            :: par_max	    ! maximum parameter values
	real*8,dimension(17), intent(in)		            :: par_min	    ! minimum parameter values
	logical,dimension(17), intent(in)		            :: optim	    ! minimum parameter values
	real*8,dimension(:), allocatable, intent(in)	    :: sumax_data 	! sumax data

	real*8,dimension(17)             		:: par_val	! minimum parameter values
	real*8,dimension(17)             		:: paramset	! set of parameters
	real*8,allocatable, dimension(:)		:: q		! modelled discharge
	real*8,allocatable, dimension(:)		:: qval		! validation discharge
    real*8,dimension(17) 				    :: r		    ! random number array
	real*8, dimension(28)            	    :: EC    	! evalutation criteria
	real*8, dimension(28)                   :: ECval    	! evalutation criteria
	real*8						            :: NSE		! Nash-Sutcliffe efficienct
	real*8						            :: LNSE		! log NSE
	real*8						            :: VE		! Volume error
	real*8						            :: FDC_NSE	! Nash of flow duration curve
	real*8						            :: NSE_val	! NSE for validation	
	real*8						            :: LNSE_val	! log NSE for validation
	real*8						            :: VE_val	! Volume error validation
	real*8					              	:: FDC_NSE_val	! validation Nash FDC
	integer						            :: length	! length of period
	integer						            :: val_length	! length of validation period
	real*8, allocatable, dimension(:,:)	    :: output	! output matrix
	real*8, allocatable, dimension(:,:)	    :: output_val	! output matrix validation
	integer						            :: nn 		! counter
	integer, allocatable, dimension(:)		:: seed		! initial random seed
	integer						            :: t		! counter
	integer						            :: count_feasibles !number of feasible sets
	integer						            :: count_solutions !number of final solutions
	integer						            :: j		! counter
	integer						            :: i		! counter
	integer						            :: k		! counter
	integer						            :: m		! counter
	integer						            :: jj		! counter
	integer						            :: ii		! counter
	real*8						            :: bound1	! pareto_bound
	real*8						            :: bound2	! pareto_bound
	real*8						            :: bound3	! pareto_bound
	real*8						            :: bound4	! pareto_bound
	real*8, dimension(:,:), allocatable		:: par_mat	! matrix with parameters
	real*8, dimension(:,:), allocatable		:: Q_mat	! matrix with modelled Q
	real*8, dimension(:,:), allocatable		:: qval_mat	! matrix with validated Q
	real*8, dimension(:,:), allocatable		:: outputval_mat! output validation		
	real*8, dimension(:,:), allocatable		:: obj_mat	! matrix with objective functions
	real*8, dimension(:,:), allocatable		:: objval_mat	! matrix with validation objective functions
	real*8, dimension(:,:), allocatable		:: pareto_obj	! matrix with pareto objective function values
	real*8,dimension(size(sumax_data))      :: sumax_data2! sumax data
    integer                                 :: fileunit	
	real*8, dimension(:,:), allocatable		:: final_obj	! matrix with final objective function values
	real*8, dimension(:,:), allocatable		:: final_par	! matrix with final parameters
	real*8, dimension(:,:,:), allocatable	:: final_out    ! matrix with final
	real*8, dimension(:,:), allocatable		:: final_states ! matrix with final states
	logical, dimension(:), allocatable		:: temp, temp1	! temperoray arrays
	character(256)                        	:: formData	! format of data
	character(256)                        	:: formDataObj	! format of data objectives
	character(256)                        	:: formDataPar	! format of data parameters
	character(256)                        	:: formHeader	! format of data headers
	character(256)                        	:: formHeaderPar! format of headers
	real*8, dimension(:,:), allocatable		:: incon_mat	! matrix with parameters
	real*8,dimension(4)        		        :: incon_val	! initial conditions states
	integer						            :: pareto_length! number of pareto solutions
	integer						            :: uFinalStates	! file unit Q validation
	integer						            :: uObj		! file unit all objectives
	integer						            :: uObjPar	! file unit pareto objectives
	integer						            :: uObjval 	! file unit objectives validation
	integer						            :: uParam	! file unit parameters
	integer						            :: uParamPar	! file unit pareto parameters
	integer						            :: uSig 	! file unit objectives validation
	integer						            :: uSigVal 	! file unit objectives validation
	integer						            :: uQ 		! file unit Q
	integer						            :: uQPar	! file unit Q
	integer						            :: uQval	! file unit Q validation
    real*8                                  :: fIt
    real                                    :: start, finish

call cpu_time(start)


uQ=10
uQval=12
uQPar=13
uParam=14
uParamPar=15
uObj=16
uObjval=17
uObjPar=18
uSig=19
uSigVal=20
uFinalStates=21


length=size(prec_data( iw_start:iv_end ))

allocate( q( length) )
allocate( Q_mat( length, int(Iterations)) )
allocate( par_mat(16, int(Iterations)) )
allocate( obj_mat(4, int(Iterations)) )
allocate( incon_mat(4, int(Iterations)) )

call init_random_seed()


!start monte carlo runs

count_feasibles = 0
fIt = 0.1
paramset=par_ini

 do nn=1, Iterations

	call random_number(r)

	!paramset=r*(par_max-par_min)+par_min	! random parameter set
forall(ii=1:size(paramset,1),optim(ii) .eqv. .TRUE.) paramset(ii)=r(ii)*(par_max(ii)-par_min(ii))+par_min(ii) 

        if(dyn_mode .eqv. .TRUE.) then 
           if(read_data .eqv. .TRUE.) then 
	      call model(paramset,incon, prec_data( iw_start:iv_end ) &
	      ,temp_data(iw_start:iv_end),etp_data(iw_start:iv_end), dem, cellsize, output,sumax_data)
            else
               call gen_sumax(paramset(15:17), paramset(4), ichange_start, ichange_end, prec_data( iw_start:iv_end ), sumax_data2 )

               call model(paramset,incon, prec_data( iw_start:iv_end ) &
	      ,temp_data(iw_start:iv_end),etp_data(iw_start:iv_end), dem, cellsize, output,sumax_data2)


            end if
         else

           call model(paramset,incon, prec_data( iw_start:iv_end ) &
	   ,temp_data(iw_start:iv_end),etp_data(iw_start:iv_end), dem, cellsize, output)
        end if
  

	q=output(1,:)

	!determine objectives
	 call NashEfficiency(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), NSE)
	 call LogNashEfficiency(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), LNSE)
	 call VolumeError(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), VE)
	 call FDCNashEfficiency(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), FDC_NSE)

	!keep parameter set when all objectives are bigger than 0, with maximum size of 5000

!	if( (NSE .gt. zero_dp) .and. (LNSE .gt. zero_dp) .and. (VE .gt. zero_dp) .and. (FDC_NSE .gt. zero_dp) ) then
	count_feasibles = count_feasibles + 1

	par_mat(:,count_feasibles)     = paramset		!every column is one run
	Q_mat(:, count_feasibles) = output(1,:)
	!output_mat(2:5,:, count_feasibles) = output(10:13,:)
        incon_mat(1:4, count_feasibles) = output(8:11,ic_end)
	obj_mat(:,count_feasibles)= (/NSE, LNSE, VE, FDC_NSE/)
        
!	end if


if(nn .eq. int(fIt*Iterations)) then
print *, "Iterations finished: ", nn
print *, "Feasible solutions : ", count_feasibles 
print *, ""
fIt = fIt + 0.1
end if



end do



print *, "Max NSE:"
print *, maxval(obj_mat(1,1:count_feasibles))
print *, "Max LNSE:"
print *, maxval(obj_mat(2,1:count_feasibles))
print *, "Max VE:"
print *, maxval(obj_mat(3,1:count_feasibles))
print *, "Max FDC_NSE:"
print *, maxval(obj_mat(4,1:count_feasibles))


print *, "Feasible parameter sets:"
print *, count_feasibles
print *, "Creating pareto front..."


!now create pareto front
allocate( temp ( count_feasibles))
allocate( temp1 ( count_feasibles))
allocate( pareto_obj(4, int(0.5*count_feasibles)) )

m=1
do i=1, count_feasibles
    temp(:)= .TRUE.
    do j=1,4
	temp1=(1-obj_mat (j,1:count_feasibles) .lt. 1-obj_mat (j,i)) 
	temp = ((temp1 .eqv. .TRUE.) .and. (temp .eqv. .TRUE.))
    end do
	if (any(temp) .eqv. .FALSE.) then
     	pareto_obj(1:4,m)    = obj_mat(1:4,i)
        m=m+1
    end if
end do

pareto_length=m-1

print *, "Number of pareto points:"
print *, pareto_length

! find space spanned by pareto front and save results to file

bound1=maxval(1-pareto_obj(1,1:pareto_length) )
bound2=maxval(1-pareto_obj(2,1:pareto_length) )
bound3=maxval(1-pareto_obj(3,1:pareto_length) )
bound4=maxval(1-pareto_obj(4,1:pareto_length) )

deallocate(pareto_obj)


open(uParam, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Parameters.txt")), status='unknown', action='write')
open(uQ, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Q.txt")), status='unknown', action='write')
open(uFinalStates, file=trim(adjustl(output_dir_cal)) // trim(adjustl("FinalStates.txt")), status='unknown', action='write')
open(uObj, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Objectives.txt")), status='unknown', action='write')

write(uParam,'(17(1X a8))')  "Meltfactor", "Tthresh", "Imax", "Sumax", "beta", &
                              "Kf", "Ks", "LP", "D", "Pmax", "Tlagf" ,"Tlags" ,"InfMax", "Kof","a", "b","sumax_min"

write(formHeader, *) '(2X a10, 2X a10, 2X a10, 2X a10)'
write(formDataObj, *) '(2X f10.3, 2X f10.3, 2X f10.3, 2X f10.3)'

write(uObj,formHeader) 'NSE','LNSE','VE','LNSE_FDC'


count_solutions=0
do i=1, count_feasibles
!   if( (1-obj_mat(1,i) .le. bound1) .and. &
!       (1-obj_mat(2,i) .le. bound2) .and. &
!       (1-obj_mat(3,i) .le. bound3) .and. &
!       (1-obj_mat(4,i) .le. bound4)) then
        count_solutions   = count_solutions+1
        write(uObj,formDataObj)  obj_mat(1:4,i)
        write(uParam,'(17(1X f10.2))')  par_mat(:,i)
        write(uQ,'(100000(1X f6.2))')  Q_mat(:,i)
        write(uFinalStates,*) incon_mat(:,i)
!   end if
end do

 close(uParam)
 close(uQ)
 close(uFinalStates)
 close(uObj)

print *, ""
print *, "Number of solutions in pareto space:"
print *, count_solutions


deallocate( Q_mat )
deallocate( par_mat    )
deallocate( obj_mat    )
deallocate( incon_mat  )
deallocate( temp       )
deallocate( temp1      )

 call cpu_time(finish)
 print '("Calibration finished in = ",f6.3," seconds.")',finish-start

 call cpu_time(start)

!*************************************************************************************
!post evalutation
!*************************************************************************************

print *, ""
print *, "Determine signatures  and validate..."
print *, ""

!allocate( EC(count_solutions, 26) )
!allocate( ECval(count_solutions, 26) )
val_length=size(prec_data( iv_start:iv_end ))

allocate( qval( val_length) )
allocate( objval_mat( 4, count_solutions) )

open(uParam, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Parameters.txt")), status='unknown', action='read')
open(uFinalStates, file=trim(adjustl(output_dir_cal)) // trim(adjustl("FinalStates.txt")), status='unknown', action='read')
open(uQ, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Q.txt")), status='unknown', action='read')

open(uQval, file=trim(adjustl(output_dir_val)) // trim(adjustl("Q.txt")), status='unknown', action='write')
open(uObjval, file=trim(adjustl(output_dir_val)) // trim(adjustl("Objectives.txt")), status='unknown', action='write')
 open(uSig, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Signatures.txt")), status='unknown', action='write')
 open(uSigVal, file=trim(adjustl(output_dir_val)) // trim(adjustl("Signatures.txt")), status='unknown', action='write')

write(formHeader, *) '(2X a10, 2X a10, 2X a10, 2X a10)'
write(formDataObj, *) '(2X f10.3, 2X f10.3, 2X f10.3, 2X f10.3)'
write(uObjval,formHeader) 'NSE','LNSE','VE','LNSE_FDC'
write(uSig,'(28(A12))') 'Q_MA', 'AC', 'AC_low', 'AC_high', 'RLD', 'DLD', 'Q5', 'Q50', 'Q95', 'Q5_low', 'Q50_low', 'Q95_low', &
                        'Q5_high', 'Q50_high', 'Q95_high', 'Peaks', 'Peaks_low', 'Peaks_high', 'Qpeak10','Qpeak50', &
                        'Qlow_peak10', 'Qlow_peak50','Qhigh_peak10', 'Qhigh_peak50', 'SFDC', 'LFR', 'FDC_serie', 'AC_serie'
write(uSigVal,'(28(A12))') 'Q_MA', 'AC', 'AC_low', 'AC_high', 'RLD', 'DLD', 'Q5', 'Q50', 'Q95', 'Q5_low', 'Q50_low', 'Q95_low',&
                           'Q5_high','Q50_high', 'Q95_high', 'Peaks', 'Peaks_low', 'Peaks_high', 'Qpeak10','Qpeak50', & 
                           'Qlow_peak10', 'Qlow_peak50', 'Qhigh_peak10', 'Qhigh_peak50', 'SFDC', 'LFR','FDC_serie', 'AC_serie'


read(uParam,*) !skip header


fIt = 0.1



   


do i=1, count_solutions



if(i .eq. int(fIt*count_solutions)) then
print *, "Post-evaluations finished : ", i
print *, "Total solutions           : ", count_solutions 
print *, ""
fIt = fIt + 0.1
end if



   read(uParam,*) par_val
   read(uFinalStates,*) incon_val
   read(uQ,*) q

   !evaluate signatures for calibration period

   call eval_signatures(q(iw_end+1:ic_end), Qobs_data(iw_end+1:ic_end), prec_data( iw_end+1:ic_end ), dates_data(iw_end+1:ic_end), &
                     .FALSE., EC )
   !validate



   call model(par_val,incon_val, prec_data( iv_start:iv_end ) &
	,temp_data(iv_start:iv_end),etp_data(iv_start:iv_end), dem, cellsize, output_val)

	qval=output_val(1,:)

	 call NashEfficiency(Qobs_data(iv_start:iv_end),qval, NSE_val)
	 call LogNashEfficiency(Qobs_data(iv_start:iv_end),qval, LNSE_val)
	 call VolumeError(Qobs_data(iv_start:iv_end),qval, VE_val)
	 call FDCNashEfficiency(Qobs_data(iv_start:iv_end),qval, FDC_NSE_val)

	
        objval_mat(:,i)= (/NSE_val, LNSE_val, VE_val, FDC_NSE_val/)
	!outputval_mat(:,i)=qval

   !evaluate signatures validation

         call eval_signatures(qval, Qobs_data(iv_start:iv_end), prec_data( iv_start:iv_end ), dates_data(iv_start:iv_end), &
                     .FALSE., ECval )

  	write(uQval,'(100000(1X f6.2))')  qval
  	write(uObjval,formDataObj)  (/NSE_val, LNSE_val, VE_val, FDC_NSE_val/)
  	write(uSig,'(28(1X f8.3))')  EC
  	write(uSigVal,'(28(1X f8.3))')  ECval

end do



 close(uParam)
 close(uQ)
 close(uFinalStates)

 close(uQval)
 close(uObjval)
 close(uSig)
 close(uSigval)

deallocate( q          )
             
 call cpu_time(finish)
 print '("Post-evaluation finished in = ",f6.3," seconds.")',finish-start




print *, "Optimization finished!"

end subroutine



end module mo_multi_obj
