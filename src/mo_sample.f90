MODULE mo_sample



IMPLICIT NONE

CONTAINS

	subroutine sample(par_ini, par_max,par_min,incon, prec_data,temp_data,etp_data,&
                             dem, cellsize,Qobs_data, dates_data, optim, sumax_data, imax_data)
		
    USE mo_model
    USE mo_readdata
    USE mo_objectives
    USE mo_eval_signatures
    USE mo_init_random_seed

	real*8,dimension(4), intent(in)		                :: incon	    ! initial conditions states
	real*8,dimension(:), allocatable, intent(in)	    :: etp_data	    ! evaporation data
	real*8,dimension(:), allocatable, intent(in)	    :: prec_data 	! precipitation data
	real*8,dimension(:), allocatable, intent(in)	    :: temp_data 	! temperature data
    real*8, dimension(:,:), allocatable, intent(in)     :: dem          ! elevation
    real*8,intent(in)                                   :: cellsize     ! size of grid dem
	real*8,dimension(:), allocatable, intent(in)	    :: Qobs_data	! observed discharge
	character*10,dimension(:), allocatable, intent(in)	:: dates_data	! observed discharge
    real*8,dimension(14) 				                :: r		    ! random number array
	real*8,dimension(14), intent(in)		            :: par_ini	    ! maximum parameter values
	real*8,dimension(14), intent(in)		            :: par_max	    ! maximum parameter values
	real*8,dimension(14), intent(in)		            :: par_min	    ! minimum parameter values
	logical,dimension(14), intent(in)		            :: optim	    ! parameters to optimize
	real*8,dimension(:), allocatable, intent(in)	    :: sumax_data 	! sumax data
	real*8,dimension(:), allocatable,  intent(in)	    :: Imax_data 	! interception cap. data


	real*8, allocatable, dimension(:,:)	                :: output       ! output matrix
	real*8, allocatable, dimension(:,:)	                :: output_val   ! output matrix validation
	integer						                        :: nn           ! counter
	integer, allocatable, dimension(:)		            :: seed		    ! initial random seed
	integer						                        :: count_feasibles !number of feasible sets
	integer						                        :: count_solutions !number of final solutions
	real*8, dimension(:),   allocatable		            :: Ebest	    ! matrix with best evaporation
	real*8, dimension(28)            	                :: EC    	    ! evalutation criteria
	real*8, dimension(28)                               :: ECval        ! evalutation criteria
	real*8						                        :: FDC_NSE	    ! Nash of flow duration curve
	real*8						                        :: FDC_NSE_val  ! validation Nash FDC
    integer                                             :: fileunit	    ! number of file
	real*8, dimension(:,:), allocatable		            :: final_obj	! matrix with final objective function values
	real*8, dimension(:,:), allocatable		            :: final_par	! matrix with final parameters
	real*8, dimension(:,:,:), allocatable		        :: final_out    ! matrix with final
	real*8, dimension(:,:), allocatable		            :: final_states ! matrix with final states
    real                                                :: finish       ! end time
	logical, dimension(:), allocatable		            :: temp, temp1	! temperoray arrays
	character(256)                        		        :: formData	    ! format of data
	character(256)                        		        :: formDataObj	! format of data objectives
	character(256)                        		        :: formDataPar	! format of data parameters
	character(256)                        		        :: formHeader	! format of data headers
	character(256)                        		        :: formHeaderPar! format of headers
	integer						                        :: i		    ! counter
	real*8,dimension(size(sumax_data))                  :: Imax_data2 ! interception cap. data
	integer						                        :: ii		    ! counter
	real*8,dimension(4)        		                    :: incon_val	! initial conditions states
	integer						                        :: j		    ! counter
	integer						                        :: jj		    ! counter
	integer						                        :: k		    ! counter
	real*8						                        :: KGE	        ! Kling-gupta efficiency
	real*8						                        :: KGEmax	    ! Highest Kling-gupta efficiency
	real*8						                        :: KGE_val	  ! Kling-gupta efficiency validation
	integer						                        :: length	  ! length of period
	real*8						                        :: LKGE		  ! log KGE
	real*8						                        :: LKGE_corr  ! log KGE corrected
	real*8						                        :: LKGE_val	  ! log KGE validation
	real*8						                        :: LKGE_valcorr ! log KGE corr. validation
	real*8						                        :: LNSE		  ! log NSE
	real*8						                        :: LNSE_corr  ! log NSE corrected
	real*8						                        :: LNSE_val	  ! log NSE for validation
	real*8						                        :: LNSE_valcorr  ! log NSE for validation
	integer						                        :: m		    ! counter
	real*8						                        :: NSE          ! Nash-Sutcliffe efficiency
	real*8						                        :: NSE_val	    ! NSE for validation	
	real*8, dimension(:,:), allocatable		            :: outputval_mat! output validation		
	real*8, dimension(:,:), allocatable		            :: obj_mat	    ! matrix with objective functions
	real*8, dimension(:,:), allocatable		            :: objval_mat	! matrix with validation objective functions
	real*8,dimension(14)             		            :: paramset	    ! set of parameters
	real*8,dimension(14)             		            :: parbest      ! best set of parameters
	integer						                        :: pareto_length! number of pareto solutions
	real*8, dimension(:,:), allocatable		            :: pareto_obj	! matrix with pareto objective function values
	real*8,dimension(14)             		            :: par_val	    ! minimum parameter values
	real*8, dimension(:,:), allocatable		            :: par_mat	    ! matrix with parameters
	real*8,allocatable, dimension(:)		            :: q		    ! modelled discharge
	real*8, dimension(:),   allocatable		            :: Qbest	    ! matrix with modelled Q
	real*8,allocatable, dimension(:)		            :: qval		    ! validation discharge
	real*8, dimension(:,:), allocatable		            :: qval_mat	    ! matrix with validated Q

    real*8                                              :: r2           ! random number
	real*8						                        :: snow_obj   ! snow objective
    real                                                :: start        ! starttime
	real*8,dimension(size(sumax_data))                  :: sumax_data2! sumax data
	integer						                        :: t	      ! counter
	integer						                        :: uEbest_cal ! file unit E calibration
	integer						                        :: uEbest_val ! file unit E validation
	integer						                        :: uFinalStates	! file unit Q validation
	integer						                        :: uObj		! file unit all objectives
	integer						                        :: uObjPar	! file unit pareto objectives
	integer						                        :: uObjval 	! file unit objectives validation
	integer						                        :: uParam	! file unit parameters
	integer						                        :: uParamPar	! file unit pareto parameters
	integer						                        :: uQ 		! file unit Q
	integer						                        :: uQPar	! file unit Q
	integer						                        :: uQval	! file unit Q validation
	integer						                        :: uQbest_cal ! file unit Q calibration
	integer						                        :: uQbest_val ! file unit Q validation
	integer						                        :: val_length ! length of validation period
	real*8						                        :: VE		  ! Volume error
	real*8						                        :: VE_val	  ! Volume error validation
    real*8                                              :: fIt


   call cpu_time(start)


   uQ=10
   uQval=12
   uQPar=13
   uParam=14
   uParamPar=15
   uObj=16
   uObjval=17
   uObjPar=18
   uFinalStates=21
   uQbest_cal=22
   uQbest_val=23
   uEbest_cal=24
   uEbest_val=25

   !length of full timeseries
   length=size(prec_data( iw_start:iv_end ))

   !allocate timeseries
   allocate( q( length) )
   allocate( Qbest( length ) )
   allocate( Ebest( length ) )
   allocate( par_mat(14, int(Iterations)) )
   allocate( obj_mat(7, int(Iterations)) )
   allocate( objval_mat( 7, int(Iterations)) )

   !initialize random seed
   call init_random_seed()


   !start monte carlo runs
   count_feasibles = 0
   fIt = 0.1
   paramset=par_ini
   KGEmax = -9999


   do nn=1, Iterations

	call random_number(r)

	! random parameter set based on random number and max min values
    forall(ii=1:size(paramset,1),optim(ii) .eqv. .TRUE.) paramset(ii)=r(ii)*(par_max(ii)-par_min(ii))+par_min(ii) 


        if(dyn_mode .eqv. .TRUE.) then 

        paramset(4)=r2*(par_max(4)-paramset(14))+paramset(14)

           !reading in root zone data
           if(read_data .eqv. .TRUE.) then 
              !make a timeserie of Imax
              call gen_sumax(paramset(12:14), paramset(3), ichange_start, ichange_end, prec_data( iw_start:iv_end ), imax_data2 )

	       call model(paramset,incon, prec_data( iw_start:iv_end ) &
	      ,temp_data(iw_start:iv_end),etp_data(iw_start:iv_end), dem, cellsize, output,sumax_data, imax_data2)
           
           ! else make a timeserie of sumax
           else               
               call gen_sumax(paramset(12:14), paramset(4), ichange_start, ichange_end, prec_data( iw_start:iv_end ), sumax_data2 )

               !call the model
               call model(paramset,incon, prec_data( iw_start:iv_end ) &
	                 ,temp_data(iw_start:iv_end),etp_data(iw_start:iv_end), dem, cellsize, output,sumax_data2)
           end if
         else

           !call the model without dynamic root zone values
           call model(paramset,incon, prec_data( iw_start:iv_end ) &
	   ,temp_data(iw_start:iv_end),etp_data(iw_start:iv_end), dem, cellsize, output)
        end if 

	q=output(1,:)

	!determine objectives
	 call NashEfficiency(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), NSE)
	 call LogNashEfficiency(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), LNSE)
	 call LogNashEfficiency_corr(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), LNSE_corr)
	 call VolumeError(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), VE)
	 call logKlingGupta(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), LKGE)
	 call KlingGupta(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), KGE)
	 call logKlingGupta_corr(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), LKGE_corr)

     if(KGE .gt. KGEmax) then
        Qbest = output(1,:)
        Ebest = output(5,:) + output(6,:)
        parbest = paramset
     end if

	par_mat(:, nn)     = paramset		!every column is one run


	obj_mat(:, nn)     = (/NSE, LNSE, VE, KGE, LKGE, LKGE_corr, LNSE_corr/)
        
	!determine objectives validation period
	 call NashEfficiency(Qobs_data(iv_start:iv_end),q(iv_start:iv_end), NSE_val)
	 call LogNashEfficiency(Qobs_data(iv_start:iv_end),q(iv_start:iv_end), LNSE_val)
	 call LogNashEfficiency_corr(Qobs_data(iv_start:iv_end),q(iv_start:iv_end), LNSE_valcorr)
	 call VolumeError(Qobs_data(iv_start:iv_end),q(iv_start:iv_end), VE_val)
	 call KlingGupta(Qobs_data(iv_start:iv_end),q(iv_start:iv_end), KGE_val)
     call logKlingGupta(Qobs_data(iv_start:iv_end),q(iv_start:iv_end), LKGE_val)
	 call logKlingGupta_corr(Qobs_data(iv_start:iv_end),q(iv_start:iv_end), LKGE_valcorr)
	
     objval_mat(:,nn)= (/NSE_val, LNSE_val, VE_val, KGE_val, LKGE_val, LKGE_valcorr, LNSE_valcorr/)


if(nn .eq. int(fIt*Iterations)) then
print *, "Iterations finished: ", nn
print *, "Total iterations   : ", Iterations
print *, ""
fIt = fIt + 0.1
end if


end do


print *, "Max NSE:"
print *, maxval(obj_mat(1,1:Iterations))
print *, "Max LNSE:"
print *, maxval(obj_mat(2,1:Iterations))
print *, "Max VE:"
print *, maxval(obj_mat(3,1:Iterations))
print *, "Max KGE:"
print *, maxval(obj_mat(4,1:Iterations))
print *, "Max LKGE:"
print *, maxval(obj_mat(5,1:Iterations))
print *, "Max LKGE corrected:"
print *, maxval(obj_mat(6,1:Iterations))
print *, "Max LNSE corrected:"
print *, maxval(obj_mat(7,1:Iterations))


!write calibration data to files
open(uParam, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Parameters.txt")), status='unknown', action='write')
open(uQbest_cal, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Qbest.txt")), status='unknown', action='write')
open(uEbest_cal, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Ebest.txt")), status='unknown', action='write')
open(uObj, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Objectives.txt")), status='unknown', action='write')
open(uObjval, file=trim(adjustl(output_dir_val)) // trim(adjustl("Objectives.txt")), status='unknown', action='write')

write(uParam,'(14(1X a8))')  "Meltfactor", "Tthresh", "Imax", "Sumax", "beta", &
                              "Kf", "Ks", "LP", "D", "Pmax","alpha","a", "b", "sumax_min"

write(formHeader, *) '(2X a10, 2X a10, 2X a10, 2X a10, 2X a10,2X a10, 2X a10, 2X a10)'
write(formDataObj, *) '(2X f10.3, 2X f10.3, 2X f10.3, 2X f10.3, 2X f10.6,2X f10.6, 2X f10.3, 2X f10.3)'

write(uObj,formHeader) 'NSE','LNSE','VE','KGE', 'LKGE', 'LKGE_corr', 'LNSE_corr', 'snow_obj'

do i=1, size(Qbest)
   write(uQbest_cal,*)  Qbest(i)
   write(uEbest_cal,*)  Ebest(i)
end do

do i=1, Iterations

        count_solutions   = count_solutions+1
        write(uObj,formDataObj)  obj_mat(:,i)
  	    write(uObjval,formDataObj)  objval_mat(:,i)
        write(uParam,'(14(1X f10.5))')  par_mat(:,i)

end do

 close(uParam)
 close(uQbest_cal)
 close(uEbest_cal)
 close(uObj)



deallocate( par_mat    )
deallocate( obj_mat    )
deallocate( objval_mat    )


 call cpu_time(finish)
 print '("Calibration finished in = ",f12.3," seconds.")',finish-start

 call cpu_time(start)




print *, "Optimization finished!"

end subroutine




end module mo_sample
