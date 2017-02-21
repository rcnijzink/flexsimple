MODULE mo_sample



IMPLICIT NONE

CONTAINS

	subroutine sample(par_ini, par_max,par_min,incon, prec_data,temp_data,etp_data,&
                             dem, cellsize,Qobs_data, dates_data, optim, sumax_data, imax_data, snow_data)
		
USE mo_model
USE mo_readdata
USE mo_objectives
USE mo_eval_signatures
USE mo_init_random_seed

	real*8,dimension(4), intent(in)		        :: incon	! initial conditions states
	real*8,dimension(:), allocatable, intent(in)	:: etp_data	! evaporation data
	real*8,dimension(:), allocatable, intent(in)	:: prec_data 	! precipitation data
	real*8,dimension(:), allocatable, intent(in)	:: temp_data 	! temperature data
        real*8, dimension(:,:), allocatable, intent(in)     :: dem
        real*8,intent(in)                       :: cellsize
	real*8,dimension(:), allocatable, intent(in)	:: Qobs_data	! observed discharge
	character*10,dimension(:), allocatable, intent(in)	:: dates_data	! observed discharge
        real*8,dimension(14) 				:: r		! random number array
	real*8,dimension(14), intent(in)		:: par_ini	! maximum parameter values
	real*8,dimension(14), intent(in)		:: par_max	! maximum parameter values
	real*8,dimension(14), intent(in)		:: par_min	! minimum parameter values
	logical,dimension(14), intent(in)		:: optim	! minimum parameter values
	real*8,dimension(:), allocatable, intent(in)	:: sumax_data 	! sumax data
	real*8,dimension(:), allocatable, intent(in)	:: Imax_data 	! sumax data
	real*8,dimension(:), allocatable, intent(in)	:: snow_data 	! sumax data
	real*8,dimension(size(sumax_data))              :: sumax_data2 	! sumax data
	real*8,dimension(size(sumax_data))              :: Imax_data2 	! sumax data

	real*8,dimension(14)             		:: par_val	! minimum parameter values
	real*8,dimension(14)             		:: paramset	! set of parameters
	real*8,allocatable, dimension(:)		:: q		! modelled discharge
	real*8,allocatable, dimension(:)		:: qval		! validation discharge
	real*8, dimension(28)            	        :: EC    	! evalutation criteria
	real*8, dimension(28)                           :: ECval    	! evalutation criteria
	real*8						:: NSE		! Nash-Sutcliffe efficienct
	real*8						:: snow_obj	! snow objective
	real*8						:: LNSE		! log NSE
	real*8						:: LNSE_corr		! log NSE
	real*8						:: KGE		! log NSE
	real*8						:: KGE_val	! log NSE
	real*8						:: LKGE		! log NSE
	real*8						:: LKGE_corr		! log NSE
	real*8						:: LKGE_val	! log NSE
	real*8						:: VE		! Volume error
	real*8						:: FDC_NSE	! Nash of flow duration curve
	real*8						:: NSE_val	! NSE for validation	
	real*8						:: LNSE_val	! log NSE for validation
	real*8						:: VE_val	! Volume error validation
	real*8						:: FDC_NSE_val	! validation Nash FDC
	integer						:: length	! length of period
	integer						:: val_length	! length of validation period
	real*8, allocatable, dimension(:,:)	        :: output	! output matrix
	real*8, allocatable, dimension(:,:)	        :: output_val	! output matrix validation
	integer						:: nn 		! counter
	integer, allocatable, dimension(:)		:: seed		! initial random seed
	integer						:: t		! counter
	integer						:: count_feasibles !number of feasible sets
	integer						:: count_solutions !number of final solutions
	integer						:: j		! counter
	integer						:: i		! counter
	integer						:: k		! counter
	integer						:: m		! counter
	integer						:: jj		! counter
	integer						:: ii		! counter
	real*8						:: bound1	! pareto_bound
	real*8						:: bound2	! pareto_bound
	real*8						:: bound3	! pareto_bound
	real*8						:: bound4	! pareto_bound
	real*8, dimension(:,:), allocatable		:: par_mat	! matrix with parameters
	real*8, dimension(:,:), allocatable		:: Q_mat	! matrix with modelled Q
	real*8, dimension(:,:), allocatable		:: qval_mat	! matrix with validated Q
	real*8, dimension(:,:), allocatable		:: outputval_mat! output validation		
	real*8, dimension(:,:), allocatable		:: obj_mat	! matrix with objective functions
	real*8, dimension(:,:), allocatable		:: objval_mat	! matrix with validation objective functions
	real*8, dimension(:,:), allocatable		:: pareto_obj	! matrix with pareto objective function values
!	real*8, dimension(:,:), allocatable		:: pareto_par	! matrix with pareto parameters
!	real*8, dimension(:,:,:), allocatable		:: pareto_output! matrix with pareto output
        real*8                                          :: r2
        integer                                         :: fileunit	
	real*8, dimension(:,:), allocatable		:: final_obj	! matrix with final objective function values
	real*8, dimension(:,:), allocatable		:: final_par	! matrix with final parameters
	real*8, dimension(:,:,:), allocatable		:: final_out    ! matrix with final
	real*8, dimension(:,:), allocatable		:: final_states ! matrix with final states
	logical, dimension(:), allocatable		:: temp, temp1	! temperoray arrays
	character(256)                        		:: formData	! format of data
	character(256)                        		:: formDataObj	! format of data objectives
	character(256)                        		:: formDataPar	! format of data parameters
	character(256)                        		:: formHeader	! format of data headers
	character(256)                        		:: formHeaderPar! format of headers
	real*8, dimension(:,:), allocatable		:: incon_mat	! matrix with parameters
	real*8,dimension(4)        		        :: incon_val	! initial conditions states
	integer						:: pareto_length! number of pareto solutions
	integer						:: uFinalStates	! file unit Q validation
	integer						:: uObj		! file unit all objectives
	integer						:: uObjPar	! file unit pareto objectives
	integer						:: uObjval 	! file unit objectives validation
	integer						:: uParam	! file unit parameters
	integer						:: uParamPar	! file unit pareto parameters
	integer						:: uSig 	! file unit objectives validation
	integer						:: uSigVal 	! file unit objectives validation
	integer						:: uQ 		! file unit Q
	integer						:: uQPar	! file unit Q
	integer						:: uQval	! file unit Q validation
        real*8                                          :: fIt

real :: start, finish

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
allocate( par_mat(14, int(Iterations)) )
allocate( obj_mat(8, int(Iterations)) )
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

	!call random_number(r2)

        !if(optim(7) .eqv. .TRUE.) then
        !   paramset(7)=r2*(par_max(7)-paramset(6))+paramset(6)
        !   call random_number(r2)
        !end if

        if(dyn_mode .eqv. .TRUE.) then 

        paramset(4)=r2*(par_max(4)-paramset(14))+paramset(14)

           if(read_data .eqv. .TRUE.) then 
              !make a timeserie of Imax
              call gen_sumax(paramset(12:14), paramset(3), ichange_start, ichange_end, prec_data( iw_start:iv_end ), imax_data2 )
              !make a time serie of the meltfactor
              !call gen_sumax(paramset(12:14), paramset(1), ichange_start, ichange_end, prec_data( iw_start:iv_end ), melta_data )

	      call model(paramset,incon, prec_data( iw_start:iv_end ) &
	      ,temp_data(iw_start:iv_end),etp_data(iw_start:iv_end), dem, cellsize, output,sumax_data, imax_data2)
            else
               !make a timeserie of sumax
               call gen_sumax(paramset(12:14), paramset(4), ichange_start, ichange_end, prec_data( iw_start:iv_end ), sumax_data2 )
               !call the model
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
	 call LogNashEfficiency_corr(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), LNSE_corr)
	 call VolumeError(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), VE)
	 call logKlingGupta(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), LKGE)
	 call KlingGupta(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), KGE)
	 call logKlingGupta_corr(Qobs_data(iw_end+1:ic_end),q(iw_end+1:ic_end), LKGE_corr)
	 call snow_objective(snow_data(iw_end+1:ic_end),output(12,iw_end+1:ic_end), snow_obj)


	par_mat(:, nn)     = paramset		!every column is one run
	Q_mat(:, nn)       = output(1,:)
	!output_mat(2:5,:, count_feasibles) = output(10:13,:)
        incon_mat(1:4, nn) = output(8:11,ic_end)
	obj_mat(:, nn)     = (/NSE, LNSE, VE, KGE, LKGE, LKGE_corr, LNSE_corr, snow_obj/)
        

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
print *, "Max snow objective:"
print *, maxval(obj_mat(8,1:Iterations))


open(uParam, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Parameters.txt")), status='unknown', action='write')
open(uQ, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Q.txt")), status='unknown', action='write')
open(uFinalStates, file=trim(adjustl(output_dir_cal)) // trim(adjustl("FinalStates.txt")), status='unknown', action='write')
open(uObj, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Objectives.txt")), status='unknown', action='write')

write(uParam,'(14(1X a8))'),  "Meltfactor", "Tthresh", "Imax", "Sumax", "beta", &
                              "Kf", "Ks", "LP", "D", "Pmax","alpha","a", "b", "sumax_min"

write(formHeader, *) '(2X a10, 2X a10, 2X a10, 2X a10, 2X a10,2X a10, 2X a10, 2X a10)'
write(formDataObj, *) '(2X f10.3, 2X f10.3, 2X f10.3, 2X f10.3, 2X f10.6,2X f10.6, 2X f10.3, 2X f10.3)'

write(uObj,formHeader), 'NSE','LNSE','VE','KGE', 'LKGE', 'LKGE_corr', 'LNSE_corr', 'snow_obj'


do i=1, Iterations

        count_solutions   = count_solutions+1
        write(uObj,formDataObj),  obj_mat(1:8,i)
        write(uParam,'(14(1X f10.5))'),  par_mat(:,i)
        write(uQ,'(100000(1X f6.2))'),  Q_mat(:,i)
        write(uFinalStates,*), incon_mat(:,i)

end do

 close(uParam)
 close(uQ)
 close(uFinalStates)
 close(uObj)


deallocate( Q_mat )
deallocate( par_mat    )
deallocate( obj_mat    )
deallocate( incon_mat  )


 call cpu_time(finish)
 print '("Calibration finished in = ",f6.3," seconds.")',finish-start

 call cpu_time(start)

!*************************************************************************************
!post evalutation
!*************************************************************************************

print *, ""
print *, "Determine signatures  and validate..."
print *, ""


val_length=size(prec_data( iv_start:iv_end ))

allocate( qval( val_length) )
allocate( objval_mat( 5, Iterations) )

!files to read
 open(uParam, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Parameters.txt")), status='unknown', action='read')
 open(uFinalStates, file=trim(adjustl(output_dir_cal)) // trim(adjustl("FinalStates.txt")), status='unknown', action='read')
 open(uQ, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Q.txt")), status='unknown', action='read')

!files to write
 open(uQval, file=trim(adjustl(output_dir_val)) // trim(adjustl("Q.txt")), status='unknown', action='write')
 open(uObjval, file=trim(adjustl(output_dir_val)) // trim(adjustl("Objectives.txt")), status='unknown', action='write')
 open(uSig, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Signatures.txt")), status='unknown', action='write')
 open(uSigVal, file=trim(adjustl(output_dir_val)) // trim(adjustl("Signatures.txt")), status='unknown', action='write')

!headers
write(formHeader, *) '(2X a10, 2X a10, 2X a10, 2X a10, 2X a10)'
write(formDataObj, *) '(2X f10.3, 2X f10.3, 2X f10.3, 2X f10.3, 2X f10.3)'
write(uObjval,formHeader), 'NSE','LNSE','VE','KGE' , 'LKGE'
write(uSig,'(28(A12))'), 'Q_MA', 'AC', 'AC_low', 'AC_high', 'RLD', 'DLD', 'Q5', 'Q50', 'Q95', 'Q5_low', 'Q50_low', 'Q95_low', &
                        'Q5_high', 'Q50_high', 'Q95_high', 'Peaks', 'Peaks_low', 'Peaks_high', 'Qpeak10','Qpeak50', &
                        'Qlow_peak10', 'Qlow_peak50','Qhigh_peak10', 'Qhigh_peak50', 'SFDC', 'LFR', 'FDC_serie', 'AC_serie'
write(uSigVal,'(28(A12))'), 'Q_MA', 'AC', 'AC_low', 'AC_high', 'RLD', 'DLD', 'Q5', 'Q50', 'Q95', 'Q5_low', 'Q50_low', 'Q95_low',&
                           'Q5_high','Q50_high', 'Q95_high', 'Peaks', 'Peaks_low', 'Peaks_high', 'Qpeak10','Qpeak50', & 
                           'Qlow_peak10', 'Qlow_peak50', 'Qhigh_peak10', 'Qhigh_peak50', 'SFDC', 'LFR','FDC_serie', 'AC_serie'


read(uParam,*) !skip

fIt = 0.1

do i=1, Iterations

if(i .eq. int(fIt*Iterations)) then
   print *, "Post-evaluations finished : ", i
   print *, "Total solutions           : ", Iterations
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
	 call KlingGupta(Qobs_data(iv_start:iv_end),qval, KGE_val)
         call logKlingGupta(Qobs_data(iv_start:iv_end),qval, LKGE_val)

	
        objval_mat(:,i)= (/NSE_val, LNSE_val, VE_val, KGE_val, LKGE_val/)
	!outputval_mat(:,i)=qval

   !evaluate signatures validation

         call eval_signatures(qval, Qobs_data(iv_start:iv_end), prec_data( iv_start:iv_end ), dates_data(iv_start:iv_end), &
                     .FALSE., ECval )

  	write(uQval,'(100000(1X f6.2))'),  qval
  	write(uObjval,formDataObj),  (/NSE_val, LNSE_val, VE_val, KGE_val, LKGE_val/)
  	write(uSig,'(28(1X f8.3))'),  EC
  	write(uSigVal,'(28(1X f8.3))'),  ECval

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




end module mo_sample
