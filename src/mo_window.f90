MODULE mo_moving_window


IMPLICIT NONE

CONTAINS

	subroutine moving_window(par_ini, par_max,par_min,incon, prec_data,temp_data,etp_data,&
                             dem, cellsize,Qobs_data, dates_data, optim)
		
USE mo_model
USE mo_readdata
USE mo_objectives
USE mo_eval_signatures
USE mo_signatures
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
	real*8,dimension(14)             		:: par_val	! minimum parameter values
	real*8,dimension(14)             		:: paramset	! set of parameters
	real*8,allocatable, dimension(:)		:: q		! modelled discharge
	real*8,allocatable, dimension(:)		:: qval		! validation discharge
	real*8, dimension(26)            	        :: EC    	! evalutation criteria
	real*8, dimension(26) 	                        :: ECval    	! evalutation criteria
	real*8						:: NSE		! Nash-Sutcliffe efficienct
	real*8						:: LNSE		! log NSE
	real*8						:: VE		! Volume error
	real*8						:: FDC_NSE	! Nash of flow duration curve
	real*8						:: NSE_val	! NSE for validation	
	real*8						:: LNSE_val	! log NSE for validation
	real*8						:: VE_val	! Volume error validation
	real*8						:: FDC_NSE_val	! validation Nash FDC
	integer						:: length	! length of period
	integer						:: val_length	! length of validation period
	real*8, allocatable, dimension(:,:)	        :: output	! output matrix
	real*8, allocatable, dimension(:,:)	        :: output_ini	! output matrix
	real*8, allocatable, dimension(:,:)	        :: output_val	! output matrix validation
	integer						:: nn 		! counter
	integer, allocatable, dimension(:)		:: seed		! initial random seed
	integer						:: t		! counter
	integer, allocatable, dimension(:)		:: count_feasibles !number of feasible sets
	integer						:: count_solutions !number of final solutions
	integer						:: j		! counter
	integer						:: i		! counter
	integer						:: k		! counter
	integer						:: m		! counter
	integer						:: jj		! counter
	integer						:: ii		! counter
        integer                                         :: iWindow      ! window counter
	real*8						:: bound1	! pareto_bound
	real*8						:: bound2	! pareto_bound
	real*8						:: bound3	! pareto_bound
	real*8						:: bound4	! pareto_bound
        integer                                         :: end_ind
        integer                                         :: endyear
        integer                                         :: extra
	real*8,dimension(:), allocatable        	:: prec_in 	! precipitation data
	real*8,dimension(:), allocatable        	:: etp_in 	! precipitation data
	real*8,dimension(:), allocatable        	:: temp_in 	! precipitation data
	real*8, dimension(:,:), allocatable		:: par_mat	! matrix with parameters
	real*8, dimension(:,:), allocatable		:: Q_mat	! matrix with modelled Q
	real*8, dimension(:,:), allocatable		:: qval_mat	! matrix with validated Q
	real*8, dimension(:,:), allocatable		:: outputval_mat! output validation		
	real*8, dimension(:,:), allocatable		:: obj_mat	! matrix with objective functions
	real*8, dimension(:,:), allocatable		:: pareto_obj	! matrix with pareto objective function values
	real*8, dimension(:,:), allocatable		:: pareto_par	! matrix with pareto objective function values
	real*8, dimension(:,:), allocatable		:: pareto_q	! matrix with pareto objective function values
        integer                                         :: start_ind
        integer                                         :: startyear
        character(10)                                    :: temp_date
        character(10)                                    :: temp_date2
        integer                                         :: fileunit	
	real*8, dimension(:,:), allocatable		:: final_obj	! matrix with final objective function values
	real*8, dimension(:,:), allocatable		:: final_par	! matrix with final parameters
	real*8, dimension(:,:,:), allocatable		:: final_out    ! matrix with final
	real*8, dimension(:,:), allocatable		:: final_states ! matrix with final states
	logical, dimension(:), allocatable		:: temp, temp1	! temperoray arrays
	character(256)                        		:: formData	! format of data
	character(256)                        		:: formDataObj	! format of data objectives
	character(256)                        		:: formDataSig	! format of data objectives
	character(256)                        		:: formDataPar	! format of data parameters
	character(256)                        		:: formHeader	! format of data headers
	character(256)                        		:: formHeaderSig! format of signature headers
	character(256)                        		:: formHeaderPar! format of headers
	real*8, dimension(:,:,:), allocatable		:: incon_mat	! matrix with parameters
	real*8,dimension(4)        		        :: incon_val	! initial conditions states
        integer                                         :: nWindow
	integer						:: pareto_length! number of pareto solutions
	integer						:: uFinalStates	! file unit Q validation
	integer						:: uObj		! file unit all objectives
	integer						:: uParam	! file unit parameters
	integer						:: uQ 		! file unit Q
	integer						:: uSig 		! file unit signatures
	integer						:: uEa 		! file unit signatures
	integer						:: uPeff 		! file unit signatures
	integer						:: uEi 		! file unit signatures
	integer						:: uSu 		! file unit signatures
        real*8                                          :: fIt
        character(10)                                   :: window_string


        real*8                                          :: RCmod
        real*8                                          :: var_p
        real*8                                          :: var_q
        real*8                                          :: var_qmod
        real*8                                          :: var_mod
        real*8                                          :: var
        real*8,dimension(:,:), allocatable              :: FDC
        real*8                                          :: KGE
        real*8                                          :: LKGE
        real*8                                          :: Q5
        real*8                                          :: Q95
        real*8                                          :: Q50
        real*8                                          :: SFDC
        real*8                                          :: RC
        real*8                                          :: VR
        real*8                                          :: HPC
        real*8                                          :: LFR
        real*8                                          :: RLD
        real*8                                          :: DLD
        real*8                                          :: SPDC
        real*8                                          :: BFI
        real*8                                          ::Qpeak10
        real*8                                          ::Qpeak50
        integer                                         :: SFDCcount
        integer                                         :: RCcount
        integer                                         :: VRcount
        integer                                         :: HPCcount
        integer                                         :: LFRcount
        real*8                                         :: WB
        real*8                                          :: maxBound
        real*8                                          :: minBound
        real*8                                          :: r2

       real*8,dimension(3000)                           :: Qtest
       real*8,dimension(3000)                           :: Qsimtest
       real*8,dimension(3000)                           :: Ptest
       real*8,dimension(2)                              :: tmp
       real*8,dimension(2)                           :: dummy
     integer                                            :: imax


!initialize


temp_date=dates_data(ic_start)
temp_date2=dates_data(ic_end)

read( temp_date(1:4) ,'(I4)') startyear
read(temp_date2(1:4) ,'(I4)') endyear 


nWindow = ceiling(real(endyear-startyear,8)/real(window))


uParam=14
uQ=uParam+1+nWindow
uObj=uQ+1+nWindow
uSig=uObj+1+nWindow
!uEa=uSig+10+nWindow
!uEi=uEa+1+nWindow
!uSu=uEi+1+nWindow
!uPeff=uSnow+1+nWindow

length=size(prec_data( iw_start:ic_end ))
allocate( q( length) )
allocate( count_feasibles( nWindow  ) )
allocate(prec_in (length+(iw_end-iw_start) ))
allocate(etp_in (length+(iw_end-iw_start) ))
allocate(temp_in (length+(iw_end-iw_start) ))

count_feasibles = 0
fIt = 0.1
paramset=par_ini

call init_random_seed()



window = window*365
print *, 'Windows:', nWindow


minBound = 0.7_8
maxBound = 1.3_8

SFDCcount =0
HPCcount=0
VRcount=0
RCcount=0
LFRcount=0

!start monte carlo runs

write(formHeader, *) '(2X a10, 2X a10, 2X a10, 2X a10, 2X a10, 2X a10)'
write(formHeaderSig, *) '(2X a10, 2X a10, 2X a10, 2X a10, 2X a10, 2X a10, 2X a10, 2X a10, 2X a10)'
write(formDataObj, *) '(2X f10.3, 2X f10.3, 2X f10.3, 2X f10.3, 2X f10.3, 2X f10.3)'
write(formDataSig, *) '(2X f10.5, 2X f10.5, 2X f10.5, 2X f10.5, 2X f10.5, 2X f10.5, 2X f10.5,2X f10.5, 2X f10.5)'


do iWindow =1, nWindow
           write(window_string,'(I0)') iWindow
           open(uParam+iWindow, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Parameters_")) // trim(adjustl(window_string))&
                  // trim(adjustl(".txt")), status='unknown', action='write')
           open(uQ+iWindow, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Q_"))// trim(adjustl(window_string))&
                  // trim(adjustl(".txt")), status='unknown', action='write')
           open(uObj+iWindow, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Objectives_")) // trim(adjustl(window_string))  &
             // trim(adjustl(".txt")),  status='unknown', action='write')
           open(uSig+iWindow, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Signatures_")) // trim(adjustl(window_string))  &
             // trim(adjustl(".txt")),  status='unknown', action='write')
!           open(uEa+iWindow, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Ea_"))// trim(adjustl(window_string))&
!             // trim(adjustl(".txt")),  status='unknown', action='write')
!           open(uEi+iWindow, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Ei_"))// trim(adjustl(window_string))&
!             // trim(adjustl(".txt")),  status='unknown', action='write')
!           open(uSu+iWindow, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Su_"))// trim(adjustl(window_string))&
!             // trim(adjustl(".txt")),  status='unknown', action='write')

 !           open(uPeff+iWindow, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Peff_"))// trim(adjustl(window_string))&
 !            // trim(adjustl(".txt")),  status='unknown', action='write')



           write(uParam+iwindow,'(14(1X a8))')  "Meltfactor", "Tthresh", "Imax", "Sumax", "beta", &
                              "Kf", "Ks", "LP", "D", "Pmax", "alpha","a", "b", "sumax_min"


           write(uObj+iWindow,formHeader) 'NSE','LNSE','VE', 'KGE', 'LKGE', 'WB'
           write(uSig+iWindow,formHeaderSig) 'RC','VarRatio','SFDC','HPC','Q90Q50','SPDC', 'BFI','RLD','DLD'

end do




 do nn=1, Iterations

	call random_number(r)

        forall(ii=1:size(paramset,1),optim(ii) .eqv. .TRUE.) paramset(ii)=r(ii)*(par_max(ii)-par_min(ii))+par_min(ii) 

        !if(optim(7) .eqv. .TRUE.) then
	!   call random_number(r2)
        !   paramset(7)=r2*(par_max(7)-paramset(6))+paramset(6)
        !end if

        !run the whole model
        !add artificial warmup period
        if(iw_start .eq. ic_start) then
           prec_in = (/ prec_data( iw_start:iw_end ), prec_data( iw_start:ic_end ) /)
           etp_in  = (/ etp_data( iw_start:iw_end ) , etp_data( iw_start:ic_end ) /)
           temp_in = (/ temp_data( iw_start:iw_end ), temp_data( iw_start:ic_end ) /)
           ic_end  =ic_end+iw_end 


        else
           prec_in =  prec_data( iw_start:ic_end )
           etp_in  =  etp_data( iw_start:ic_end ) 
           temp_in =  temp_data( iw_start:ic_end ) 
        end if



!	call model(paramset,output_ini(8:11,iw_end), prec_data( iw_start:ic_end ) &
!	,temp_data(iw_start:ic_end),etp_data(iw_start:ic_end), dem, cellsize, output)

        call model(paramset,incon, prec_in, temp_in, etp_in, dem, cellsize, output)


	q=output(1,:)
        extra=0
	!determine objectives for each window period
        iWindow = 1
        do jj= ic_start-1, ic_end-window, window

          start_ind = jj+1+extra
          end_ind   = jj+window+extra


          !check for uneven years
          temp_date=dates_data(start_ind)
          temp_date2=dates_data(end_ind+1)
          if(temp_date(6:10) .ne. temp_date2(6:10)) then
             end_ind=end_ind+1
             extra=extra+1
          end if

          if( (end_ind) .gt. ic_end) then 
            end_ind = ic_end
          end if

	   call NashEfficiency(Qobs_data(start_ind:end_ind),q(start_ind:end_ind), NSE)
	   call LogNashEfficiency(Qobs_data(start_ind:end_ind),q(start_ind:end_ind), LNSE)
	   call VolumeError(Qobs_data(start_ind:end_ind),q(start_ind:end_ind), VE)
           call KlingGupta(Qobs_data(start_ind:end_ind),q(start_ind:end_ind), KGE)
           call logKlingGupta(Qobs_data(start_ind:end_ind),q(start_ind:end_ind), LKGE)

       
           call WaterBalance(q(start_ind:end_ind),output(5,start_ind:end_ind),output(6,start_ind:end_ind), &
                             prec_in(start_ind:end_ind), WB)
	  ! call FDCNashEfficiency(Qobs_data(start_ind:end_ind),q(start_ind:end_ind), FDC_NSE)

	   call FlowDurationCurves(q(start_ind:end_ind), FDC, Q5, Q50, Q95, SFDC, LFR)
           call RunoffCoeff(q(start_ind:end_ind), prec_in(start_ind:end_ind), RC)
           call VarRatio(q(start_ind:end_ind), prec_in(start_ind:end_ind), VR)
           call HighPulsCount(q(start_ind:end_ind), HPC)
           call BaseFlowIndex(q(start_ind:end_ind), BFI)
           call PeakDistribution(q(start_ind:end_ind), SPDC, Qpeak10, Qpeak50)
           call Limb_densities  (q(start_ind:end_ind), RLD, DLD)

	   !keep parameter set when all signatures are between specified bounds

           !******** test mode *********

!          open(unit=30, file='/host/Documents/Sumax/Testdata_bristol/Qobs.txt'  , action='read')
!          open(unit=31, file='/host/Documents/Sumax/Testdata_bristol/Prec.txt'  , action='read')
!          open(unit=32, file='/host/Documents/Sumax/Testdata_bristol/Qsim.txt'  , action='read')

!           do i=1, 3000
!              read(30, *) Qtest(i)
!              read(31, *) Ptest(i)
!              read(32, *) Qsimtest(i)
!           end do

!          close(30)
!          close(31)
!          close(32)


!	   call NashEfficiency(Qtest,Qsimtest, NSE)
!	   call LogNashEfficiency(Qtest,Qsimtest, LNSE)
!	   call VolumeError(Qtest,Qsimtest, VE)
!	   call FDCNashEfficiency(Qtest,Qsimtest, FDC_NSE)
!          call KlingGupta(Qtest,Qsimtest, KGE)
!           call logKlingGupta(Qtest,Qsimtest, LKGE)

!	   call FlowDurationCurves(Qtest, FDC, Q5, Q50, Q95, SFDC, LFR)
!           call RunoffCoeff(Qtest, Ptest, RC)
!           call VarRatio(Qtest, Ptest, VR)
!           call HighPulsCount(Qtest, HPC)
!           call BaseFlowIndex(Qtest, BFI)
!           call PeakDistribution(Qtest, SPDC, Qpeak10, Qpeak50)
!           call Limb_densities  (Qtest, RLD, DLD)


           !***** end test **********



           write(uParam+iwindow,'(14(1X f8.2))')  paramset         
           write(uQ+iwindow,*)  output(1,start_ind:end_ind)
           write(uObj+iwindow,formDataObj) (/NSE, LNSE, VE, KGE, LKGE, WB/)
!           write(uEa+iWindow,*) output(5,start_ind:end_ind)
           write(uSig+iwindow,formDataSig) (/RC, VR, SFDC, HPC,LFR, SPDC, BFI, RLD, DLD/)
!           write(uEi+iWindow,*) output(6,start_ind:end_ind)
!           write(uSu+iWindow,*) output(9,start_ind:end_ind)
 !          write(uPeff+iWindow,*) output(4,start_ind:end_ind)
 
        iWindow = iWindow+1 



        end do

        if(nn .eq. int(fIt*Iterations)   ) then
        print *, "Iterations finished: ", nn
        !print *, "Feasible solutions : ", count_feasibles 
        print *, ""
        fIt = fIt + 0.1
      end if

 end do

do iWindow =1, nWindow
           close(uParam+iWindow)
           close(uQ+iWindow)
           close(uObj+iWindow)
           close(uSig+iWindow)
!           close(uEa+iWindow)
!           close(uEi+iWindow)
!           close(uSu+iWindow)
   !        close(uPeff+iWindow)
end do
 


!print *, "Feasible parameter sets:"
!print *, count_feasibles
!print *, ''

!print *, 'Overview Non-Feasibles'
!print *, '      RC    ','   VR       ','        SFDC   ','      HPC  ','     LFR  '
!print *, RCcount , VRcount, SFDCcount, HPCcount, LFRcount





print *, '----------------------------'
print *, ''
print *, "Optimization finished!"

end subroutine



end module mo_moving_window
