MODULE mo_eval_signatures

 
  IMPLICIT NONE

  PUBLIC :: eval_signatures	 ! signature evaluation subroutine

  PRIVATE



CONTAINS

 


  SUBROUTINE eval_signatures(Qm, Qobs, Prec, dates, write_output, EC)

	USE mo_signatures
	USE mo_readdata
	USE mo_objectives


	IMPLICIT NONE

	
  	real*8,  dimension(:),intent(in)                        :: Qm     	        ! simulated river runoff 
	real*8,  dimension(:),intent(in)                        :: Qobs         	! observed  river flow	
	real*8,  dimension(:),intent(in)                        :: Prec         	! precipitation	
  	logical,intent(in)                                      :: write_output     	! simulated river runoff
        character*10, dimension(:), intent(in)                  :: dates
	real*8, dimension(1:28), intent(out)   	                :: EC	 		! Evaluation criteria,
 
	real*8,dimension(1:200)   	      :: AC		 	! Autocorrelation,
    	real*8		   	      :: AC_high 		! Autocorrelation high flows,
   	real*8		   	      :: AC_low 		! Autocorrelation low flows,
	real*8,dimension(1:200)   	      :: ACobs		 	! Autocorrelation,    
    	real*8		   	      :: ACobs_high 		! Autocorrelation high flows,
   	real*8		   	      :: ACobs_low 		! Autocorrelation low flows,
	real*8 			      :: BFI   		        ! baseflow index  
	real*8 			      :: BFIobs   		! baseflow index observed  
	real*8, dimension(:,:), allocatable :: d_Qmod			! simulated river flow
    	real*8		   	      :: DLD	 		! Declining limb density,
    	real*8   	      		      :: DLDobs	 		! Declining limb density,

	real*8		   	      :: RLD	 		! Rising limb density,
	real*8		   	      :: RLDobs	 		! Rising limb density,
	character(256)                        :: fName,formDataF,formDataNS,formHeader
	character(256)                        :: formDataFDC,formHeaderFDC
	character(256)                        :: formDataAC,formHeaderAC
	real*8, allocatable, dimension(:,:) :: FDC			! Flow duration curve
	real*8, allocatable, dimension(:,:) :: FDC_high		! FDC high flows
	real*8, allocatable, dimension(:,:) :: FDC_low		! FDC low flows
	real*8, allocatable, dimension(:,:) :: FDCobs			! Flow duration curve
	real*8, allocatable, dimension(:,:) :: FDCobs_high		! FDC high flows
	real*8, allocatable, dimension(:,:) :: FDCobs_low		! FDC low flows     
	real*8			            :: HPC             	! high pulse count
	real*8			            :: HPCobs           ! high pulse count observed 
    	integer                           :: iday, iS, iE,jj	! Counters
	integer                           :: iTimer, ii, gg, tt	! Counters
	real*8, allocatable, dimension(:) :: months		! FDC low flows      
	integer                           :: nTimeSteps		! Number of timesteps
	real*8 			      :: peaks			! Peak distribution
	real*8 			      :: peaks_low 		! Peak dist low flows    
	real*8 			      :: peaks_high 		! Peak dist high flows
	real*8 			      :: peaksobs		! Peak distribution
	real*8 			      :: peaksobs_low 		! Peak dist low flows    
	real*8 			      :: peaksobs_high 		! Peak dist high flows        
	real*8 			      :: SFDC   		! slope flow duration curve  
	real*8 			      :: SFDC_obs   		! slope flow duration curve          
	real*8 			      :: LFR     		! lowflowratio        
	real*8 			      :: LFR_obs     		! lowflowratio        
	real*8			      :: Q_MA	         	! mean annual discharge,
	real*8			      :: Qobs_MA	        ! mean annual discharge,
	real*8			      :: Q5	        	! 5% occurence discharge,
	real*8			      :: Q50	        	! 50% occurence discharge,
	real*8			      :: Q95	        	! 95% occurence discharge,
	real*8			      :: Qobs5	        	! 5% occurence discharge,
	real*8			      :: Qobs50	        	! 50% occurence discharge,
	real*8			      :: Qobs95	        	! 95% occurence discharge,
	real*8			      :: Qlow5	        	! 5% occurence discharge,
	real*8			      :: Qlow50	        	! 50% occurence discharge,
	real*8			      :: Qlow95	        	! 95% occurence discharge,
	real*8			      :: Qlow_obs5        	! 5% occurence discharge,
	real*8			      :: Qlow_obs50        	! 50% occurence discharge,
	real*8			      :: Qlow_obs95        	! 95% occurence discharge,
	real*8			      :: Qhigh5	        	! 5% occurence discharge,
	real*8			      :: Qhigh50        	! 50% occurence discharge,
	real*8			      :: Qhigh95        	! 95% occurence discharge,
	real*8			      :: Qhigh_obs5        	! 5% occurence discharge,
	real*8			      :: Qhigh_obs50        	! 50% occurence discharge,
	real*8			      :: Qhigh_obs95        	! 95% occurence discharge,
	real*8			      :: Qpeak10        	! 50% occurence discharge,
	real*8			      :: Qpeak50        	! 95% occurence discharge,
	real*8			      :: Qpeak10_obs        	! 50% occurence discharge,
	real*8			      :: Qpeak50_obs        	! 95% occurence discharge,
	real*8			      :: Qlowpeak10        	! 50% occurence discharge,
	real*8			      :: Qlowpeak50        	! 95% occurence discharge,
	real*8			      :: Qlowpeak10_obs        	! 50% occurence discharge,
	real*8			      :: Qlowpeak50_obs        	! 95% occurence discharge,
	real*8			      :: Qhighpeak10        	! 50% occurence discharge,
	real*8			      :: Qhighpeak50        	! 95% occurence discharge,
	real*8			      :: Qhighpeak10_obs       	! 50% occurence discharge,
	real*8			      :: Qhighpeak50_obs       	! 95% occurence discharge,
	real*8			      :: VR             	! variance ratio
	real*8			      :: VRobs             	! variance ratio

        integer                       :: usignatures
        integer                       :: uFDC
        integer                       :: uAC
real :: start, finish

usignatures=15
uFDC=22
uAC=33 


!Create series of month-values for low-high flow calculations
call getmonths(dates, months)

! Calculate signatures for simulated flows

call Q_MeanAnnual                     (Qm, Q_MA, dates) 
call Autocorrelation                  (Qm, AC)
call Autocorrelation_low              (Qm, AC_low, months)
call Autocorrelation_high             (Qm, AC_high, months)
call Limb_densities                   (Qm, RLD, DLD)
call FlowDurationCurves               (Qm, FDC,       Q5,   Q50, Q95, SFDC, LFR)
call FlowDurationCurves_low           (Qm, FDC_low,   Qlow5, Qlow50, Qlow95, months)
call FlowDurationCurves_high          (Qm, FDC_high,  Qhigh5,Qhigh50, Qhigh95, months)
call PeakDistribution                 (Qm, peaks,      Qpeak10, Qpeak50)
call PeakDistribution_low             (Qm, peaks_low,  Qlowpeak10, Qlowpeak50, months)
call PeakDistribution_high            (Qm, peaks_high, Qhighpeak10, Qhighpeak50, months)
call VarRatio                         (Qm, Prec, VR)
call HighPulsCount                    (Qm, HPC)
call BaseFlowIndex                    (Qm, BFI)

! Calculate signatures for observed flows
call Q_MeanAnnual                    (Qobs, Qobs_MA, dates)
call Autocorrelation                 (Qobs, ACobs)
call Autocorrelation_low             (Qobs, ACobs_low, months)
call Autocorrelation_high            (Qobs, ACobs_high, months)
call Limb_densities                  (Qobs, RLDobs, DLDobs)
call FlowDurationCurves              (Qobs, FDCobs,      Qobs5,      Qobs50,       Qobs95, SFDC_obs, LFR_obs)
call FlowDurationCurves_low          (Qobs, FDCobs_low,  Qlow_obs5,  Qlow_obs50,   Qlow_obs95, months)
call FlowDurationCurves_high         (Qobs, FDCobs_high, Qhigh_obs5, Qhigh_obs50,  Qhigh_obs95, months)
call PeakDistribution                (Qobs, peaksobs,    Qpeak10_obs,Qpeak50_obs)
call PeakDistribution_low            (Qobs, peaksobs_low, Qlowpeak10_obs, Qlowpeak50_obs, months)
call PeakDistribution_high           (Qobs, peaksobs_high, Qhighpeak10_obs, Qhighpeak50_obs, months)
call VarRatio                         (Qobs, Prec, VRobs)
call HighPulsCount                   (Qobs, HPCobs)
call BaseFlowIndex                   (Qobs, BFIobs)

! Calculate evaluation criteria/efficiencies
call Rel_err(Qobs_MA, Q_MA, EC(1))
call Rel_err(ACobs(1), AC(1), EC(2))
call Rel_err(ACobs_low,AC_low, EC(3))
call Rel_err(ACobs_high, AC_high, EC(4))
call Rel_err(RLDobs, RLD, EC(5))
call Rel_err(DLDobs, DLD, EC(6))

call Rel_err(Qobs5, Q5, EC(7))
call Rel_err(Qobs50, Q50, EC(8))
call Rel_err(Qobs95, Q95, EC(9))

call Rel_err(Qlow_obs5, Qlow5, EC(10))
call Rel_err(Qlow_obs50, Qlow50, EC(11))
call Rel_err(Qlow_obs95, Qlow95, EC(12))

call Rel_err(Qhigh_obs5, Qhigh5, EC(13))
call Rel_err(Qhigh_obs50, Qhigh50, EC(14))
call Rel_err(Qhigh_obs95, Qhigh95,EC(15))


call Rel_err(peaksobs, peaks, EC(16))
call Rel_err(peaksobs_low, peaks_low,EC(17))
call Rel_err(peaksobs_high, peaks_high, EC(18))

call Rel_err(Qpeak10_obs, Qpeak10, EC(19))
call Rel_err(Qpeak50_obs, Qpeak50, EC(20))

call Rel_err(Qlowpeak10_obs, Qlowpeak10, EC(21))
call Rel_err(Qlowpeak50_obs, Qlowpeak50, EC(22))

call Rel_err(Qhighpeak10_obs, Qhighpeak10, EC(23))
call Rel_err(Qhighpeak50_obs, Qhighpeak50, EC(24))

call Rel_err(SFDC_obs, SFDC, EC(25))
call Rel_err(LFR_obs, LFR, EC(26))

call Rel_err(VRobs, VR, EC(27))
call Rel_err(HPCobs, HPC, EC(28))

call Rel_err(HPCobs, HPC, EC(28))

call NashEfficiency(FDCobs(:,2), FDC(:,2), EC(27))
call NashEfficiency(ACobs, AC, EC(28))

if( write_output) then

print *, "   Q_MA:     ", Q_MA 
print *, "   AC:       ", AC(1) 
print *, "   AClow:    ", AC_low
print *, "   AChigh:   ", AC_high
print *, "   RLD:      ", RLD
print *, "   DLD:      ", DLD
print *, "   Peaks:    ", peaks
print *, "   Peakslow: ", peaks_low
print *, "   Peakshigh:", peaks_high
print *, "   Q5:       ", Q5
print *, "   Q50:      ", Q50
print *, "   Q95:      ", Q95




!write signatures to file
fName=  trim(adjustl(Input_dir))//trim("Signatures.txt")
 open(usignatures, file=fName, status='unknown', action='write')
 write(formHeader, *) '(2X,   a10 , 2X,  a10, 2X,  a10, 2X,  a10, 2X, a10    )'
 write(formDataF, *) '(2X,   a10 , 2X,  f10.2, 2X,  f10.2, 2X, f10.2,  a10    )'
 write(formDataNS, *) '(2X,   a10 , 2X,  a10, 2X,  a10, 2X, a10, f10.2    )'

 write(usignatures,formHeader) 'Signature ', 'Obs ', 'Mod ', 'Rel.Err(-)','NSE(-)'

 write(usignatures,formDataF) 'Q_MA', Qobs_MA, Q_MA, EC(1),'-'
 write(usignatures,formDataF) 'AC', ACobs(1), AC(1), EC(2),'-'
 write(usignatures,formDataF) 'AC_low', ACobs_low, AC_low, EC(3),'-'
 write(usignatures,formDataF) 'AC_high', ACobs_high, AC_high, EC(4),'-'
 write(usignatures,formDataF) 'RLD', RLDobs, RLD, EC(5),'-'
 write(usignatures,formDataF) 'DLD', DLDobs, DLD, EC(6),'-'

 write(usignatures,formDataF) 'Q5', Qobs5, Q5, EC(7),'-'
 write(usignatures,formDataF) 'Q50', Qobs50, Q50, EC(8),'-'
 write(usignatures,formDataF) 'Q95', Qobs95, Q95, EC(9),'-'

 write(usignatures,formDataF) 'Q5_low', Qlow_obs5, Qlow5, EC(10),'-'
 write(usignatures,formDataF) 'Q50_low', Qlow_obs50, Qlow50, EC(11),'-'
 write(usignatures,formDataF) 'Q95_low', Qlow_obs95, Qlow95, EC(12),'-'

 write(usignatures,formDataF) 'Q5_high', Qhigh_obs5, Qhigh5, EC(13),'-'
 write(usignatures,formDataF) 'Q50_high', Qhigh_obs50, Qhigh50, EC(14),'-'
 write(usignatures,formDataF) 'Q95_high', Qhigh_obs95, Qhigh95, EC(15),'-'

 write(usignatures,formDataF) 'Peaks', peaksobs, peaks, EC(16),'-'
 write(usignatures,formDataF) 'Peaks_low', peaksobs_low, peaks_low, EC(17),'-'
 write(usignatures,formDataF) 'Peaks_high', peaksobs_high, peaks_high, EC(18),'-'

 write(usignatures,formDataF) 'Qpeak10', Qpeak10_obs, Qpeak10, EC(19),'-'
 write(usignatures,formDataF) 'Qpeak50', Qpeak50_obs, Qpeak50, EC(20),'-'

 write(usignatures,formDataF) 'Qlow_peak10', Qlowpeak10_obs, Qlowpeak10, EC(21),'-'
 write(usignatures,formDataF) 'Qlow_peak50', Qlowpeak50_obs, Qlowpeak50, EC(22),'-'

 write(usignatures,formDataF) 'Qhigh_peak10', Qhighpeak10_obs, Qhighpeak10, EC(23),'-'
 write(usignatures,formDataF) 'Qhigh_peak50', Qhighpeak50_obs, Qhighpeak50, EC(24),'-'

 write(usignatures,formDataNS) 'FDC_serie', '-', '-', '-',EC(25)
 write(usignatures,formDataNS) 'AC_serie', '-', '-', '-',EC(26)
 close(usignatures)

fName=  trim(adjustl(input_dir))//trim("FDC.txt")
 open(uFDC, file=fName, status='unknown', action='write')
 write(formHeaderFDC, *) '(2X,   a10 , 2X, a10, 2X, a2, 2X,  a10, 2X,  a10    )'
 write(formDataFDC, *) '(2X,   f10.2 , 2X,   f10.2, 2X, a2, 2X, f10.2, 2X,  f10.2)'


 write(uFDC,formHeaderFDC) 'Observed','','|','', 'Modelled'
 write(uFDC,formHeaderFDC) 'Percentage(%) ', 'Q(m3/s)','|', 'Percentage(%) ', 'Q(m3/s)'

do jj=1, size(FDCobs(:,1)) 
write(uFDC,formDataFDC) FDCobs(jj,1), FDCobs(jj,2), '|', FDC(jj,1), FDC(jj,2)
end do 
 close(uFDC)


fName=  trim(adjustl(input_dir))//trim("AC.txt")
 open(uac, file=fName, status='unknown', action='write')
 write(formHeaderAC, *) '(2X,   a10 , 2X, a10, 2X, a2, 2X,  a10, 2X,  a10    )'
 write(formDataAC, *) '(2X,   i10 , 2X,  f10.2, 2X, a2, 2X, i10, 2X,  f10.2)'


 write(uAC,formHeaderAC) 'Observed','','|','', 'Modelled'
 write(uAC,formHeaderAC) 'Tlag(days) ', 'AC(-)','|', 'Tlag(days) ', 'AC(-)'

do jj=1,200 
write(uAC,formDataAC) jj, ACobs(jj), '|', jj, AC(jj)
end do
 close(uac)

end if






 
  END SUBROUTINE eval_signatures



END MODULE mo_eval_signatures

