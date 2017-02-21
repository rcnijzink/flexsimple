MODULE mo_signatures


  IMPLICIT NONE


 PUBLIC :: Q_MeanAnnual			! Autocorrelation function
 PUBLIC :: Autocorrelation		! Autocorrelation function
 PUBLIC :: Autocorrelation_low		! Autocorrelation for low flows
 PUBLIC :: Autocorrelation_high		! Autocorrelation for high flows
 PUBLIC :: Limb_densities		! Rising and declining lomb densities
 PUBLIC :: FlowDurationCurves		! Flow duration curves
 PUBLIC :: FlowDurationCurves_low	! Flow duration curves low flows
 PUBLIC :: FlowDurationCurves_high	! Flow duration curves high flows
 PUBLIC :: PeakDistribution		! Peak distribution parameter
 PUBLIC :: PeakDistribution_low		! Peaks for low flows
 PUBLIC :: PeakDistribution_high	! Peaks for high flows
 PUBLIC :: RunoffCoeff          	! Runoff coefficient
 PUBLIC :: VarRatio                  	! Variance Ratio
 PUBLIC :: HighPulsCount                ! High Puls Count
 PUBLIC :: BaseFlowIndex                ! High Puls Count


  ! ------------------------------------------------------------------

CONTAINS


 
SUBROUTINE Q_MeanAnnual(Q_serie, Q_MA, dates_data)

     IMPLICIT NONE

 	real*8, dimension(:),intent(in)     :: Q_serie       ! simulated  daily river flow,
 	character*10, dimension(:),intent(in)     :: dates_data       ! simulated  daily river flow,    	    	

	integer                           :: nTimeSteps    ! number of timesteps
	real*8			      :: Qsum          ! sum of discharge 
	real*8, allocatable, dimension(:)   :: Q_Ann         ! Annual discharge
	real*8, intent(out)		      :: Q_MA          ! Mean annual discharge
	integer                           :: ii,jj         ! Counters
	character*10			                :: temp_date
	character*10			                :: temp2_date
	character*10			                :: temp_month
        real*8                                          :: month


allocate(Q_Ann(size(dates_data)/365+2))


Q_Ann = 0_8
Qsum = 0_8
jj=0
temp2_date = dates_data(1)
do ii=1, size(Q_serie(:))
	temp_date = dates_data (ii)
		if( temp_date(1:4) .eq. temp2_date(1:4) ) then
		Qsum=Qsum+Q_serie(ii)
			if(ii .eq. size(Q_serie(:))) then
			jj=jj+1
			Q_Ann(jj)=Qsum
			end if
		else
		jj=jj+1
		Q_Ann(jj)=Qsum
		Qsum=Q_serie(ii)
		end if
		
 temp2_date = temp_date

end do

	Q_MA=sum(Q_Ann(1:jj))/jj

deallocate(Q_Ann)
END SUBROUTINE Q_MeanAnnual

!-------------------------------------------------------------------------------
!         Autocorrelation

 
SUBROUTINE Autocorrelation(Q_serie, AC)

    IMPLICIT NONE

 	real*8, dimension(:),intent(in)           :: Q_serie       	! Simulated  daily river flow,
	real*8, dimension(1:200),intent(out)      :: AC		! Autocorrelation serie
	real*8			      :: Qmean		! Mean discharge
	real*8			      :: AC_up, AC_down	! Nominator and denominator
	integer                           :: ii,tlag	! Counters
	real*8, dimension(:), allocatable    :: tmp, tmp2	! Nominator and denominator


Qmean=sum(Q_serie)/size(Q_serie,1)
AC_up=0_8
AC_down=sum( (Q_serie - Qmean)**2)
allocate( tmp( size(Q_serie,1)))
allocate( tmp2( size(Q_serie,1)))



do tlag=1, 200

   tmp(1:size(Q_serie,1)-tlag)  = Q_serie(1:size(Q_serie,1)-tlag) - Qmean
   tmp2(1:size(Q_serie,1)-tlag) = Q_serie(1+tlag: size(Q_serie,1) )  - Qmean
  
   AC_up = sum(tmp(1:size(Q_serie,1)-tlag) * tmp2(1:size(Q_serie,1)-tlag))

   AC(tlag)=AC_up/AC_down
   AC_up = 0_8

   tmp   = 0_8
   tmp2  = 0_8
end do


!do tlag=1, 200

!   do ii=1, (size(Q_serie,1)-tlag)
!      AC_up=AC_up + (Q_serie(ii)-Qmean)*(Q_serie(ii+tlag)-Qmean)
!   end do

!AC(tlag)=AC_up/AC_down
!AC_up=0_8
!end do



END SUBROUTINE Autocorrelation

!-------------------------------------------------------------------------------
!         Autocorrelation_low


SUBROUTINE Autocorrelation_low(Q_serie, AC_low, months_serie)

USE mo_readdata

    IMPLICIT NONE

 	real*8, dimension(:),intent(in)           :: Q_serie       	! River flow,
 	real*8, dimension(:),intent(in)           :: months_serie         ! Dates,    	            
	real*8, intent(out)		          :: AC_low		! Autocorrelation 1 day low flows
	real*8			                  :: AC_down	        ! Denominator
	real*8			                  :: AC_up	        ! Nominator 
        integer, dimension(size(Q_serie,1) )      :: ind_low
	integer                                   :: jj                 ! Counters
        real*8                                    :: month
        real*8                                    :: month1
	integer                                   :: n_low		! Number of low flow days
	integer                                   :: n_ind		! Indeces of low flow days
	integer                                   :: tlag	        ! Lag time
	character*10			          :: temp_date
	character*10			          :: temp2_date
	character*10			          :: temp_month
	real*8			                  :: Qmean 	        ! Mean discharge
	real*8			                  :: Qsum 	        ! Discharge sum

AC_down = 0_8
AC_up   = 0_8
Qsum    = 0_8
n_low   = 0
ind_low = 0

!calculate the total mean of the low flows

if( mm_low .lt. mm_high) then

   do jj=1, size(Q_serie)
      !temp_date=dates_data(jj)
      !temp_month = temp_date(6:7)
      !read(temp_month,*) month

      if( (months_serie(jj) .ge. mm_low) .and. (months_serie(jj) .lt. mm_high) ) then

        n_low          = n_low+ 1
 	ind_low(n_low) = jj 
        Qsum           = Qsum + Q_serie(jj)

      end if

   end do

   Qmean=Qsum/real(n_low,8)

   !calculate AClow
   AC_down= sum( (Q_serie(ind_low(1:n_low) )-Qmean)**2)
   AC_up  = sum( (Q_serie(ind_low(1:n_low) )-Qmean)*(Q_serie(ind_low(1:n_low)+1)-Qmean))

   AC_low=AC_up/AC_down

else

   do jj=1, size(Q_serie)
      !temp_date=dates_data(jj)
      !temp_month = temp_date(6:7)
      !read(temp_month,*) month

      if( (months_serie(jj) .ge. mm_low) .or. (months_serie(jj) .lt. mm_high) ) then

        n_low          = n_low+ 1
 	ind_low(n_low) = jj 
        Qsum           = Qsum + Q_serie(jj)

       end if

   end do

   Qmean=Qsum/real(n_low,8)

   !calculate AClow
   AC_down= sum( (Q_serie(ind_low(1:n_low) )-Qmean)**2)
   AC_up  = sum( (Q_serie(ind_low(1:n_low) )-Qmean)*(Q_serie(ind_low(1:n_low)+1)-Qmean))

   AC_low=AC_up/AC_down


end if





END SUBROUTINE Autocorrelation_low

!-------------------------------------------------------------------------------
!     NAME
!         Autocorrelation_high


SUBROUTINE Autocorrelation_high(Q_serie, AC_high, months_serie)

USE mo_readdata

    IMPLICIT NONE

 	real*8, dimension(:),intent(in)           :: Q_serie    ! River flow,  
 	real*8, dimension(:),intent(in)           :: months_serie ! Dates,        	          
	real*8, intent(out)		          :: AC_high	! Autocorrelation 1 day low flows

	real*8			                  :: AC_down	! Denominator
	real*8			                  :: AC_up   	! Nominator 
        integer, dimension(size(Q_serie,1) )      :: ind_high
	integer                                   :: jj	        ! Counters
        real*8                                    :: month
        real*8                                    :: month1
	integer                                   :: n_high	! Number of high flow days
	integer                                   :: tlag	! Counters
	character*10			          :: temp_date
	character*10			          :: temp2_date
	character*10			          :: temp_month
	real*8			                  :: Qmean  	! Mean discharge	        
        real*8			                  :: Qsum 	! Discharge sum

AC_down  = 0_8
AC_up    = 0_8
Qsum     = 0_8
n_high   = 0
ind_high = 0

!calculate the total mean of the high flows

if(mm_low .lt. mm_high) then

   !calculate Qmean and indeces high flow
   do jj=1, size(Q_serie)
      !temp_date=dates_data(jj)
      !temp_month = temp_date(6:7)
      !read(temp_month,*) month

      if( (months_serie(jj) .lt. mm_low) .or. (months_serie(jj) .ge. mm_high) ) then

        n_high           = n_high + 1
 	ind_high(n_high) = jj 
        Qsum             = Qsum + Q_serie(jj)

      end if

   end do

   Qmean=Qsum/real(n_high,8)

   !calculate AChigh
   AC_down= sum( (Q_serie(ind_high(1:n_high) )-Qmean)**2)
   AC_up  = sum( (Q_serie(ind_high(1:n_high) )-Qmean)*(Q_serie(ind_high(1:n_high)+1)-Qmean))

   AC_high=AC_up/AC_down

else

   !calculate Qmean and indeces high flow
   do jj=1, size(Q_serie)
     ! temp_date=dates_data(jj)
     ! temp_month = temp_date(6:7)
     ! read(temp_month,*) month

      if( (months_serie(jj) .lt. mm_low) .and. (months_serie(jj) .ge. mm_high) ) then

        n_high           = n_high + 1
 	ind_high(n_high) = jj 
        Qsum             = Qsum + Q_serie(jj)

       end if
 
   end do

   Qmean=Qsum/real(n_high,8)

   !calculate AChigh
   AC_down= sum( (Q_serie(ind_high(1:n_high) )-Qmean)**2)
   AC_up  = sum( (Q_serie(ind_high(1:n_high) )-Qmean)*(Q_serie(ind_high(1:n_high)+1)-Qmean))

   AC_high=AC_up/AC_down

end if


END SUBROUTINE Autocorrelation_high

!-------------------------------------------------------------------------------
!         Limb_densities


SUBROUTINE Limb_densities(Q_serie, RLD, DLD)

    IMPLICIT NONE

 	real*8, dimension(:),intent(in)   :: Q_serie       	! River flow,    
	real*8, intent(out)		  :: RLD, DLD	! Rising and declining limb density
	real*8                            :: t_rise,t_fall	! Counters
	integer                           :: n_peak, jj     ! Number of peaks
	integer                           :: nn           ! Number of peaks
	integer                           :: nn_max           ! Number of peaks
	integer                           :: flat           ! Number of peaks
	integer,dimension(size(Q_serie))  :: mask_peak           ! Number of peaks
	logical                           :: foundpeak      ! Number of peaks
	logical                           :: found_rise      ! Number of peaks
	logical                           :: found_fall      ! Number of peaks

t_fall=0_8
t_rise=0_8

!calculate the total mean of the high flows

n_peak=0
foundpeak  = .FALSE.
found_rise = .FALSE.
found_fall = .FALSE.

flat = 1
!calculate the total mean of the high flows and count peaks
mask_peak = 0

!find maxima
do jj=2, size(Q_serie)-1
	if( (Q_serie(jj-1) .lt. Q_serie(jj)) .and. (Q_serie(jj+1) .lt. Q_serie(jj)) ) then
        n_peak=n_peak+1
        foundpeak = .TRUE.
        mask_peak(jj) = 1
        end if

       if(jj+20 .le. size(Q_serie)) then
          nn_max = 20
       else
          nn_max = size(Q_serie) - jj
       end if

              do nn=2, nn_max
                 if((Q_serie(jj-1) .lt. Q_serie(jj)) .and. (Q_serie(jj+nn-1) .eq. Q_serie(jj)) .and. (foundpeak .eqv. .FALSE.) )then

                     flat = flat + 1
                     if( (Q_serie(jj+nn) .lt. Q_serie(jj)) .and. (flat .eq. nn ) ) then
                               n_peak=n_peak+1
                               foundpeak = .TRUE.
                               mask_peak(jj:(jj+nn-1)) = 1
                     end if
                  end if
              end do

        flat = 1
        foundpeak = .FALSE.

end do

!find trise and tfall

do jj=1, size(Q_serie)-1
!do jj=1, size(Q_serie)-1
	if( (Q_serie(jj+1) .gt. Q_serie(jj)) ) then
        t_rise = t_rise+1
        found_rise = .TRUE.
        end if

        if( (Q_serie(jj+1) .lt. Q_serie(jj)) ) then
        t_fall = t_fall+1
        found_fall = .TRUE.
        end if

       if(jj+20 .le. size(Q_serie)) then
          nn_max = 20
       else
          nn_max = size(Q_serie) - jj
       end if

              do nn=2, nn_max
                 if((Q_serie(jj-1) .lt. Q_serie(jj)) .and. (Q_serie(jj+nn-1) .eq. Q_serie(jj)) &
                     .and.(found_rise .eqv. .FALSE.) .and. (flat .eq. (nn-1) )                  )then

                     flat = flat + 1
                     if( (Q_serie(jj+nn) .gt. Q_serie(jj))  ) then
                               t_rise = t_rise+(nn-1)
                               found_rise = .TRUE.
                     end if
                  end if
                  if((Q_serie(jj-1) .gt. Q_serie(jj)) .and. (Q_serie(jj+nn-1) .eq. Q_serie(jj)) &
                     .and.(found_fall .eqv. .FALSE.) .and. (flat .eq. (nn-1) )                  )then

                     flat = flat + 1
                     if( (Q_serie(jj+nn) .lt. Q_serie(jj))  ) then
                               t_fall = t_fall+(nn-1)
                               found_fall = .TRUE.
                     end if
                  end if


              end do

        flat = 1
        found_rise = .FALSE.
        found_fall = .FALSE.
end do

!t_fall=0_8
!flat = 1
!found_fall = .FALSE.

!find tfall
!do jj=1, size(Q_serie)-1
!	if( (Q_serie(jj+1) .lt. Q_serie(jj)) ) then
!        t_fall = t_fall+1
!        found_fall = .TRUE.
!        end if

!       if(jj+20 .le. size(Q_serie)) then
!          nn_max = 20
!       else
!          nn_max = size(Q_serie) - jj
!       end if

!              do nn=2, nn_max
!                 if((Q_serie(jj-1) .gt. Q_serie(jj)) .and. (Q_serie(jj+nn-1) .eq. Q_serie(jj)) &
!                     .and.(found_fall .eqv. .FALSE.) .and. (flat .eq. (nn-1) )                  )then

!                     flat = flat + 1
!                     if( (Q_serie(jj+nn) .lt. Q_serie(jj))  ) then
!                               t_fall = t_fall+(nn-1)
!                               found_fall = .TRUE.
!                     end if
!                  end if
!              end do

!        flat = 1
!        found_fall = .FALSE.

!end do


RLD=real(n_peak,8)/real(t_rise,8)
DLD=real(n_peak,8)/real(t_fall,8)




END SUBROUTINE Limb_densities
!-------------------------------------------------------------------------------
!         FlowDurationCurves


SUBROUTINE FlowDurationCurves(Q_serie, FDC_serie, Q5, Q50, Q95, slope, lowflow_ratio)

	
    use mo_sort
    IMPLICIT NONE

 	real*8, dimension(:),intent(in)   	:: Q_serie      ! River flow,
 	real*8,intent(out)   			:: Q5           ! 5% river flow, 
 	real*8,intent(out)   			:: Q50          ! 50% river flow, 
 	real*8,intent(out)   			:: Q95          ! 95% river flow,     
	real*8, allocatable, dimension(:,:),intent(out)	:: FDC_serie	! ordened Q
        real*8, optional,intent(out)            :: slope        !slope between Q66 and Q33
        real*8, optional,intent(out)            :: lowflow_ratio!ratio between Q90 and Q50

	real*8, allocatable,dimension(:)      :: per,Q_sort  	! percentages and sorted Q
	integer                               :: ii		! Counter
	integer                               :: n		! Size of array
	real*8, allocatable,dimension(:)      :: Q_quantiles	!Q with p=10% and 50%
	real*8,dimension(3)		      :: p		!Q with p=90%
	real*8,dimension(2)		      :: pslope	!Q with p=90%
	real*8, allocatable,dimension(:)      :: Qslope         !Q with p=33% and 66%
	real*8,dimension(2)		      :: plow	        !p=90% and 50%
	real*8, allocatable,dimension(:)      :: Qlow           !Q with p=50% and 90%

       

       	allocate( per(size(Q_serie)) ) 
       	allocate( Q_sort(size(Q_serie)) ) 
	allocate( FDC_serie(size(Q_serie),2)) 

lowflow_ratio = 0_8

Q_sort=Q_serie

call sort(Q_sort)

n=size(Q_sort)
Q_sort=Q_sort(n:1:-1)


do ii=1, size(Q_serie)
per(ii)=real(ii,8)/size(Q_serie)*100
end do

p=(/0.05_8,0.5_8,0.95_8/)
call Quantile(Q_sort,p,Q_quantiles)

Q5=Q_quantiles(1)
Q50=Q_quantiles(2)
Q95=Q_quantiles(3)

FDC_serie(:,1)=per
FDC_serie(:,2)=Q_sort

pslope =(/0.33_8, 0.66_8/)

call Quantile(Q_sort,pslope, Qslope)


slope = (log(Qslope(1)) - log(Qslope(2))) / ((pslope(2) - pslope(1))*100_8)


plow = (/0.5_8, 0.9_8/)
call Quantile(Q_sort,plow, Qlow)

lowflow_ratio = Qlow(2)/Qlow(1)



       	deallocate( per ) 
       	deallocate( Q_sort ) 
        

END SUBROUTINE FlowDurationCurves

!-------------------------------------------------------------------------------
!         FlowDurationCurves_low



SUBROUTINE FlowDurationCurves_low(Q_serie, FDClow_serie, Qlow5, Qlow50, Qlow95, months_serie)

        USE mo_readdata
    	USE mo_sort         		

    IMPLICIT NONE

 	real*8, dimension(:),intent(in)                 :: Q_serie     ! River flow,
 	real*8, dimension(:),intent(in)                 :: months_serie! Dates,        	          
 	real*8,intent(out)   			        :: Qlow5       ! River flow, 
 	real*8,intent(out)   			        :: Qlow50      ! River flow, 
 	real*8,intent(out)   			        :: Qlow95      ! River flow,    
	real*8, allocatable, dimension(:,:),intent(out)	:: FDClow_serie	! ordened Q 
	
	integer                                         :: jj	! Counters
	integer                                         :: kk	! Counters
        real*8                                          :: month
	integer                                         :: n! Array sizes
	integer                                         :: n_low! Array sizes
	real*8,dimension(3)		                :: p		!Q with p=90%
	real*8, allocatable,dimension(:)                :: per  	! percentages and sorted Q
	character*10			                :: temp_date
	character*10			                :: temp2_date
	character*10			                :: temp_month
        real*8,dimension(:),allocatable                 :: Qlow         ! discharge in low flow period
	real*8,dimension(size(Q_serie,1))               :: Qlow_tmp     ! discharge in low flow period
	real*8, allocatable,dimension(:)                :: Q_sort  	! sorted Q
	real*8, allocatable,dimension(:)                :: Q_quantiles	!Q with p=10% and 50%


n_low=0



if(mm_low .lt. mm_high) then
!determine number of low flow days

do jj=1, size(Q_serie)
!temp_date=dates_data(jj)
!temp_month = temp_date(6:7)
!   read(temp_month,*) month
	if( (months_serie(jj) .ge. mm_low) .and. (months_serie(jj) .lt. mm_high) ) then
	n_low=n_low+1
	Qlow_tmp(n_low)=Q_serie(jj)
	end if

end do

!Create low flow array
allocate(Qlow(n_low))

Qlow = Qlow_tmp(1:n_low)


else

do jj=1, size(Q_serie)
!temp_date=dates_data(jj)
!temp_month = temp_date(6:7)
!   read(temp_month,*) month
	if( (months_serie(jj) .ge. mm_low) .or. (months_serie(jj) .lt. mm_high) ) then
	n_low=n_low+1
	Qlow_tmp(n_low)=Q_serie(jj)
	end if

end do



!Create low flow array
allocate(Qlow(n_low))
Qlow = Qlow_tmp(1:n_low)


end if


       	allocate( per(n_low) ) 
       	allocate( Q_sort(n_low) ) 
	allocate( FDClow_serie(n_low,2)) 


!Create low flow duration curve
Q_sort = Qlow
call sort(Q_sort)



do kk=1, size(Q_sort)
per(kk)=real(kk)/size(Q_sort)*100
end do

n=size(Q_sort)
Q_sort=Q_sort(n:1:-1)

FDClow_serie(:,1)=per
FDClow_serie(:,2)=Q_sort

p=(/0.05,0.5,0.95/)
call Quantile(Q_sort,p,Q_quantiles)

Qlow5=0
Qlow50=0
Qlow95=0


Qlow5=Q_quantiles(1)
Qlow50=Q_quantiles(2)
Qlow95=Q_quantiles(3)




       	deallocate( per ) 
       	deallocate( Q_sort ) 



END SUBROUTINE FlowDurationCurves_low

!-------------------------------------------------------------------------------
!         FlowDurationCurves_high


SUBROUTINE FlowDurationCurves_high(Q_serie, FDChigh_serie, Qhigh5, Qhigh50, Qhigh95, months_serie)
       
        USE mo_readdata
    	USE mo_sort         		

    IMPLICIT NONE

 	real*8, dimension(:),intent(in)                 :: Q_serie      ! River flow,
 	real*8, dimension(:),intent(in)                 :: months_serie   ! Dates, 	          
 	real*8,intent(out)   		                :: Qhigh5       ! 5% high river flow, 
 	real*8,intent(out)   		                :: Qhigh50      ! 50% high river flow, 
 	real*8,intent(out)   		                :: Qhigh95      ! 95% high river flow,    
	real*8, allocatable, dimension(:,:),intent(out)	:: FDChigh_serie !ordened Q 
	integer                                         :: jj    	! counters
 	integer                                         :: kk	        ! counters
        real*8                                          :: month
	integer                                         :: n             ! Array sizes
	integer                                         :: n_high        ! Array sizes
	real*8,dimension(3)		                :: p		! probablity of occurence
	real*8, allocatable,dimension(:)                :: per   	! percentages and sorted Q
	character*10			                :: temp_date
	character*10			                :: temp2_date
	character*10			                :: temp_month
	real*8,dimension(:),allocatable                 :: Qhigh
	real*8,dimension(size(Q_serie,1))               :: Qhigh_tmp
	real*8, allocatable,dimension(:)                :: Q_sort           ! sorted Q
	real*8, allocatable,dimension(:)                :: Q_quantiles	! Q with p=5%, p=10% and 50%


n_high=0


if(mm_low .lt. mm_high) then

!determine number of high flow days
do jj=1, size(Q_serie)
!temp_date=dates_data(jj)
!temp_month = temp_date(6:7)
!   read(temp_month,*) month
	if( (months_serie(jj) .lt. mm_low) .or. (months_serie(jj) .ge. mm_high) ) then
	n_high=n_high+1
	Qhigh_tmp(n_high)=Q_serie(jj)
	end if

end do


!Create high flow array
allocate(Qhigh(n_high))
Qhigh = Qhigh_tmp(1:n_high)



else

!determine number of high flow days
do jj=1, size(Q_serie)
!temp_date=dates_data(jj)
!temp_month = temp_date(6:7)
!   read(temp_month,*) month
	if( (months_serie(jj) .lt. mm_low) .and. (months_serie(jj) .ge. mm_high) ) then
	n_high=n_high+1
	Qhigh_tmp(n_high)=Q_serie(jj)
	end if

end do


!Create high flow array
allocate(Qhigh(n_high))
Qhigh = Qhigh_tmp(1:n_high)


end if

       	allocate( per(size(Qhigh)) ) 
       	allocate( Q_sort(size(Qhigh)) ) 
	allocate( FDChigh_serie(size(Qhigh),2)) 


!Create high flow duration curve


Q_sort=Qhigh

call sort(Q_sort)


do kk=1, size(Q_sort)
per(kk)=real(kk)/size(Q_sort)*100
end do

n=size(Q_sort)
Q_sort=Q_sort(n:1:-1)

p=(/0.05,0.5,0.95/)
call Quantile(Q_sort,p,Q_quantiles)

Qhigh5=Q_quantiles(1)
Qhigh50=Q_quantiles(2)
Qhigh95=Q_quantiles(3)


FDChigh_serie(:,1)=per
FDChigh_serie(:,2)=Q_sort



       	deallocate( per ) 
       	deallocate( Q_sort ) 

END SUBROUTINE FlowDurationCurves_high

!-------------------------------------------------------------------------------
!         PeakDistribution


SUBROUTINE PeakDistribution(Q_serie, peaks, Qpeak10, Qpeak50)

    	USE mo_sort

    IMPLICIT NONE

 	real*8, dimension(:),intent(in)     :: Q_serie       	! River flow,    
	real*8, intent(out)		      :: peaks   	! slope of peak flow duration curve
	real*8,dimension(:),allocatable     :: Qpeak       	! peak distribution discharge
	real*8, intent(out)	    	      :: Qpeak10       	! 10% peak distribution discharge
	real*8, intent(out)	    	      :: Qpeak50       	! 50% peak distribution discharge
	integer                           :: n_peak         ! Number of peaks
	integer                           :: nn_max         ! Number of peaks

	real*8, allocatable,dimension(:)    :: per     	! percentages
	real*8, allocatable,dimension(:)    :: Q_sort         ! sorted Q
	integer                           :: n	        ! Array sizes
        integer                           :: jj, ii,kk, nn  ! counters
	real*8, allocatable, dimension(:,:) :: FDCpeak	! all FDC
	real*8, allocatable,dimension(:)    :: Q_quantiles	! Q with p=10% and 50%
	real*8,dimension(2)		      :: p		! probility array
        logical                             :: foundpeak
        real*8                              :: flat

n_peak=0
foundpeak = .FALSE.
flat = 1
!calculate the total mean of the high flows and count peaks

do jj=2, size(Q_serie)-1
	if( (Q_serie(jj-1) .lt. Q_serie(jj)) .and. (Q_serie(jj+1) .lt. Q_serie(jj)) ) then
        n_peak=n_peak+1
        foundpeak = .TRUE.
        end if

       if(jj+20 .le. size(Q_serie)) then
          nn_max = 20
       else
          nn_max = size(Q_serie) - jj
       end if

              do nn=2, nn_max
                 if((Q_serie(jj-1) .lt. Q_serie(jj)) .and. (Q_serie(jj+nn-1) .eq. Q_serie(jj)) .and. (foundpeak .eqv. .FALSE.) )then
                     flat = flat + 1
                     if( (Q_serie(jj+nn) .lt. Q_serie(jj)) .and. (flat .eq. nn ) ) then
                               n_peak=n_peak+1
                               foundpeak = .TRUE.
                     end if
                  end if
              end do
        flat = 1
        foundpeak = .FALSE.

end do


allocate(Qpeak(n_peak))

! find peaks
kk=0
do ii=2, size(Q_serie)-1
	if( (Q_serie(ii-1) .lt. Q_serie(ii)) .and. (Q_serie(ii+1) .lt. Q_serie(ii)) ) then
	kk=kk+1	
	Qpeak(kk)=Q_serie(ii)
        foundpeak = .TRUE.
	end if

      if(ii+20 .le. size(Q_serie)) then
          nn_max = 20
       else
          nn_max = size(Q_serie) - ii
       end if

        do nn=2, nn_max
                 if((Q_serie(ii-1) .lt. Q_serie(ii)) .and.  (Q_serie(ii+nn-1)) .eq. Q_serie(ii) .and. (foundpeak .eqv. .FALSE.))then
                     flat = flat + 1
                     if( (Q_serie(ii+nn) .lt. Q_serie(ii) ).and. (flat .eq. nn )) then
                       kk=kk+1	
	               Qpeak(kk)=Q_serie(ii)
                       foundpeak = .TRUE.
                     end if
                  end if
         end do

        foundpeak = .FALSE.
        flat = 1



end do

       	allocate( per(size(Qpeak)) ) 
       	allocate( Q_sort(size(Qpeak)) ) 


!Create peak flow duration curve


Q_sort=Qpeak

call sort(Q_sort)

do nn=1, size(Q_sort)
per(nn)=real(nn)/size(Q_sort)*100
end do

n=size(Q_sort)
Q_sort=Q_sort(n:1:-1)

!calculate slope between 10 and 50 percent quantiles
p=(/0.1,0.5/)
call Quantile(Q_sort,p,Q_quantiles)


peaks=(log(Q_quantiles(1))-log(Q_quantiles(2)))/( (p(2)-p(1))*100  )


Qpeak10=Q_quantiles(1)
Qpeak50=Q_quantiles(2)


       	deallocate( per ) 
       	deallocate( Q_sort ) 
	deallocate( Qpeak)


END SUBROUTINE PeakDistribution

!-------------------------------------------------------------------------------
!         PeakDistribution_low


SUBROUTINE PeakDistribution_low(Q_serie, peaks_low, Qpeak10, Qpeak50, month_serie)

        USE mo_readdata
    	USE mo_sort

    IMPLICIT NONE

 	real*8, dimension(:),intent(in)     :: Q_serie       	! River flow,   
 	real*8, dimension(:),intent(in)     :: month_serie      ! simulated  daily river flow,    	    	           
	real*8, intent(out)		      :: peaks_low   	! slope of peak flow duration curve
	real*8,dimension(:),allocatable     :: Qpeak       	! Peak discharge low flows
	real*8, intent(out)	    	      :: Qpeak10       	! 10% low peak distribution discharge
	real*8, intent(out)	    	      :: Qpeak50       	! 50% low peak distribution discharge
	integer                           :: n_peak         ! Number of peaks
        integer                           :: jj, ii,kk, nn  ! counters

	real*8, allocatable,dimension(:)    :: per,Q_sort  	! percentages and sorted Q
	integer                           :: n	        ! Array sizes
	real*8, allocatable,dimension(:)    :: Q_quantiles	! Q with p=10% and 50%
	real*8,dimension(2)		      :: p		! array of probabilities
	character*10			                :: temp_date
	character*10			                :: temp2_date
	character*10			                :: temp_month
        real*8                                          :: month



n_peak=0
kk=0







if(mm_low .lt. mm_high) then

do jj=2, size(Q_serie)-1
!temp_date=dates_data(jj)
!temp_month = temp_date(6:7)
!   read(temp_month,*) month
	if( (Q_serie(jj-1) .le. Q_serie(jj)) .and. (Q_serie(jj+1) .le. Q_serie(jj)) &
! 	.and. (month .ge. mm_low) .and. (month .lt. mm_high) ) then
 	.and. (month_serie(jj) .ge. mm_low) .and. (month_serie(jj) .lt. mm_high) ) then
	n_peak=n_peak+1
	end if

end do

allocate(Qpeak(n_peak))


do ii=2, size(Q_serie)-1
!temp_date=dates_data(ii)
!temp_month = temp_date(6:7)
!read(temp_month,*) month
	if( (Q_serie(ii-1) .le. Q_serie(ii)) .and. (Q_serie(ii+1) .le. Q_serie(ii)) &
!	.and. (month .ge. mm_low) .and. (month .lt. mm_high) ) then
.and. (month_serie(ii) .ge. mm_low) .and. (month_serie(ii) .lt. mm_high) ) then
	kk=kk+1	
	Qpeak(kk)=Q_serie(ii)
	end if

end do

else

do jj=2, size(Q_serie)-1
!temp_date=dates_data(jj)
!temp_month = temp_date(6:7)
!   read(temp_month,*) month
	if( (Q_serie(jj-1) .le. Q_serie(jj)) .and. (Q_serie(jj+1) .le. Q_serie(jj)) &
! 	.and. ((month .ge. mm_low) .or. (month .lt. mm_high)) ) then
 	.and. ((month_serie(jj) .ge. mm_low) .or. (month_serie(jj) .lt. mm_high)) ) then
	n_peak=n_peak+1
	end if

end do

allocate(Qpeak(n_peak))




do ii=2, size(Q_serie)-1
!temp_date=dates_data(ii)
!temp_month = temp_date(6:7)
!read(temp_month,*) month
	if( (Q_serie(ii-1) .le. Q_serie(ii)) .and. (Q_serie(ii+1) .le. Q_serie(ii)) &
! 	.and. ((month .ge. mm_low) .or. (month .lt. mm_high)) ) then
 	.and. ((month_serie(jj) .ge. mm_low) .or. (month_serie(jj) .lt. mm_high)) ) then
	kk=kk+1	
	Qpeak(kk)=Q_serie(ii)
	end if

end do

end if






       	allocate( per(size(Qpeak)) ) 
       	allocate( Q_sort(size(Qpeak)) ) 



!Create peak flow duration curve

Q_sort=Qpeak

call sort(Q_sort)

do nn=1, size(Q_sort)
per(nn)=real(nn)/size(Q_sort)*100
end do

n=size(Q_sort)
Q_sort=Q_sort(n:1:-1)


!calculate slope between 10 and 50 percent quantiles
p=(/0.1,0.5/)
call Quantile(Q_sort,p,Q_quantiles)

peaks_low=(Q_quantiles(1)-Q_quantiles(2))/(0.9-0.5)

Qpeak10=Q_quantiles(1)
Qpeak50=Q_quantiles(2)

       	deallocate( per ) 
       	deallocate( Q_sort ) 
	deallocate( Qpeak )


END SUBROUTINE PeakDistribution_low

!-------------------------------------------------------------------------------
!         PeakDistribution_high

SUBROUTINE PeakDistribution_high(Q_serie, peaks_high, Qpeak10, Qpeak50, months_serie)

        USE mo_readdata
    	USE mo_sort

    IMPLICIT NONE   	

 	real*8, dimension(:),intent(in)     :: Q_serie       	! River flow   
 	real*8, dimension(:),intent(in)     :: months_serie       ! simulated  daily river flow,    	    	           
	real*8, intent(out)		      :: peaks_high   	! slope of peak flow duration curve
	real*8,dimension(:),allocatable     :: Qpeak       	! Peak discharge low flows
	real*8, intent(out)	    	      :: Qpeak10       	! 10% high peak distribution discharge
	real*8, intent(out)	    	      :: Qpeak50       	! 50% high peak distribution discharg
	integer                           :: n_peak         ! Number of peaks
        integer                           :: jj, ii,kk, nn  ! counters

	real*8, allocatable,dimension(:)    :: per,Q_sort  	! percentages and sorted Q
	integer                           :: n	        ! Array sizes
	real*8, allocatable,dimension(:)    :: Q_quantiles	! Q with p=10% and 50%
	real*8,dimension(2)		      :: p		! array of probabilities
	character*10			                :: temp_date
	character*10			                :: temp2_date
	character*10			                :: temp_month
        real*8                                          :: month

n_peak=0
kk=0


if( mm_low .lt. mm_high) then
! count peaks
do jj=2, size(Q_serie)-1
!temp_date=dates_data(jj)
!temp_month = temp_date(6:7)
!read(temp_month,*) month
	if( (Q_serie(jj-1) .le. Q_serie(jj)) .and. (Q_serie(jj+1) .le. Q_serie(jj)) &
 	.and. ((months_serie(jj) .lt. mm_low) .or. (months_serie(jj) .ge. mm_high)) ) then
	n_peak=n_peak+1
	end if

end do

allocate(Qpeak(n_peak))


do ii=2, size(Q_serie)-1
!temp_date=dates_data(ii)
!temp_month = temp_date(6:7)
!read(temp_month,*) month
	if( (Q_serie(ii-1) .le. Q_serie(ii)) .and. (Q_serie(ii+1) .le. Q_serie(ii)) &
	.and. ((months_serie(ii) .lt. mm_low) .or. (months_serie(ii) .ge. mm_high)) ) then
	kk=kk+1	
	Qpeak(kk)=Q_serie(ii)
	end if

end do

else

! count peaks
do jj=2, size(Q_serie)-1
!temp_date=dates_data(jj)
!temp_month = temp_date(6:7)
!read(temp_month,*) month
	if( (Q_serie(jj-1) .le. Q_serie(jj)) .and. (Q_serie(jj+1) .le. Q_serie(jj)) &
 	.and. ((months_serie(jj) .lt. mm_low) .and. (months_serie(jj) .ge. mm_high)) ) then
	n_peak=n_peak+1
	end if

end do

allocate(Qpeak(n_peak))


do ii=2, size(Q_serie)-1
!temp_date=dates_data(ii)
!temp_month = temp_date(6:7)
!read(temp_month,*) month
	if( (Q_serie(ii-1) .le. Q_serie(ii)) .and. (Q_serie(ii+1) .le. Q_serie(ii)) &
	.and. ((months_serie(ii) .lt. mm_low) .and. (months_serie(ii) .ge. mm_high)) ) then
	kk=kk+1	
	Qpeak(kk)=Q_serie(ii)
	end if

end do

end if



       	allocate( per(size(Qpeak)) ) 
       	allocate( Q_sort(size(Qpeak)) ) 



!Create peak flow duration curve


Q_sort=Qpeak

call sort(Q_sort)

do nn=1, size(Q_sort)
per(nn)=real(nn)/size(Q_sort)*100
end do

n=size(Q_sort)
Q_sort=Q_sort(n:1:-1)

!calculate slope between 10 and 50 percent quantiles
p=(/0.1,0.5/)
call Quantile(Q_sort,p,Q_quantiles)

peaks_high=(Q_quantiles(1)-Q_quantiles(2))/(0.9-0.5)

Qpeak10=Q_quantiles(1)
Qpeak50=Q_quantiles(2)


       	deallocate( per ) 
       	deallocate( Q_sort ) 
	deallocate( Qpeak )

END SUBROUTINE PeakDistribution_high

!-------------------------------------------------------------------------------
!         RunoffCoeff


SUBROUTINE RunoffCoeff(Q_serie, P_serie, RC)

	
    use mo_sort
    IMPLICIT NONE

 	real*8, dimension(:),intent(in)   	:: Q_serie      ! River flow,
 	real*8, dimension(:),intent(in)   	:: P_serie      ! River flow,
 	real*8,intent(out)   			:: RC           ! 5% river flow, 
 	

       
          RC    = sum(Q_serie) /  (sum(P_serie))


END SUBROUTINE RunoffCoeff

!-------------------------------------------------------------------------------
!         VarRatio


SUBROUTINE VarRatio(Q_serie, P_serie, VR)


    IMPLICIT NONE

 	real*8, dimension(:),intent(in)   	:: Q_serie      ! River flow,
 	real*8, dimension(:),intent(in)   	:: P_serie      ! River flow,
 	real*8,intent(out)   			:: VR           ! 5% river flow, 
        real*8                                  :: var_p        ! variance precipitation
        real*8                                  :: var_q        ! variance discharge
 	

          var_p     = sum( (P_serie  - (sum( P_serie )/real( size( P_serie,1) ,8) )  )**2)
          var_q     = sum( (Q_serie  - (sum( Q_serie )/real( size( Q_serie,1) ,8) )  )**2)
       
          VR    = var_q /  var_p


END SUBROUTINE VarRatio

!-------------------------------------------------------------------------------
!         HighPulsCount

SUBROUTINE HighPulsCount(Q_serie, HPC)

    use mo_sort

    IMPLICIT NONE

 	real*8, dimension(:),intent(in)   	:: Q_serie      ! River flow,
 	real*8,intent(out)   			:: HPC           ! 5% river flow, 
        real*8                                  :: eventCount   ! variance precipitation
        real*8, allocatable,dimension(:)        :: tmp
        real*8                                  :: thresh
        integer                                 :: i
        integer                                 :: n
        real*8, dimension(size(Q_serie,1))        :: Q_sort
 	
eventCount = 0
HPC = 0

Q_sort = Q_serie
call sort(Q_sort)

n=size(Q_sort)
Q_sort=Q_sort(n:1:-1)


call Quantile(Q_sort,(/0.5_8, 0.5_8/),tmp) 
thresh = tmp(1) *3_8

	  do i = 1, size(Q_serie,1)
		  if(Q_serie(i) .gt. thresh) then
			  if(Q_serie(i-1) .lt. thresh) then
				  eventCount = eventCount + 1_8
			  end if
		  end if
	  end do

         if(eventCount .gt. 0) HPC =real(eventCount,8)/ (real(size(Q_serie,1),8)/365.25_8)



END SUBROUTINE HighPulsCount
!-------------------------------------------------------------------------------
!         Baseflow Index

SUBROUTINE BaseFlowIndex(Q_serie, BFI)

    use mo_sort

    IMPLICIT NONE

 	real*8, dimension(:),intent(in)   	:: Q_serie      ! River flow,
 	real*8,intent(out)   			:: BFI          ! Base flow index
        
        integer                                 :: i
        integer                                 :: n
        integer                                 :: length
        real*8                                 :: a
        real*8                                 :: c1
        real*8                                 :: c2
        real*8, dimension(size(Q_serie,1))        :: Q_sort
        real*8, dimension(size(Q_serie,1))        :: baseflow
 	
a = 0.925

length = size(Q_serie)

  baseflow(1) = 0.5*Q_serie(1)
  

  
  
  do i=2, length
    c1 = a*baseflow(i-1)
    c2 = (1-a)*0.5*( Q_serie(i)+Q_serie(i-1) )
    baseflow(i) = c1+c2
    if (baseflow(i) .gt. Q_serie(i)) then
      baseflow(i) = Q_serie(i)
    end if
  end do


BFI = sum(baseflow) / sum(Q_serie)



END SUBROUTINE BaseFlowIndex

!-------------------------------------------------------------------------------
!         Quantile



SUBROUTINE Quantile(x,p,Q)


IMPLICIT NONE
	real*8, dimension(:),intent(in)		::x   ! input serie
	real*8, dimension(:),intent(in)		::p   ! probability array
	real*8, dimension(:),allocatable,intent(out)	::Q   ! quantiles
	real*8,dimension(:),allocatable		::m,g ! temp variables
	integer					::n   ! length of x
	integer,dimension(:),allocatable		::j   ! indexes


allocate(m(size(p)))
allocate(g(size(p)))
allocate(Q(size(p)))
allocate(j(size(p)))

n=size(x)

m=1-p
j=floor(n*p+m)
g=n*p+m-real(j,8)
Q=(1-g)*x(j)+g*x(j+1) 

deallocate(m)
deallocate(g)
deallocate(j)



END SUBROUTINE Quantile


!-------------------------------------------------------------------------------
!         getMonths


SUBROUTINE getMonths(dates_data, months_serie)


IMPLICIT NONE

 	character*10, dimension(:),intent(in)     :: dates_data       ! simulated  daily river
 	real*8,dimension(:), allocatable,intent(out)          :: months_serie       ! simulated 
	character*10			                :: temp_date
	character*10			                :: temp_month 
        real*8                                          :: month
        integer                                         :: ii

allocate(months_serie(size(dates_data,1)))

do ii=1, size(dates_data)
temp_date=dates_data(ii)
temp_month = temp_date(6:7)
read(temp_month,*) month
months_serie(ii) = month

end do

END SUBROUTINE getMonths


END MODULE mo_signatures
