MODULE mo_objectives

IMPLICIT NONE

CONTAINS


subroutine NashEfficiency(Qo,Qm, NSE)

		real*8,dimension(:), intent(in)	                :: Qo
		real*8,dimension(:), intent(in)	                :: Qm
		real*8						:: Q_av
		real*8, intent(out)	                        :: NSE
        real*8,dimension(:), allocatable                :: Qo_corr
        real*8,dimension(:), allocatable                :: Qm_corr

        Qo_corr = PACK(Qo, isnan(Qo) .eqv. .False. )
        Qm_corr = PACK(Qm, isnan(Qo) .eqv. .False. )

		Q_av=sum(Qo_corr)/size(Qo_corr)
		NSE = 1-sum( (Qm_corr - Qo_corr)**2) / sum( (Qo_corr - Q_av)**2)



end subroutine


!---------------------------------------------------------------------------------------

subroutine LogNashEfficiency(Qo,Qm, LNSE)

		real*8,dimension(:), intent(in)	                :: Qo
		real*8,dimension(:), intent(in)	                :: Qm
		real*8, intent(out)	                        :: LNSE

		real*8						:: Q_av
		real*8						:: Up
		real*8						:: Down
		real*8						:: SumQ

		real*8,dimension(size(Qo,1))                    :: logQo
		real*8,dimension(size(Qm,1))                    :: logQm
		integer						:: i
		integer						:: length, length2
        real*8,dimension(:), allocatable                :: Qo_corr
        real*8,dimension(:), allocatable                :: Qm_corr

        Qo_corr = PACK(Qo, isnan(Qo) .eqv. .False. )
        Qm_corr = PACK(Qm, isnan(Qo) .eqv. .False. )

do i = 1, size(Qo_corr,1)
   if(Qo_corr(i) .gt. tiny(Qo_corr)) then
      logQo(i)  = log(Qo_corr(i)) 
   else
      logQo(i)  = log(0.0001) 
   end if
end do


do i = 1, size(Qm_corr,1)
   if(Qm_corr(i) .gt. tiny(Qm_corr)) then
      logQm(i)  = log(Qm_corr(i)) 
   else
      logQm(i)  = log(0.0001) 
   end if
end do


Q_av=sum(logQo)/size(logQo)
LNSE = 1_8-sum( (logQm - logQo)**2) / sum( (logQo - Q_av)**2)


end subroutine

!---------------------------------------------------------------------------------------

subroutine LogNashEfficiency_corr(Qo,Qm, LNSE)

		real*8,dimension(:), intent(in)	                :: Qo
		real*8,dimension(:), intent(in)	                :: Qm
		real*8, intent(out)	                        :: LNSE

		real*8						:: Q_av
		real*8						:: Up
		real*8						:: Down
		real*8						:: SumQ
		real*8						:: n_obs

		real*8,dimension(size(Qo,1))                    :: logQo
		real*8,dimension(size(Qm,1))                    :: logQm
		integer						:: k
		integer						:: i
		integer						:: length, length2
        real*8,dimension(:), allocatable                :: Qo_corr
        real*8,dimension(:), allocatable                :: Qm_corr

Qo_corr = PACK(Qo, isnan(Qo) .eqv. .False. )
Qm_corr = PACK(Qm, isnan(Qo) .eqv. .False. )

k=1
do i = 1, size(Qo_corr,1)
   if(Qo_corr(i) .gt. tiny(Qo_corr)) then
      logQo(k)  = log(Qo_corr(k)) 
      logQm(k)  = log(Qm_corr(k)) 
      n_obs= real(k,8)
      k=k+1
   end if
  
end do




Q_av=sum(logQo(1:int(n_obs)))/n_obs
LNSE = 1_8-sum( (logQm - logQo)**2) / sum( (logQo - Q_av)**2)


end subroutine



!---------------------------------------------------------------------------------------

subroutine VolumeError(Qo,Qm, VE)

		real*8,dimension(:), intent(in)	                :: Qo
		real*8,dimension(:), intent(in)	                :: Qm
		real*8						:: Q_av
		real*8						:: Up
		real*8						:: Down
		real*8						:: SumQ
		real*8, intent(out)	                        :: VE
        real*8,dimension(:), allocatable                :: Qo_corr
        real*8,dimension(:), allocatable                :: Qm_corr

        Qo_corr = PACK(Qo, isnan(Qo) .eqv. .False. )
        Qm_corr = PACK(Qm, isnan(Qo) .eqv. .False. )

        VE=1- abs( sum(Qm_corr)-sum(Qo_corr) )/sum( Qo_corr ) 

end subroutine

!---------------------------------------------------------------------------------------

subroutine FDCNashEfficiency(Qo,Qm, FDC_NSE)

USE mo_sort

		real*8,dimension(:), intent(in)	                :: Qo
		real*8,dimension(:), intent(in)	                :: Qm
		real*8						:: Q_av
		integer, allocatable,dimension(:)               :: ind_sort_mod	! indexes
		integer, allocatable,dimension(:)               :: ind_sort_obs	! indexes
		real*8, allocatable,dimension(:)      		:: Q_sort_mod  	! sorted Q
		real*8, allocatable,dimension(:)      		:: Q_sort_obs  	! sorted Q
		integer						:: m,n
		real*8, intent(out)	                        :: FDC_NSE
        real*8,dimension(:), allocatable                :: Qo_corr
        real*8,dimension(:), allocatable                :: Qm_corr

        Qo_corr = PACK(Qo, isnan(Qo) .eqv. .False. )
        Qm_corr = PACK(Qm, isnan(Qo) .eqv. .False. )

allocate( ind_sort_mod(size(Qm_corr)) )
allocate( ind_sort_obs(size(Qo_corr)) )

allocate( Q_sort_mod(size(Qm_corr)) )
allocate( Q_sort_obs(size(Qo_corr)) )



Q_sort_obs=Qo_corr
Q_sort_mod=Qm_corr

call sort(Q_sort_obs)
call sort(Q_sort_mod)

n=size(Q_sort_obs)
Q_sort_obs=Q_sort_obs(n:1:-1)

m=size(Q_sort_mod)
Q_sort_mod=Q_sort_mod(m:1:-1)



call LogNashEfficiency(Q_sort_obs, Q_sort_mod, FDC_NSE)


end subroutine

!---------------------------------------------------------------------------------------

subroutine Rel_err( Qobs, Qm, F)
 
    IMPLICIT NONE

    real*8, intent(in)    				:: Qm, Qobs
    real*8, intent(out)                                	:: F


      F = 1-abs(1-(Qm/Qobs))

end subroutine
!---------------------------------------------------------------------------------------

subroutine logKlingGupta(Qo,Qm, LKGE)

		real*8,dimension(:), intent(in)	                :: Qo
		real*8,dimension(:), intent(in)	                :: Qm

		real*8, intent(out)	                        :: LKGE

                real*8                                          :: alpha
                real*8                                          :: beta
                real*8                                          :: cov
                real*8                                          :: cor
                real*8                                          :: ED
                integer                                         :: length
		real*8,dimension(size(Qo,1))                    :: logQo
		real*8,dimension(size(Qm,1))                    :: logQm
                integer                                         :: n
                integer                                         :: i
		real*8						:: Qo_mean
                real*8						:: Qm_mean
                real*8                                          :: sd_m
                real*8                                          :: sd_o
                real*8                                          :: sumQ
                real*8                                          :: tmp
                real*8                                          :: tmp2
                real*8                                          :: tmp3
                real*8,dimension(:), allocatable                :: Qo_corr
                real*8,dimension(:), allocatable                :: Qm_corr

Qo_corr = PACK(Qo, isnan(Qo) .eqv. .False. )
Qm_corr = PACK(Qm, isnan(Qo) .eqv. .False. )

do i = 1, size(Qo_corr,1)
   if(Qo_corr(i) .gt. tiny(Qo_corr)) then
      logQo(i)  = log(Qo_corr(i)) 
   else
      logQo(i)  = log(0.0001) 
   end if
end do


do i = 1, size(Qm_corr,1)
   if(Qm_corr(i) .gt. tiny(Qm_corr)) then
      logQm(i)  = log(Qm_corr(i)) 
   else
      logQm(i)  = log(0.0001) 
   end if
end do





Qo_mean=sum(logQo)/size(logQo,1)

Qm_mean=sum(logQm)/size(logQm,1)


n=size(logQo)

cov = (1/(real(n,8) - 1) ) * sum( (logQo-Qo_mean) * (logQm-Qm_mean))

sd_o = sqrt((1/(real(n,8) - 1) )* sum( (logQo-Qo_mean)**2))
sd_m = sqrt((1/(real(n,8) - 1) )* sum( (logQm-Qm_mean)**2))


cor  = cov/(sd_o*sd_m)

alpha= sd_m/sd_o

beta=Qm_mean/Qo_mean

ED=sqrt( (cor-1)**2 + (alpha-1)**2 + (beta-1)**2    )
  
LKGE=1-ED


end subroutine

!---------------------------------------------------------------------------------------

subroutine logKlingGupta_corr(Qo,Qm, LKGE)

		real*8,dimension(:), intent(in)	                :: Qo
		real*8,dimension(:), intent(in)	                :: Qm

		real*8, intent(out)	                        :: LKGE

                real*8                                          :: alpha
                real*8                                          :: beta
                real*8                                          :: cov
                real*8                                          :: cor
                real*8                                          :: ED
                integer                                         :: n_obs
		real*8,dimension(size(Qo,1))                    :: logQo
		real*8,dimension(size(Qm,1))                    :: logQm
                integer                                         :: k
                integer                                         :: n
                integer                                         :: i
		real*8						:: Qo_mean
                real*8						:: Qm_mean
                real*8                                          :: sd_m
                real*8                                          :: sd_o
                real*8                                          :: sumQ
                real*8                                          :: tmp
                real*8                                          :: tmp2
                real*8                                          :: tmp3
                real*8                                          :: length
                real*8,dimension(:), allocatable                :: Qo_corr
                real*8,dimension(:), allocatable                :: Qm_corr

Qo_corr = PACK(Qo, isnan(Qo) .eqv. .False. )
Qm_corr = PACK(Qm, isnan(Qo) .eqv. .False. )


k=1
do i = 1, size(Qo_corr,1)
   if(Qo(i) .gt. tiny(Qo_corr)) then
      logQo(k)  = log(Qo_corr(k)) 
      logQm(k)  = log(Qm_corr(k)) 
      n_obs= real(k,8)
      k=k+1
   end if
  
end do

!k=1
!do i = 1, size(Qo,1)
!   if(Qm(i) .gt. tiny(Qm)) then
!      logQm(k)  = log(Qm(k)) 
!      length2= real(k,8)
!   end if
      
!end do





Qo_mean=sum(logQo(1:int(n_obs)))/n_obs !size(logQo,1)

Qm_mean=sum(logQm(1:int(n_obs)))/n_obs !size(logQm,1)


n=size(logQo)

cov = (1/(real(n,8) - 1) ) * sum( (logQo-Qo_mean) * (logQm-Qm_mean))

sd_o = sqrt((1/(real(n,8) - 1) )* sum( (logQo-Qo_mean)**2))
sd_m = sqrt((1/(real(n,8) - 1) )* sum( (logQm-Qm_mean)**2))


cor  = cov/(sd_o*sd_m)

alpha= sd_m/sd_o

beta=Qm_mean/Qo_mean

ED=sqrt( (cor-1)**2 + (alpha-1)**2 + (beta-1)**2    )
  
LKGE=1-ED


end subroutine


subroutine KlingGupta(Qo,Qm, KGE)

		real*8,dimension(:), intent(in)	                :: Qo
		real*8,dimension(:), intent(in)	                :: Qm

		real*8, intent(out)	                        :: KGE

                real*8                                          :: alpha
                real*8                                          :: beta
                real*8                                          :: cov
                real*8                                          :: cor
                real*8                                          :: ED
                integer                                         :: n
		real*8						:: Qo_mean
                real*8						:: Qm_mean
                real*8                                          :: sd_m
                real*8                                          :: sd_o
                real*8,dimension(:), allocatable                :: Qo_corr
                real*8,dimension(:), allocatable                :: Qm_corr

Qo_corr = PACK(Qo, isnan(Qo) .eqv. .False. )
Qm_corr = PACK(Qm, isnan(Qo) .eqv. .False. )

Qo_mean=sum(Qo_corr )/size(Qo_corr)
Qm_mean=sum(Qm_corr)/size(Qm_corr)


n=size(Qo_corr)

cov = (1/(real(n,8) - 1) ) * sum( (Qo_corr-Qo_mean) * (Qm_corr-Qm_mean))

sd_o = sqrt((1/(real(n,8) - 1) )* sum( (Qo_corr-Qo_mean)**2))
sd_m = sqrt((1/(real(n,8) - 1) )* sum( (Qm_corr-Qm_mean)**2))


cor  = cov/(sd_o*sd_m)

alpha= sd_m/sd_o

beta=Qm_mean/Qo_mean

ED=sqrt( (cor-1)**2 + (alpha-1)**2 + (beta-1)**2    )
  
KGE=1-ED




end subroutine

subroutine WaterBalance(Qm, Ea, Ei, Prec, WB)

		real*8,dimension(:), intent(in)	                :: Qm
		real*8,dimension(:), intent(in)	                :: Ea
		real*8,dimension(:), intent(in)	                :: Ei
		real*8,dimension(:), intent(in)	                :: Prec

		real*8, intent(out)	                        :: WB

                real*8                                          :: EaP
                real*8                                          :: EiP
                real*8                                          :: RC
             


RC = sum(Qm)/sum(prec)
EaP = sum(Ea)/sum(prec)
EiP = sum(Ei)/sum(prec)

WB = EaP + EiP + RC



end subroutine


!---------------------------------------------------------------------------------------

subroutine snow_objective(Snow_obs,Snow_mod, snow_efficiency)

		real*8,dimension(:), intent(in)	                :: Snow_obs
		real*8,dimension(:), intent(in)	                :: Snow_mod

		real*8, intent(out)	                        :: snow_efficiency

                integer                                         :: i
                real*8                                          :: snowdays
		real*8,dimension(size(Snow_obs,1))              :: Snow_tmp


snowdays = 0_8

do i = 1, size(Snow_mod,1)

   if(  Snow_obs(i) .gt. tiny(Snow_obs)    ) then
   snowdays = snowdays + 1_8
   end if

   if( (Snow_mod(i) .gt. tiny(Snow_mod) ) .and.  ( Snow_obs(i) .gt. tiny(Snow_obs) )    ) then
      Snow_tmp(i)  = 1_8
   else
      Snow_tmp(i)  = 0_8
   end if
end do


snow_efficiency = 1_8 -  (snowdays - sum(Snow_tmp))/snowdays




end subroutine



END MODULE mo_objectives
