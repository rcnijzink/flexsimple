MODULE mo_model

IMPLICIT NONE

CONTAINS

      subroutine model(param, incon, prec, airt, ep, dem, cellsize, output, sumax_serie,Imax_serie)

USE mo_readdata

       	real*8, dimension(11), intent(in)	             :: param	! model parameters
      	real*8, dimension(4),intent(in)		             :: incon	! initial states	
      	real*8,dimension(:), intent(in)		             :: prec	! precipation series
	real*8,dimension(:), intent(in)		             :: airt	! temperature series
	real*8,dimension(:), intent(in)		             :: ep	! evaporation series
        real*8, dimension(:,:), allocatable, intent(in)      :: dem
        real*8,intent(in)                                    :: cellsize
      	real*8,allocatable, dimension(:,:), intent(out)	     :: output ! Output
	real*8,dimension(:), intent(in), optional            :: sumax_serie ! sumax series
	real*8,dimension(:), intent(in), optional            :: Imax_serie ! sumax series
	!real*8,dimension(:), intent(in), optional            :: Mmelt_serie ! sumax series


        integer 				:: tmax	        ! number of timesteps
	integer 				:: i		! counter
	integer 				:: it		! counter of timestep

     !forcing
     real*8 					:: temp		! temperature at timestep it
     real*8 					:: Pdt  	! precipitation at timestep it
     real*8 					:: dt	        ! water balance

     !model parameters
     real*8					::  Meltfactor  
     real*8					::  Tthresh
     real*8					::  Imax 
     real*8					::  Sumax        
     real*8					::  beta 
     real*8					::  Kf                      
     real*8					::  Ks               
     real*8					::  LP        
     real*8					::  D
     real*8					::  Pmax
     real*8					::  alpha

    !states
     real*8					:: Si
     real*8					:: Su
     real*8					:: Su_accent
     real*8					:: Sf
     real*8					:: Ss
     real*8					:: Ss_in
     real*8					:: Ss_corr

    !fluxes
     real*8                  		        :: Pedt	! precipation series
     real*8 					:: R	! runoff at timestep it
     real*8, dimension(:), allocatable	        :: Rs 	! rain at timestep it
     real*8, dimension(:), allocatable	        :: Rf 	! rain at timestep it
     real*8, dimension(:), allocatable	        :: Peff 	! rain at timestep it
     real*8, dimension(:), allocatable	        :: SnowRest 	! rain at timestep it
     real*8 					:: Kt 	        ! rain at timestep it
     real*8 					:: Pr 	        ! rain at timestep it      
     real*8 					:: Pr_su        ! rain at timestep it 
     real*8 					:: Pedt_active 	! Potential evaporation     
     real*8, dimension(:), allocatable	        :: Weights_f 	! rain at timestep it
     real*8 					:: Eidt         ! interception evaporation
     real*8 					:: Eadt	        ! actual evaporation
     real*8 					:: Epdt 	! Potential evaporation
     real*8 					:: Eipdt	! degree day factor

     real*8 					:: Qf	! degree day factor
     real*8 					:: Qr	! degree day factor
     real*8 					:: Qs	! degree day factor
     real*8 					:: Qs_tot	! degree day factor
     real*8 					:: Qm	! degree day factor
     real*8 					:: Qsum	! degree day factor

    !water balance
     real*8 					:: Flux_out	! sum of out
     real*8 					:: Flux_in	! sum of in
     real*8 					:: Storage_in	! storage t=1
     real*8 					:: Storage_out	! storage last t
     real*8 					:: WB	        ! water balance



!Initialize model parameters

tmax=size(prec)
allocate( Rs( tmax) )
allocate( Rf( tmax) )
allocate( output(12,tmax) )

      
  Meltfactor  = param(1) 
  Tthresh     = param(2)  
  Imax        = param(3)
  Sumax       = param(4)
  beta        = param(5)
  Kf          = param(6)
  Ks          = param(7)
  LP          = param(8)
  D           = param(9)
  Pmax        = param(10)
  alpha       = param(11)


!Sumax     = sumax_loc/ (1+beta)  !maximum storage in catchment equals (1+beta) times maximum occuring value of soil moisture capacity

!Initialize the states of the model

  
  Si = incon(1)
  Su = incon(2)
  Ss = incon(3)
  Sf = incon(4)
Su_accent = zero_dp
Ss_corr = zero_dp
Rs = zero_dp
Pr = zero_dp

  !call snow module
   if(snow_flag .eqv. .TRUE.) then
        call snowmod(airt, prec, dem, cellsize, Tthresh, Meltfactor, Peff, SnowRest, .TRUE.) 
    else 
   Peff = prec
   end if

dt =1_8

do it=1,tmax
       
       Pdt =Peff(it)
       temp=airt(it)
       Epdt=ep(it)

       if(present(sumax_serie)) then
          Sumax = sumax_serie(it)

       end if

       if(present(Imax_serie)) then
          Imax = Imax_serie(it)
       end if

       !if(present(Mmelt_serie)) then
       !   Meltfactor = Mmelt_serie(it)
       !end if

       if(Pdt .le. epsilon(zero_dp)) then
       Pdt = zero_dp
       end if


!---------------------------------------------------
! simple 
!---------------------------------------------------

! Interception reservoir  

  Si=Si+Pdt 

  Pedt=maxval( (/zero_dp, (Si-Imax) /) )
  Si=Si-Pedt

  Eidt=minval( (/Epdt   ,Si/) )
  Si=Si-Eidt
  Eipdt=maxval( (/zero_dp, (Epdt-Eidt) /) )


   if((Pedt) .gt. zero_dp) then

      if( (Pedt +Su_accent) .lt. ((1+beta)*Sumax) ) then
         R  = Pedt  - Sumax + Su + Sumax * (1- (Pedt  + Su_accent)/( (1+beta)*Sumax ) )**(1+beta) !R=P-E-Wm+W+Wm[1-(P-E+a)/W'mm]^(1+b)
      end if

      if( (Pedt +Su_accent) .ge. ((1+beta)*Sumax) ) then
         R = Pedt  - (Sumax - Su)
      end if

   else

      R = zero_dp

   end if

   Su = Su + (Pedt  - R)



!actual evaporation/transpiration
  if(Su .lt. (Sumax*LP)) then
    Kt=Su/(Sumax*LP)
  else
    Kt=1
  end if 


  Eadt=Eipdt*Kt
  Eadt=minval( (/Eadt, Su/) )
  Su = Su - Eadt



!percolation from soil moisture to groundwater
   Pr = Pmax * ( Su/Sumax)
   Pr = minval((/Su, Pr/))
   Su = Su - Pr

   Su_accent = ((1+beta)*Sumax )*(1- (1-Su/Sumax)**(1/(1+beta)))


! Fast Reservoir -  non-linear
  Rf(it)=(1-D)*R

  Sf = Sf + Rf(it)
  Qf = (Sf**alpha) / Kf
  Qf = minval( (/Qf,Sf/) )
  Sf = Sf - Qf

! Recharge to slow reservoir
  Rs(it)=D*R


!---------------------------------------------------
! Slow (shared) reservoir 
!---------------------------------------------------

  Ss = Ss + Rs(it) + Pr


  !Qs = zero_dp
  !Qs = (Ss**Nss) / Ks 
  Qs = Ss / Ks 
  Qs = minval( (/Qs,Ss/) )

  Ss = Ss - Qs



!  Ss_in = Ss + Rs(it) + Pr
!   Ss = (Ss_in*exp(- 1_8/Ks)) - (  (L*Ks) * ( 1_8-exp(- 1_8/Ks) ) )

!if(Ss .gt. zero_dp) then
!   Qs_tot = Ss_in - Ss 
!   Qs = maxval( (/zero_dp, (Qs_tot - L)   /) )
!else
!   Qs_tot = L
!   Ss = Ss_in - L
!   Qs = zero_dp
!end if



!---------------------------------------------------
! Results
!---------------------------------------------------

 Qm = Qs  + Qf



 !write fluxes to output
       output(1,it)  = Qm 
       output(2,it)  = Qs
       output(3,it)  = Pr !Qf
       output(4,it)  = Peff(it)
       output(5,it)  = Eadt
       output(6,it)  = Eidt
       output(7,it)  = Rs(it)

!write states to output

       output(8,it)=Si
       output(9,it)=Su
       output(10,it)=Sf
       output(11,it)=Ss
       output(12,it)=SnowRest(it)

end do


! check waterbalances
! WB 
      Flux_out   =     sum(output(1,:))  + &
                       sum(output(5,:))  + & 
                       sum(output(6,:))  

      Flux_in    =     sum(Peff)

      Storage_in =    incon(1)  +&
                      incon(2)  +&
                      incon(3)  +&
                      incon(4) 

      Storage_out =  output(8,tmax) +&
                     output(9,tmax) +&
                     output(10,tmax)+&
                     output(11, tmax) 
                     

      WB = Flux_in - Flux_out - (Storage_out-Storage_in) !- L*real(tmax,8) 
      
if( (WB .gt. 0.00001) .or. (WB .lt. -0.00001) ) then
print *, "WARNING: waterbalance not closed"
print *, "Waterbalance:                  " , WB
print *, Flux_out, Flux_in, Storage_in, Storage_out, Sumax
stop
end if


      end subroutine model

!------------------------------------------------------------------------------------------------

      subroutine weighfun_triangle(Tlag, Weights)

       real*8, intent(in) 				:: Tlag
       real*8, dimension(:), allocatable, intent(out)	:: Weights 
       integer                                          :: i, j
       integer                                          :: nlag


nlag = ceiling(Tlag)

allocate( Weights(nlag) )


do j=1, nlag
Weights(j) = real(j,8) / real(sum( (/(i,i=1,nlag)/)),8)
end do


      end subroutine weighfun_triangle


!------------------------------------------------------------------------------------------------

      subroutine snowmod(Temp, Prec, dem, cellsize, Tthresh, Meltfactor, Peff, SnowRest, zones)

use mo_readdata

        real*8, dimension(:), intent(in)                    :: Temp
        real*8, dimension(:), intent(in)                    :: Prec
        real*8, dimension(:,:), intent(in)                  :: dem
        real*8, intent(in)                                  :: Tthresh                             
        real*8, intent(in)                                  :: Meltfactor                             
        real*8, intent(in)                                  :: cellsize
        logical, intent(in)                                 :: zones 
        integer                                             :: station_zone
        real*8, dimension(:), allocatable                   :: CatElev
        real*8, dimension(:), allocatable                   :: Area
        real*8, dimension(:), allocatable                   :: Weights
        real*8, dimension(:), allocatable                   :: h
        real*8, dimension(:), allocatable                   :: Elevation
        real*8                                              :: LapseRate 
        real*8                                              :: MinElev
        real*8                                              :: MaxElev
        integer                                             :: i, j, n
        integer                                             :: nElev
        integer                                             :: Area_cat
        real*8, dimension(:,:), allocatable                 :: SnowEq
        real*8, dimension(:,:), allocatable                 :: Melt
        real*8, dimension(:,:), allocatable                 :: TempMat
        real*8, dimension(:,:), allocatable                 :: Pmat 
        real*8, dimension(:), allocatable, intent(out)      :: Peff
        real*8, dimension(:), allocatable, intent(out)      :: SnowRest

LapseRate=-0.0064

if(zones .eqv. .TRUE.) then



!determine weights and height zones

!1D array of DEM
Elevation = pack(dem, dem /= -9999)

MinElev=minval((/Elevation,Elev_station/))
MaxElev=maxval((/Elevation,Elev_station/))


!number of zones, 100m elevation difference per zone
nElev=(ceiling(MaxElev/real(100,8))-floor(MinElev/real(100,8)))



!find the zone where the meteostation is located
if( Elev_station .gt. MinElev ) then
station_zone=ceiling(Elev_station/real(100,8))-floor(MinElev/real(100,8))
else
station_zone = 1 
end if

!create elevation categories
allocate( CatElev( nElev+1) )


CatElev(1) = floor(MinElev/real(100,8))*real(100,8)
do i=2, nElev+1
CatElev(i) = CatElev(i-1) + real(100,8)
end do



!determine the weights for each zone
allocate( Area( nElev) )
allocate( Weights( nElev) )

do j=1, nElev
Area(j)=0
do i=1, size(Elevation)  
  if( (Elevation(i) .gt. CatElev(j)) .and. (Elevation(i) .le. CatElev (j+1)) ) then
  Area(j) = Area(j) + cellsize**2
  end if
end do
end do

Weights=Area/sum(Area)



deallocate( Area )
allocate( h( nElev) )
h(:) = 0


!determine relative height from meteo station
do i=Station_zone, 1, -1
  if(i==station_zone) then
    h(i) = 0
  else
    h(i) = h(i+1) - real(100,8)
  end if
end do

do j = (station_zone+1) , nElev
      h(j)=h(j-1)+ real(100,8)
end do



else

nElev = 1
allocate( h( nElev) )
allocate( Weights( nElev) )

h = 0
Weights = 1


end if


 allocate( SnowEq( size(Temp),nElev  ) )
 allocate( Melt(size(Temp),nElev ) )
 allocate( TempMat( size(Temp),nElev ) )
 allocate( Pmat(size(Temp),nElev ) )

do n = 1, nElev
TempMat(:,n) = Temp(:) + h(n)*LapseRate
end do



SnowEq(:,:) = 0
Melt(:,:) = 0

do j = 1, nElev
  do i = 1, size(Temp)
     if( TempMat(i,j) .lt. Tthresh) then
        if(i .gt. 1) then
               SnowEq(i,j) = Prec(i) + SnowEq(i-1,j)
               Pmat(i,j) = 0
        else
               SnowEq(i,j) = Prec(i)
               Pmat(i,j)=0
        end if   
      else
        if(i .gt. 1) then
               Melt(i,j)    = Meltfactor * (TempMat(i,j) - Tthresh)
               Melt(i,j)    = minval( (/Melt(i,j) , SnowEq(i-1,j)/) )
               SnowEq(i,j)  = SnowEq(i-1,j) - Melt(i,j)
               Pmat(i,j)    = Prec(i)
        else
               SnowEq(i,j)  = 0
               Pmat(i,j)    = Prec(i)
        end if
      end if
  end do
end do




allocate( Peff     (size(temp) ))
allocate( SnowRest (size(temp) ))

Peff(:)      = 0
SnowRest(:)  = 0

do i = 1, nElev
Peff(:)        = Peff(:) + Pmat(:,i) * Weights(i) + Melt(:,i) * Weights(i)
SnowRest(:)    = SnowRest(:) + SnowEq(:,i) * Weights(i)
end do



 deallocate( SnowEq )
 deallocate( Melt )
 deallocate( TempMat )
 deallocate( Pmat )







      end subroutine snowmod

END MODULE mo_model