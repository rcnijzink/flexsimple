PROGRAM simple_model
! Program that runs the flex concept with optimization options

!Written by: R.C. Nijzink june 2014
!References: 
!Adjustments: 


USE mo_model
USE mo_readdata
USE mo_multi_obj
USE mo_objectives
USE mo_moving_window
USE mo_sample
USE mo_resample

IMPLICIT NONE


  
	integer                                 :: itsteps		         ! counter
	real*8, allocatable,dimension(:)	    :: etp_data,prec_data, Qobs_data,q, temp_data, sumax_data
	real*8, allocatable,dimension(:)	    :: imax_data
	real*8, allocatable,dimension(:)	    :: pmax_data
	real*8, allocatable,dimension(:)	    :: snow_data
	real*8, allocatable, dimension(:,:)	    :: output
	character*10, allocatable, dimension(:)	:: dates_data
	real*8,dimension(16)			        :: param
	real*8,dimension(16)		            :: param_max
	real*8,dimension(16)			        :: param_min
	logical,dimension(16)			        :: optim
	real*8,dimension(4)			            :: incon

        real*8, dimension(:,:), allocatable     :: dem
        real*8, dimension(:,:), allocatable     :: ndvi
        integer, dimension(:,:), allocatable    :: facc
        integer, dimension(:,:), allocatable    :: fdir
        real*8, dimension(:,:), allocatable     :: slope
        real*8                                  :: cellsize
        real*8                                  :: xllcorner
        real*8                                  :: yllcorner
        real*8                                  :: Area_thresh
        real*8                                  :: Slope_thresh
        real*8                                  :: Hand_thresh
        real*8                                  :: nodata
        integer                                 :: nCols
        integer                                 :: nRows
        integer                                 :: i, jj
        integer                                 :: fileunit


print *, "---------------------------------------"
print *, "               FLEX-simple               "
print *, "              version 0.1              "
print *, "             4 bucket model           "
print *, "---------------------------------------"



print *, ""
print *, "   Reading input data ..."
call read_config()
call read_forcing(etp_data,prec_data, Qobs_data, temp_data, dates_data, sumax_data, imax_data, pmax_data, snow_data)

if(snow_zones .eqv. .True.) then
call read_dem(dem, nRows, nCols, nodata, cellsize, xllcorner, yllcorner)
end if
call read_param(param, param_max, param_min, incon, optim)



! --run the model without calibration--
if(Optim_flag .eqv. .FALSE.) then

print *, ""
print *, "   Starting FLEXsimple ..."

call model(param(1:12),incon, prec_data, temp_data, etp_data, dem, cellsize, output)

print *, ""
print *, "   Saving output ..."
do i=1, 11

fileunit=1


open(fileunit, file=trim(adjustl(output_dir)) // trim(fluxes_states_names(i)), status='unknown', action='write')
	do jj=1, size(output(1,:)) 
  	write(fileunit,'((f6.2))') output(i,jj)
	end do
 close(fileunit)



end do

print *, ""
print *, "   FLEXsimple finished!"


! --Start optimization--
else

select case(Method)

!if( Method .eq. 1) then
case(1)
print *, ""
print *, "Start Monte Carlo optimization ..."

call multi_obj(param,param_max,param_min,incon, prec_data,temp_data, etp_data,&
               dem, cellsize, Qobs_data, dates_data, optim, sumax_data)


case(2)
print *, ""
print *, "   Start moving window calibration with random values ..."
call moving_window(param,param_max,param_min,incon, prec_data,temp_data,etp_data,&
                 dem, cellsize, Qobs_data, dates_data, optim)

case(3)
print *, ""
print *, "   Start random run calibration ..."
call sample(param,param_max,param_min,incon, prec_data,temp_data, etp_data,&
               dem, cellsize, Qobs_data, dates_data, optim, sumax_data, imax_data)

case(4)
print *, ""
print *, "   Start calibration with resampled sets ..."
call resample(param,param_max,param_min,incon, prec_data,temp_data, etp_data,&
               dem, cellsize, Qobs_data, dates_data, optim, sumax_data)


end select




print *, ""
print *, "  FLEXsimple Finished!"
end if






END PROGRAM
