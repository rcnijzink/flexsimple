MODULE mo_readdata

IMPLICIT NONE


	logical,public				::Optim_flag
	logical,public				::snow_flag
	logical,public				::snow_zones
	logical,public				::dyn_mode
        logical,public	                        ::read_data
        logical,public	                        ::optim_growth
        character*10,public                     ::change_start
        character*10,public                     ::change_end

	real*8,public				::mm_low
	real*8,public				::mm_high
	real*8,public				::Elev_station
	integer,public				::Method
    character*100               ::param_file    
	integer,public				::Iterations
	integer,public				::window
	integer,public				::N_samples
	integer,public				::maxout
	integer,public				::warmup
	character*10,public			::startdate
	character*10,public			::enddate
	character*10,public			::cal_start
	character*10,public			::val_start
	character*10,public			::cal_end
	character*10,public			::val_end
	character*10,public			::warmup_start
	character*10,public			::warmup_end
	integer,public				:: ichange_start
	integer,public				:: ichange_end
	integer,public				:: ic_start
	integer,public				:: ic_end
	integer,public				:: iv_start
	integer, public				:: iv_end
	integer,public				:: iw_start
	integer, public				:: iw_end
	character*30,public			::input_dem
	character*30,public			::input_facc
	character*30,public			::input_fdir
	character*30,public			::input_ndvi
	character*30,public			::input_slope
	character*150,public			::input_dir
	character*150,public			::output_dir
	character*150,public			::output_dir_cal
	character*150,public			::output_dir_val
	character*20,dimension(15),public	::fluxes_states_names
        real*8, public                          ::zero_dp = real(0,8)






CONTAINS

      	subroutine read_forcing(etp_data,prec_data, Qobs_data, temp_data, dates_data, sumax_data,Imax_data,Pmax_data, Snow_data)

  	IMPLICIT NONE

  	character                                               :: filename	        ! name of input file
    	integer                                                 :: fileunit1	        ! file unit
    	integer                                                 :: fileunit2	        ! file unit
    	integer                                                 :: fileunit3	        ! file unit
    	integer                                                 :: fileunit4	        ! file unit
	integer				                        :: fileunit5		! file unit
	integer				                        :: fileunit6		! file unit
	integer				                        :: fileunit7		! file unit
	integer				                        :: fileunit8		! file unit
	integer				                        :: fileunit9		! file unit
	integer                                                 :: i	         	! counter
	integer				                        :: length		! length of modelling period
	real*8,dimension(:),allocatable,intent(out)             :: etp_data,prec_data, Qobs_data, temp_data, sumax_data
	real*8,dimension(:),allocatable,intent(out)             :: Imax_data
	real*8,dimension(:),allocatable,intent(out)             :: Pmax_data
	real*8,dimension(:),allocatable,intent(out)             :: Snow_data
	character*10,dimension(:),allocatable,intent(out)	:: dates_data
	character*10			:: simdate


simdate=startdate

fileunit1=61
fileunit2=62
fileunit3=63
fileunit4=64
fileunit5=65
fileunit6=66
fileunit7=67
fileunit8=68
fileunit8=69

open(unit=fileunit1, file=trim(adjustl(input_dir)) // trim(adjustl("Date.txt")), action='read')

! find indexes of calibration start and end and determine the total length to be modeled.
length=0
	do while (simdate .ne. enddate )
	 	read(fileunit1, *) simdate
		length=length+1
		if (simdate .eq. cal_start) then
		ic_start=length
		end if
		if (simdate .eq. cal_end) then
		ic_end=length
		end if
		if (simdate .eq. val_start) then
		iv_start=length
		end if
		if (simdate .eq. val_end) then
		iv_end=length
		end if
		if (simdate .eq. warmup_start) then
		iw_start=length
		end if
		if (simdate .eq. warmup_end) then
		iw_end=length
		end if
                if (simdate .eq. change_start) then
		ichange_start=length
		end if
                if (simdate .eq. change_end) then
		ichange_end=length
		end if

         end do
 	
 close(fileunit1)

open(unit=fileunit1, file=trim(adjustl(input_dir)) // trim(adjustl("Date.txt")), action='read')
open(unit=fileunit2, file=trim(adjustl(input_dir)) // trim(adjustl("Etp.txt")), action='read', status='old')
open(unit=fileunit3, file=trim(adjustl(input_dir)) // trim(adjustl("Prec.txt")), action='read', status='old')
open(unit=fileunit4, file=trim(adjustl(input_dir)) // trim(adjustl("Qobs.txt")), action='read', status='old')

if(snow_flag .eqv. .TRUE.) then
open(unit=fileunit5, file=trim(adjustl(input_dir)) // trim(adjustl("Temp.txt")), action='read', status='old')
!open(unit=fileunit9, file=trim(adjustl(input_dir)) // trim(adjustl("Swe.txt")), action='read', status='old')
end if

if(dyn_mode .eqv. .TRUE.) then
open(unit=fileunit6, file=trim(adjustl(input_dir)) // trim(adjustl("Sumax.txt")), action='read', status='old')
open(unit=fileunit7, file=trim(adjustl(input_dir)) // trim(adjustl("Imax.txt")), action='read', status='old')
open(unit=fileunit8, file=trim(adjustl(input_dir)) // trim(adjustl("Pmax.txt")), action='read', status='old')
end if



allocate( etp_data  ( length) )
allocate( prec_data ( length) )
allocate( Qobs_data ( length) )
allocate( temp_data ( length) )
allocate( dates_data     ( length) )
allocate( sumax_data     ( length) )
allocate( Imax_data     ( length) )
allocate( Pmax_data     ( length) )
allocate( Snow_data     ( length) )

	do i=1, length
                read(fileunit1, *) dates_data(i)
		read(fileunit2, *) etp_data(i)
		read(fileunit3, *) prec_data(i)
		read(fileunit4, *) Qobs_data(i)

		if(snow_flag .eqv. .TRUE.) then
		   read(fileunit5, *) temp_data(i)
                 !  read(fileunit9, *) Snow_data(i)
		else
		   temp_data(i)=0
		end if

                if(dyn_mode .eqv. .TRUE.) then
		   read(fileunit6, *) sumax_data(i)
		   read(fileunit7, *) Imax_data(i)
		   read(fileunit8, *) Pmax_data(i)
		else
		   sumax_data(i)=0
		   Imax_data(i)=0
		   Pmax_data(i)=0
		end if

	end do

	 close(fileunit1)
	 close(fileunit2)
	 close(fileunit3)
 	 close(fileunit4)
		if(snow_flag .eqv. .TRUE.) then
	 	close(fileunit5)
	 !	close(fileunit9)
		end if

                if(dyn_mode .eqv. .TRUE.) then
	 	close(fileunit6)
	 	close(fileunit7)
	 	close(fileunit8)
		end if

	end subroutine

!---------------------------------------------------------------------------------------------

subroutine read_param(param, param_max, param_min, incon, optim)

  IMPLICIT NONE


	real*8					                ::  Si   
	real*8					                ::  Su
	real*8					                ::  Sf
	real*8					                ::  Ss

     real*8,dimension(4)					::  Meltfactor  
     real*8,dimension(4)					::  Tthresh    
     real*8,dimension(4)					::  Kf       
     real*8,dimension(4)					::  Ks        
     real*8,dimension(4)					::  LP   
     real*8,dimension(4)					::  Imax      
     real*8,dimension(4)					::  Sumax       
     real*8,dimension(4)					::  beta      
     real*8,dimension(4)					::  Nlagf     
     real*8,dimension(4)					::  D  
     real*8,dimension(4)					::  Pmax  
     real*8,dimension(4)					::  a
     real*8,dimension(4)					::  alpha
     real*8,dimension(4)					::  b
     real*8,dimension(4)					::  sumax_min
     real*8,dimension(7)					::  tmp 
     integer                                                    ::  i


	integer				                        :: fileunit6		! file unit

	real*8,dimension(14), intent(out)			:: param
	real*8,dimension(14), intent(out)		        :: param_max
	real*8,dimension(14), intent(out)			:: param_min
	logical,dimension(14), intent(out)			:: optim
	real*8,dimension(4), intent(out)			:: incon


namelist /parameters/   Meltfactor, Tthresh, Imax, Sumax, beta, &
     Kf, Ks, LP, D, Pmax, alpha, a,b, sumax_min

namelist /ini_states/   Si, Su, Sf, Ss
   

optim = .FALSE.

       open(90,file=param_file, delim='apostrophe')
       read(90,nml=parameters)

       read(90,nml=ini_states)
       close(90)


param=(/Meltfactor(3), Tthresh(3) , &
       Imax(3),  Sumax(3),  beta(3), Kf(3), Ks(3), LP(3), D(3), Pmax(3), alpha(3), a(3),b(3),sumax_min(3)   &
       /)

param_max=(/Meltfactor(2), Tthresh(2) , &
       Imax(2),  Sumax(2), beta(2), Kf(2), Ks(2), LP(2), D(2), Pmax(2), alpha(2), a(2), b(2),sumax_min(2) &
       /)

param_min=(/Meltfactor(1), Tthresh(1) , &
       Imax(1),  Sumax(1), beta(1), Kf(1), Ks(1), LP(1), D(1), Pmax(1), alpha(1), a(1), b(1),sumax_min(1) &
       /)

where( (/Meltfactor(4), Tthresh(4) , &
       Imax(4),  Sumax(4), beta(4), Kf(4), Ks(4), LP(4), D(4), Pmax(4), alpha(4), a(4), b(4),sumax_min(4) &
       /) .eq. 1) optim = .TRUE.

incon=(/Si, Su, Sf, Ss/)


       fluxes_states_names(1)  = "Qm" 
       fluxes_states_names(2)  = "Qs"
       fluxes_states_names(3)  = "Qf"

       fluxes_states_names(4)  = "Pedt"   
       fluxes_states_names(5)  = "Eadt"
       fluxes_states_names(6)  = "Eidt"

       fluxes_states_names(7) = "Rs"

       fluxes_states_names(8)="Si"
       fluxes_states_names(9)="Su"
       fluxes_states_names(10)="Sf"
       fluxes_states_names(11)="Ss"



fileunit6=61







end subroutine

!-------------------------------------------------------------------------------------

subroutine read_config()

  IMPLICIT NONE

        
     CHARACTER(len=100), DIMENSION(:), ALLOCATABLE :: arguments      ! array for arguments
     INTEGER                                       :: iargs          ! Counter
     INTEGER                                       :: num_args       ! Number of arguments
     CHARACTER*100                                 :: input_dir_tmp  ! Temporary inputpath 
     CHARACTER*100                                 :: output_dir_tmp ! Temporary inputpath 
     CHARACTER*100                                 :: param_file_tmp ! Temporary parameterfile 
     CHARACTER*100                                 :: config_file_tmp! Temporary configfile
     CHARACTER*100                                 :: config_file    ! configfile
     LOGICAL                                       :: change_in      ! flag to change input
     LOGICAL                                       :: change_out     ! flag to change output
     LOGICAL                                       :: change_param   ! flag to change parameterfile
     LOGICAL                                       :: change_config  ! flag to change config

    namelist /general/ snow_flag, snow_zones, Elev_station, mm_low, mm_high,                                                    &
                   startdate, enddate, warmup_start, warmup_end, cal_start, cal_end, val_start, val_end,        & 
         input_dem, input_dir, output_dir, output_dir_cal, output_dir_val
    namelist /optimization/ Optim_flag, Method, Iterations, window, dyn_mode, read_data, optim_growth, &
         change_start, change_end

     !count arguments
     num_args = command_argument_count()

     !create array to save arguments
     allocate(arguments(num_args))  

     do iargs = 1, num_args
         !loop over arguments and save them
         call get_command_argument(iargs,arguments(iargs))
     end do

     change_in  = .FALSE.
     change_out = .FALSE.
     change_config = .FALSE.

     !loop over saved arguments and check flags
     do iargs = 1, num_args

        !check inputpath
        if(arguments(iargs) .eq. "-i") then
           input_dir_tmp = arguments(iargs+1)
           change_in = .True.
        end if

        !check outputpath
        if(arguments(iargs) .eq. "-o") then
           output_dir_tmp = arguments(iargs+1)
           change_out = .True.
        end if

        !check namelist and read again if needed
        if(arguments(iargs) .eq. "-c") then
           config_file_tmp = arguments(iargs+1)
           write(*,*) "Changed config.nml:", config_file_tmp
           change_config = .True.
        end if

        if(arguments(iargs) .eq. "-p") then
           param_file_tmp = arguments(iargs+1)
           write(*,*) "Changed param.nml:", param_file_tmp
           change_param = .True.
        end if

     end do

    if(change_config .eqv. .TRUE.) then
        config_file = config_file_tmp
    else
        config_file = 'config.nml'
    end if

    if(change_param .eqv. .TRUE.) then
          param_file = param_file_tmp
        else
          param_file = 'param.nml'
    end if


       open(11,file=config_file)
       read(11,nml=general)
       read(11,nml=optimization)
       close(11)

       if(change_in .eqv. .True.) then
          input_dir = input_dir_tmp
          write(*,*) "Changed input to:", input_dir
       end if 

       if(change_out .eqv. .True.) then
          output_dir = output_dir_tmp
          output_dir_cal = trim(adjustl(output_dir))//trim(adjustl("/cal/")) 
          output_dir_val = trim(adjustl(output_dir))//trim(adjustl("/val/")) 
          write(*,*) "Changed output to:", output_dir
       end if 




end subroutine

!-------------------------------------------------------------------------------------
subroutine read_dem(dem, nRows, nCols, nodata, cellsize, xllcorner, yllcorner)

!read_dem(dem, facc, fdir, slope, NDVI, nRows, nCols, nodata, cellsize, xllcorner, yllcorner)

  	IMPLICIT NONE

    	integer                         :: fileunit1	        ! file unit
    	integer                         :: fileunit2	        ! file unit
    	integer                         :: fileunit3	        ! file unit
    	integer                         :: fileunit4	        ! file unit
    	integer                         :: fileunit5	        ! file unit
	integer                         :: i,j	         	! counter
        integer,intent(out)           :: nRows     ! number of rows of data fields: 
        integer,intent(out)           :: nCols     ! number of columns of data fields: 
        real*8,intent(out)            :: xllcorner ! header read in lower left corner
        real*8,intent(out)            :: yllcorner ! header read in lower left corner
        real*8,intent(out)                          :: cellsize  ! header read in cellsize
        real*8,intent(out)                          :: nodata    ! header read in nodata value
	real*8,dimension(:,:),allocatable,intent(out):: dem
	!integer,dimension(:,:),allocatable,intent(out):: facc
	!integer, dimension(:,:),allocatable,intent(out):: fdir
	!real*8,dimension(:,:),allocatable,intent(out):: slope
	!real*8,dimension(:,:),allocatable,intent(out):: NDVI
        character                     :: dummy


fileunit1=46
fileunit2=47
fileunit3=48
fileunit4=49
fileunit5=50


 open(unit=fileunit1, file=trim(adjustl(input_dir)) // trim(adjustl(input_dem))  , action='read')

    read (fileunit1, *) dummy, nCols
    read (fileunit1, *) dummy, nRows
    read (fileunit1, *) dummy, xllcorner
    read (fileunit1, *) dummy, yllcorner
    read (fileunit1, *) dummy, cellsize
    read (fileunit1, *) dummy, nodata
      close(fileunit1)

allocate( dem  ( nRows, nCols) )

  open (unit=fileunit1, file=trim(adjustl(input_dir)) // trim(adjustl(input_dem))  , action='read')
    ! (a) skip header
    do i=1,6
       read(fileunit1, *)
    end do
    ! (b) read data
    do i=1,nRows
       read(fileunit1, *) dem(i,:)
    end do
    close(fileunit1)

  
!allocate( facc  ( nRows, nCols) )

!  open (unit=fileunit2, file=trim(adjustl(input_dir)) // trim(adjustl(input_facc))  , action='read')
    ! (a) skip header
!    do i=1,6
!       read(fileunit2, *)
!    end do
    ! (b) read data
!    do i=1,nRows
!       read(fileunit2, *) facc(i,:)
!    end do
!   close(fileunit2)

!allocate( fdir  ( nRows, nCols) )

!  open (unit=fileunit3, file=trim(adjustl(input_dir)) // trim(adjustl(input_fdir))  , action='read')
!    ! (a) skip header
!    do i=1,6
!       read(fileunit3, *)
!    end do
!    ! (b) read data
!    do i=1,nRows
!       read(fileunit3, *) fdir(i,:)
!    end do
!    close(fileunit3)


   
!allocate( slope  ( nRows, nCols) )

!  open (unit=fileunit4, file=trim(adjustl(input_dir)) // trim(adjustl(input_slope))  , action='read')
!    ! (a) skip header
!    do i=1,6
!       read(fileunit4, *)
!    end do
!    ! (b) read data
!    do i=1,nRows
!       read(fileunit4, *) slope(i,:)
!    end do
!    close(fileunit4)

!allocate( NDVI  ( nRows, nCols) )

!  open (unit=fileunit1, file=trim(adjustl(input_dir)) // trim(adjustl(input_NDVI))  , action='read')
!    ! (a) skip header
!    do i=1,6
!       read(fileunit5, *)
!    end do
!    ! (b) read data
!    do i=1,nRows
!       read(fileunit5, *) NDVI(i,:)
!    end do
!    close(fileunit5)


end subroutine



!---------------------------------------------------------------------------------------------

subroutine gen_sumax(param,sumax, i_change_start, i_change_end, prec_serie, sumax_serie )

  IMPLICIT NONE


	real*8, dimension(3), intent(in)                        ::  param
	real*8, intent(in)			                ::  sumax
	real*8, dimension(:),intent(in)		                ::  prec_serie
	integer					                ::  i_change_start
	integer					                ::  i_change_end
	integer					                ::  k
	integer					                ::  j
	integer					                ::  i
	real*8, dimension(:), intent(inout)	                ::  sumax_serie 	! sumax data
	real*8					                ::  dSu
	real*8					                ::  Sumax_min
	real*8					                ::  a
	real*8					                ::  b
	real*8, dimension(:),allocatable	        	::  t


!allocate(sumax_serie (size(prec_serie) ) )

!initialize with equilibrium value
sumax_serie = sumax

sumax_min   = param(3)
a           = param(1)
b           = param(2)


!linear reduction during change
dSu = (sumax - sumax_min)   / real(i_change_end- i_change_start,8)


k=1_8



do j=(i_change_start+1), i_change_end
sumax_serie(j) = sumax- (dSu*real(k,8))
k=k+1

end do






allocate(t  (size(sumax_serie)-i_change_end  ) )

t(1: (size(sumax_serie)-i_change_end)   ) = (/  (real(i,8), i= 1,(size(sumax_serie)-i_change_end)  ) /)




sumax_serie( (i_change_end+1) :size(sumax_serie) )=  ((sumax - sumax_min) * (1-exp(-a*(t**b)) ) )+ sumax_min

!write test file
 open(160, file=trim(adjustl(output_dir_cal)) // trim(adjustl("Imax_serie.txt")), status='unknown', action='write')

 do j=1, size(sumax_serie)
        write(160,*)  sumax_serie(j)

 end do 



 close(160)



end subroutine




END MODULE mo_readdata
