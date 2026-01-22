program fitTAS
implicit none

integer :: i, j, k, l   ! Variables for cycles
integer :: ncalc        ! Number of theoretical spectra
integer :: points       ! Number of points of each input spectra. It is overwritten
integer :: ndata        ! Number of points of the formatted spectra
integer :: nrand        ! Input to select the number of initial random coefficients
integer :: maxstep      ! Maximum number of minimization steps per random set
integer :: counter      ! Controls the number of steps
integer :: best         ! Position of the least error fit
real*8 :: lowlim, uplim ! Spectral boundaries
real*8 :: test, lowerwf, upperwf        ! Auxiliary variables for formatting
real*8 :: test2, test2upper, test2lower, test2previo, lcont, ucont       ! More auxiliary variables for formatting
real*8 :: conver        ! Convergence criterion on error parameter
real*8 :: diff          ! Error parameter
real*8 :: prediff       ! Error parameter of the previous step
real*8 :: laplace       ! Second derivative of the error parameter
real*8 :: gradnorm      ! Norm of the average gradient
real*8, dimension(2) :: diffvar         ! First derivative of the error parameter for the previous and current steps
real*8, dimension(:), allocatable :: wavenf, intnf      ! Wavelengths and intensities of unformatted spectra
real*8, dimension(:), allocatable :: wavelength, intexp ! Wavelengths and intensities of formatted experimental spectrum
real*8, dimension(:,:), allocatable :: intcalc          ! Intensities of formatted theoretical spectra
real*8, dimension(:), allocatable :: coeff      ! Coefficients array
real*8, dimension(:), allocatable :: precoeff   ! Coefficients array of the previous step
real*8, dimension(:), allocatable :: trial      ! Auxiliary array for the guess spectrum
real*8, dimension(:), allocatable :: pretrial   ! Guess spectrum of the previous step
real*8, dimension(:), allocatable :: optdiff    ! All intermediate errors for each minimization step
real*8, dimension(:), allocatable :: gradient   ! Average gradient
real*8, dimension(:,:), allocatable :: optcoeff ! All intermediate coefficients for each minimization step
real*8, dimension(:,:), allocatable :: opttrial ! All intermediate trial spectra for each minimization step
real*8, dimension(:,:), allocatable :: alltrial        ! Matrix to store the best 20 guess spectra
real*8, dimension(:,:), allocatable :: allcoeff        ! Matrix to store the coefficients of the best 20 spectra
real*8, dimension(:), allocatable :: alldiff   ! Array to store the lowest 20 error parameters
character(len=99) :: expinput   ! Name of experimental spectrum data file
character(len=99) :: frmt1      ! Format variable - Formatted theoretical spectra matrix file header
character(len=99) :: frmt2      ! Format variable - Formatted theoretical spectra matrix file main body
character(len=99) :: frmt3      ! Format variable - Fit summary header
character(len=99) :: frmt4      ! Format variable - Fit summary main body
character(len=99) :: frmt5      ! Format variable - Coefficient optimization header
character(len=99) :: frmt6      ! Format variable - Coefficient optimization main body
character(len=99), dimension(:), allocatable :: calcfile        ! Names of all theoretical spectrum data files
character(len=99) :: specout    ! Output files for best 20 spectra

! Select the spectral window to analyze
write(*,*) "Spectral window to analyze"
write(*,*) "Minimum and maximum wavelengths, integers of 5"
do
 read(*,*) lowlim, uplim
 if (mod(lowlim,5.0d0)/=0.OR.mod(uplim,5.0d0)/=0) then
  write(*,*) "Wavelengths must be integers of 5. Please reintroduce both of them."
 else
  exit
 end if
end do
ndata=1+(int(uplim)-int(lowlim))/5
allocate(wavelength(ndata))
allocate(intexp(ndata))

! Read the experimental spectrum
write(*,*) "Experimental spectrum data file"
read(*,*) expinput
open(11,file=expinput)
points=0
do
 read(11,*,end=42)
 points=points+1
end do
42 continue
rewind(11)
allocate(wavenf(points))
allocate(intnf(points))
do i=1,points
 read(11,*) wavenf(i), intnf(i)
end do
close(11)
if (wavenf(1)>lowlim.OR.wavenf(points)<uplim) then
 write(*,*) "The provided spectrum does not contain enough information around the selected boundaries"
 write(*,*) "Please modify the data or change the spectral window"
 stop
end if

! Give the experimental spectrum the proper format for the analysis
! I could have written this part in a subroutine for further uses, but who cares
do i=1,ndata
 wavelength(i)=lowlim+5*(i-1)
 do j=1,points
  test=wavelength(i)-wavenf(j)
  if (test<0) then
   lowerwf=wavenf(j-1)
   upperwf=wavenf(j)
   exit
  end if
 end do
 do k=0,100
  lcont=k/100.0d0
  ucont=1.0d0-lcont
  test2=wavelength(i)-(lcont*lowerwf+ucont*upperwf)
  if (test2>0) then
   test2upper=test2
   test2lower=dabs(test2previo)
   if (test2upper<test2lower) then
    exit
   else
    lcont=(k-1)/100.0d0
    ucont=1.0d0-lcont
    exit
   end if
  else
   test2previo=test2
  end if
 end do
 intexp(i)=lcont*intnf(j-1)+ucont*intnf(j)
end do
intexp(:)=intexp(:)/maxval(dabs(intexp(:)))
deallocate(wavenf)
deallocate(intnf)

! Read all the theoretical spectra and give the proper format sequentially
write(*,*) "Number of theoretical spectra"
read(*,*) ncalc
write(*,'(A,I0,A)') " ", ncalc, " theoretical spectra will be read in order"
allocate(intcalc(ndata,ncalc))
allocate(calcfile(ncalc))
do l=1,ncalc
 write(*,'(A,I0,A)') " Name of theoretical spectrum data file (", l, ")"
 read(*,*) calcfile(l)
 open(11,file=calcfile(l))
 points=0
 do
  read(11,*,end=73)
  points=points+1
 end do
 73 continue
 rewind(11)
 allocate(wavenf(points))
 allocate(intnf(points))
 do j=1,points
  read(11,*) wavenf(j), intnf(j)
 end do
 close(11)
 if (wavenf(1)>lowlim.OR.wavenf(points)<uplim) then
  write(*,*) "The provided spectrum does not contain enough information around the selected boundaries"
  write(*,*) "Please modify the data or change the spectral window"
  stop
 end if
 do i=1,ndata
  do j=1,points
   test=wavelength(i)-wavenf(j)
   if (test<0) then
    lowerwf=wavenf(j-1)
    upperwf=wavenf(j)
    exit
   end if
  end do
  do k=0,100
   lcont=k/100.0d0
   ucont=1.0d0-lcont
   test2=wavelength(i)-(lcont*lowerwf+ucont*upperwf)
   if (test2>0) then
    test2upper=test2
    test2lower=dabs(test2previo)
    if (test2upper<test2lower) then
     exit
    else
     lcont=(k-1)/100.0d0
     ucont=1.0d0-lcont
     exit
    end if
   else
    test2previo=test2
   end if
  end do
  intcalc(i,l)=lcont*intnf(j-1)+ucont*intnf(j)
 end do
 deallocate(wavenf)
 deallocate(intnf)
end do

! The optimization starts here
write(*,*) "Number of initial random coefficients"
read(*,*) nrand
write(*,*) "Convergence criterion on error parameter (10**-N)"
read(*,*) conver
conver=10**(-1.0d0*conver)
write(*,*) "Maximum number of minimization steps (10**N)"
read(*,*) maxstep
maxstep=10**maxstep
allocate(coeff(ncalc))
allocate(precoeff(ncalc))
allocate(trial(ndata))
allocate(pretrial(ndata))
allocate(alltrial(ndata,nrand))
allocate(allcoeff(nrand,ncalc))
allocate(alldiff(nrand))
allocate(optdiff(2*ncalc))
allocate(gradient(ncalc))
allocate(optcoeff(2*ncalc,ncalc))
allocate(opttrial(ndata,2*ncalc))
write(frmt5,'(A,I0,A)') "(A9,A10,", ncalc, "A20)"
write(frmt6,'(A,I0,A)') "(I9,F10.5,", ncalc, "F20.2)"

do i=1,nrand
 write(specout,'(A,I0.2,A)') "fit", i, ".opt"
 open(11,file=specout)
 do j=1,ncalc
  coeff(j)=rand()
 end do
 call normalize(ncalc,coeff,1.0d0)
 call calcerror(ndata,ncalc,coeff,intcalc,intexp,trial,diff)
 diffvar(2)=diff
 laplace=diff
 counter=0
 write(11,frmt5) "#    Step", "Error", (trim(calcfile(j)),j=1,ncalc)
 write(11,frmt6) counter, diff, coeff(:)
 do while (laplace.gt.conver)
  prediff=diff
  precoeff(:)=coeff(:)
  pretrial(:)=trial(:)
  diffvar(1)=diffvar(2)
  gradnorm=0.0d0
  do j=1,ncalc
   optcoeff(j,:)=coeff(:)
   optcoeff(j,j)=optcoeff(j,j)+1.0d-2
   call calcerror(ndata,ncalc,optcoeff(j,:),intcalc,intexp,opttrial(:,j),optdiff(j))
   optcoeff(j+ncalc,:)=coeff(:)
   optcoeff(j+ncalc,j)=optcoeff(j+ncalc,j)-1.0d-2
   call calcerror(ndata,ncalc,optcoeff(j+ncalc,:),intcalc,intexp,opttrial(:,j+ncalc),optdiff(j+ncalc))
   gradient(j)=(optdiff(j)-optdiff(j+ncalc))/(2.0d-2)
   gradnorm=gradnorm+gradient(j)**2
  end do
  gradnorm=dsqrt(gradnorm)
  if (gradnorm.gt.1.0d-2) then
   gradient(:)=gradient(:)*(1.0d-2/gradnorm)
   coeff(:)=coeff(:)-gradient(:)
  else
   coeff(:)=coeff(:)-gradient(:)
  end if
  if (minval(coeff).lt.0.0d0) then
   coeff(:)=coeff(:)-minval(coeff)
  end if
  call normalize(ncalc,coeff,1.0d0)
  call calcerror(ndata,ncalc,coeff,intcalc,intexp,trial,diff)
  diffvar(2)=dabs(prediff-diff)
  laplace=dabs(diffvar(2)-diffvar(1))
  counter=counter+1
  write(11,frmt6) counter, diff, coeff(:)
  if (counter.gt.maxstep) then
   write(11,'(A,I0)') " Convergence criterion NOT met for initial random set ", i
   exit
  end if
 end do
 if (diff.lt.prediff) then
  allcoeff(i,:)=coeff(:)
  alltrial(:,i)=trial(:)
  alldiff(i)=diff
 else
  allcoeff(i,:)=precoeff(:)
  alltrial(:,i)=pretrial(:)
  alldiff(i)=prediff
 end if
 close(11)
end do

! Print all the information
open(12,file="FormattedExpSpectrum.txt")
do i=1,ndata
 write(12,'(2F15.8)') wavelength(i), intexp(i)
end do
close(12)
open(12,file="FormattedCalcSpectra.txt")
write(frmt1,'(A,I0,A)') "(", ncalc+1, "A15)"
write(frmt2,'(A,I0,A)') "(", ncalc+1, "F15.8)"
write(12,frmt1) "#", (trim(calcfile(i)),i=1,ncalc)
do i=1,ndata
 write(12,frmt2) wavelength(i), intcalc(i,:)
end do
close(12)
write(frmt3,'(A,I0,A)') "(A15, A15, ", ncalc, "A20)"
write(frmt4,'(A,I0,A)') "(A15, F15.5, ", ncalc, "F20.2)"
open(12,file="FitSummary.txt")
write(12,frmt3) "Fit file", "Fit error", (trim(calcfile(i)),i=1,ncalc)
do i=1,nrand
 write(specout,'(A,I0.2,A)') "fit", i, ".txt"
 write(12,frmt4) trim(specout), alldiff(i), allcoeff(i,:)
 open(13,file=specout)
 do j=1,ndata
  write(13,'(2F15.8)') wavelength(j), alltrial(j,i)
 end do
 close(13)
end do
best=minloc(alldiff,1)
write(specout,'(A,I0.2)') "Best: fit", best
write(12,*)
write(12,frmt4) trim(specout), alldiff(best), allcoeff(best,:)
write(12,*)

! Memory release
deallocate(wavelength)
deallocate(intexp)
deallocate(intcalc)
deallocate(calcfile)
deallocate(coeff)
deallocate(precoeff)
deallocate(trial)
deallocate(pretrial)
deallocate(alltrial)
deallocate(allcoeff)
deallocate(alldiff)
deallocate(optdiff)
deallocate(optcoeff)
deallocate(opttrial)
deallocate(gradient)

end program fitTAS

!------------------------------------
subroutine normalize(points,vector,norm)
integer, intent(in) :: points
real*8, dimension(points), intent(inout) :: vector
real*8, intent(in) :: norm
integer :: x
real*8 :: total

total=0.0d0
do x=1,points
 total=total+vector(x)
end do
vector(:)=vector(:)*(norm/total)
return

end subroutine normalize
!------------------------------------
subroutine calcerror(var1,var2,multi,theor,orig,guess,error)
integer, intent(in) :: var1, var2
real*8, dimension(var2), intent(in) :: multi
real*8, dimension(var1,var2), intent(in) :: theor
real*8, dimension(var1), intent(in) :: orig
real*8, dimension(var1), intent(out) :: guess
real*8, intent(out) :: error
integer :: x

guess=matmul(theor,multi)
guess(:)=guess(:)/maxval(dabs(guess))
error=0.0d0
do x=1,var1
 error=error+dabs(orig(x)-guess(x))
end do
return

end subroutine calcerror
!------------------------------------
