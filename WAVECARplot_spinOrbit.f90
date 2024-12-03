!Adapted by Rob Lavroff to allow plotting of spin-orbit-coupled WAVECAR (first spinor component unless modified)
!Lavroff, R. H.; Munarriz, J.; Dickerson, C. E.; Munoz, F.; Alexandrova, A. N. 
!Chemical Bonding Dictates Drastic Critical Temperature Difference in Two Seemingly Identical Superconductors. 
!2024, Proc. Natl. Acad. Sci. USA, 121, e2316101121

!!$************************* WAVECARplot*********************************
!!$
!!$   by Karol Cwieka and Kamil Czelej
!!$   (e-mail: karol.cwieka@inmat.pw.edu.pl)
!!$ 
!!$ L.A.Bumm 202-03-29 (bumm@ou.edu) and Gil Speyer (speyer@asu.edu)
!!$ 1) modified to write data outside the inner loop in lines of 5 values
!!$ 2) added OMP parallelization of the innermost loop.  
!!$
!!$   input the WAVECAR file in binary format from VASP, and write
!!$   selected real space wavefunction in a3 direction to standard output
!!$   Output file in the CHGCAR file
!!$
!!$   Compile with gfortran or ifort. Flag "-assume byterecl" is required
!!$   for ifort.
!!$
!!$   Program uses the former implementation of WaveTransPlot by 
!!$   R. M. Feenstra and M. Widom
!!$
!!$   coordinates are direct coordinates

implicit real*8 (a-h, o-z)
complex*8, allocatable :: coeff(:)
complex*16, allocatable :: cener(:)
real*8, allocatable :: occ(:)
real, allocatable :: rphi(:,:)
real, allocatable :: x_csum(:) !!bumm
integer, allocatable :: igall(:,:)
dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3),sumkg(3),vtmp(3)
dimension wk(3),xyz(3),wkpg(3),ig(3)
complex*16 csum
integer kpoint,band,nb1,nb2,nb3
character*75 wavefile, chgfile,outputfile
integer i, statchg, ngx, ngy, ngz
integer atomno, natoms, nrphi, rphicol,rphirow,rphii,rphij,xnrphi
character*150 chgrec
real chgc,chg1,chg2,chg3

print *, "print: WAVECAR plot starting"
write(6,*) "print: WAVECAR plot starting"

!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value; program
!!$*   checks for discrepancy of any results between this and VASP values)

data c/0.262465831d0/ 
!!$*   data c/0.26246582250210965422d0/ 
pi=4.*atan(1.)
natoms=0
nrphi=0
chgc=1.0

!!$ parse arguments
call parse(wavefile,chgfile,outputfile,kpoint,band,nb1,nb2,nb3)

!!$*   input
nrecl=24
open(unit=10,file=wavefile,access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) write(6,*) 'open error - iostat =',iost            
read(unit=10,rec=1) xnrecl,xnspin,xnprec
close(unit=10)
nrecl=nint(xnrecl)
nspin=nint(xnspin)
nprec=nint(xnprec) !! real*8 value, specifies the number type in the file: 45200 = complex*8; 5210 = complex*16
if(nprec.eq.45210) then
   write(0,*) '*** error - WAVECAR_double requires complex*16'
   stop
endif
if(nspin.eq.2) then
   write(0,*) '*** error - Not a spinor WAVECAR. ISPIN =',nspin
   stop
endif
open(unit=10,file=wavefile,access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) write(6,*) 'open error - iostat =',iost
read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), &
     (a3(j),j=1,3)
nwk=nint(xnwk)
nband=nint(xnband)
if (kpoint.gt.nwk) then
   write(0,*) '*** error - selected k=',kpoint,' > max k=',nwk
   stop
endif
if (band.gt.nband) then
   write(0,*) '*** error - selected band=',band,' > max band=',nband
   stop
endif
allocate(occ(nband))
allocate(cener(nband))

!!$*   compute reciprocal properties
call vcross(a2xa3,a2,a3)
Vcell=dot_product(a1,a2xa3)
a1mag=dsqrt(dot_product(a1,a1))
a2mag=dsqrt(dot_product(a2,a2))
a3mag=dsqrt(dot_product(a3,a3))
call vcross(b1,a2,a3)
call vcross(b2,a3,a1)
call vcross(b3,a1,a2)
   b1=2.*pi*b1/Vcell
   b2=2.*pi*b2/Vcell
   b3=2.*pi*b3/Vcell

b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)

phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
call vcross(vtmp,b1,b2)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
      
phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
call vcross(vtmp,b1,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
phi123=abs(asin(sinphi123))
nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
      
phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
call vcross(vtmp,b2,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
phi123=abs(asin(sinphi123))
nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
!! 2* to handle two component spinors
npmax=2*min0(npmaxA,npmaxB,npmaxC)

allocate (igall(3,npmax))
allocate (coeff(npmax))

call chg2output(chgfile,outputfile,statchg,chgrec,chgc,i,chg1,chg2,chg3,atomno,natoms,a1,a2,a3,nb1,nb2,nb3,ngx,ngy,ngz)

!!$*   Find desired wavefunction
irec=3+(kpoint-1)*(nband+1) !+(spin-1)*nwk*(nband+1)
read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
     (cener(iband),occ(iband),iband=1,nband)
nplane=nint(xnplane)

!!$*   FFT grid
if(nb1.eq.0) then
      nb1=2*nb1max
else if (nb1.eq.-1) then
      nb1=ngx-1
else if (nb2.gt.0) then
      nb1=2*nb1
else if(nb1.lt.-1) then
     call help
     stop
end if

if(nb2.eq.0) then
      nb2=2*nb2max
else if(nb2.eq.-1) then
      nb2=ngy-1
else if (nb2.gt.0) then
      nb2=2*nb2
else if(nb2.lt.-1) then
     call help
     stop
end if

if(nb3.eq.0) then
      nb3=2*nb3max
else if (nb3.eq.-1) then
      nb3=ngz-1
else if (nb3.gt.0) then
      nb3=2*nb3
else if(nb3.lt.-1) then
     call help
     stop
end if
write(13,"(/,3I5)") nb1+1, nb2+1, nb3+1

print*, 'nb1max=',nb1max                
print*, 'nb2max=',nb2max            
print*, 'nb3max=',nb3max
print*,  'a1mag= ', a1mag, 'a2mag= ', a2mag, 'a3mag= ', a3mag 
      
!!$*   Calculate plane waves
ncnt=0
do ig3=0,2*nb3max
   ig3p=ig3
   if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
   do ig2=0,2*nb2max
      ig2p=ig2
      if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
      do ig1=0,2*nb1max
         ig1p=ig1
         if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
         do j=1,3
            sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
                 (wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
         enddo
         gtot=sqrt(dot_product(sumkg,sumkg))
         etot=gtot**2/c
         if (etot.lt.ecut) then
            ncnt=ncnt+1
            igall(1,ncnt)=ig1p
            igall(2,ncnt)=ig2p
            igall(3,ncnt)=ig3p
         end if
      enddo
   enddo
enddo

if (2*ncnt.ne.nplane) then
   write(0,*) '*** error - computed ncnt=',2*ncnt, &
        ' != input nplane=',nplane
   stop
endif

irec=irec+band
read(unit=10,rec=irec) (coeff(iplane), iplane=1,nplane)

print*, 'nb1= ',nb1+1,'nb2= ',nb2+1,'nb3= ',nb3+1
nbmax=ceiling(((nb1+1)*(nb2+1)*(nb3+1))/5.)
print*,'nbmax=',nbmax
allocate (rphi(1:5,1:nbmax))
allocate (x_csum(1:nb1+1+5)) !!start array index at 1
ixx = 0

do iz=0,nb3
   z=dble(iz)/dble(1+nb3)
   xyz(3)=z
      do iy=0,nb2
         y=dble(iy)/dble(1+nb2)
         xyz(2)=y
            do  ix=0,nb1
                x=dble(ix)/dble(1+nb1)
                xyz(1)=x
                       csum=cmplx(0.,0.)
                       
!$OMP PARALLEL DO PRIVATE(wkpg,ig), REDUCTION(+:csum)
	
                          do iplane=1,ncnt
                             ig=igall(:,iplane)
                             wkpg=wk+ig
                             csum=csum+coeff(ncnt+iplane)* &
                             cdexp(2.*pi*cmplx(0.,1.)*dot_product(wkpg,xyz))
                          enddo
                          
!$OMP END PARALLEL DO
	
                       csum=csum/dsqrt(Vcell)                       
 !                      nrphi=nrphi+1.  !! seems to be a diagnostic counter

!! output z*a3mag for units of Angstroms                   
!!                       write(6,*) nrphi, sngl(z),sngl(real(csum)),sngl(dimag(csum)) !! original output
!!                       write(13,"(5E18.11E2)") sngl(real(csum)) !! original output
	!!					write(6,*) nrphi, ix, iy, iz !! bumm debug

!! ixx index offset accounts for partial lne not yet written to the oitput file
                       x_csum(ix+ixx+1) = sngl(real(csum)) 

                     
            end do
            
!! write all the values from the inner most loop one complete line at a time
            ixx_len = ixx + nb1 + 1
						ixc = 1
            do ix=1,floor(ixx_len/5.)
            	write(13,"(5E18.11E2)") x_csum(ixc:ixc+4)
            	ixc = ixc + 5;
            end do

!! If there is a partial line of data remains, shift that to the beginning of the accumulation array 
!! so it will be written at the start of the next slug of data   
            i_rem_xx = mod(ixx_len,5)         
            if (i_rem_xx > 0) then
               x_csum(1:i_rem_xx) = x_csum(ixc:ixc+i_rem_xx)
            endif
            ixx = i_rem_xx
 
      end do  
end do
!! write the last bit of data as a partial line 
if (ixx > 0) then
	write(13,"(5E18.11E2)") x_csum(1:ixx)
endif

 
close(13)
close(14)
close(15)

deallocate(igall)
deallocate(coeff)
deallocate(rphi)
stop
end program

!!$*   routine for computing vector cross-product
subroutine vcross(a,b,c)
  implicit real*8(a-h,o-z)
  dimension a(3),b(3),c(3)
  
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross      
      
!!$   parse command line arguments
subroutine parse(wavefile,chgfile,outputfile,kpoint,band,nb1,nb2,nb3)
character*75 wavefile, chgfile,outputfile
integer band,kpoint,nb1,nb2,nb3
character*20 option,value
integer iarg,narg,ia
iarg=iargc()
nargs=iarg/2
wavefile="WAVECAR"
chgfile="CHGCAR"
outputfile="KSwavefunction"
!spin = 1
kpoint = 1
band = 1
nb1= -1
nb2= -1
nb3= -1
if(iarg.ne.2*nargs) then
   call help
endif
do ia=1,nargs
   call getarg(2*ia-1,option)
   call getarg(2*ia,value)
   if(option == "-f") then
      read(value,*) wavefile
   else if(option == "-c") then
      read(value,*) chgfile
   else if(option == "-o") then
      read(value,*) outputfile
   !else if(option == "-s") then
      !read(value,*) spin
   else if(option == "-k") then
      read(value,*) kpoint
   else if(option == "-b") then
      read(value,*) band
   else if(option == "-dx") then
      read(value,*) nb1
   else if(option == "-dy") then
      read(value,*) nb2
   else if(option == "-dz") then
      read(value,*) nb3
   else if(option =="-h") then
      call help
   else
      call help
   endif
enddo
return
end subroutine parse

subroutine chg2output(chgfile,outputfile,statchg,chgrec,chgc,i,chg1,chg2,chg3,atomno,natoms,a1,a2,a3,nb1,nb2,nb3,ngx,ngy,ngz)
implicit real*8 (a-h, o-z)
dimension a1(3),a2(3),a3(3)
character*75 chgfile,outputfile
integer i, statchg, ngx, ngy, ngz
integer atomno, natoms
character*150 chgrec
real chgc,chg1,chg2,chg3

open(unit=12,file=chgfile)
open(unit=13,file=outputfile)

read(12,"(A)",iostat=statchg) chgrec
write(13,"(A)") chgrec
read(12,"(f19.14)") chgc 
write(13,"(f19.14)") chgc          
do i=1,3
read(12,"(f13.6,f12.6,f12.6)") chg1, chg2, chg3
write(13,"(f13.6,f12.6,f12.6)") chg1, chg2, chg3          
end do
read(12,"(A)") chgrec
write(13,"(A)") chgrec
   do
     read(12,"(I6)",advance="no",iostat=statchg) atomno
     natoms=natoms+atomno
     if (statchg/=0)then
         exit
     else
         write(13,"(I6)",advance="no") atomno
      end if
   end do
read(12,"(A)") chgrec
write(13,"(/,A)") trim(chgrec)
   do i=1,natoms
     read(12,"(3f10.6)") chg1, chg2, chg3
     write(13,"(3f10.6)") chg1, chg2, chg3
   end do
read(12,"(/,3I5)") ngx, ngy, ngz

close(unit=12)
return
end subroutine chg2output

subroutine help
write(6,*) 'Syntax:'
write(6,*) 'WAVECARplot -o <output file> -s <spin> -k <k-point> -b <band> -dx <dx> -dy <dy> -dz <dz>'
write(6,*) 'defaults:   -o KSwavefunction -s 1 -k 1 -b 1 -dx 0 -dy 0 -dz 0'
write(6,*) '-dx,-dy,-dz stands for FFT grid discretization; 0 (default) for code-estimated discretization;'
write(6,*) '-1 for ngx, ngy, ngz from CHGCAR; other vaues greater than 0 for custom discretization'  
write(6,*) 'output: wavefunction psi(x,y,z) with z direct coordinate on a3 axis in FFT grid'
stop
end subroutine help
