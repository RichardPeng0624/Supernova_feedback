         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    program originally writted by Elena Redaelli, A.A. 2014/2015
!!!!    and slighlty modified by FB. This vanilla version calculates the standard shock tube
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DATA
integer :: N, i, x, cc, KK
real*8 :: pi,cmpc,cmkpc,yr,kbol,mu,mp
parameter (N=500)
parameter (KK=5000)
parameter(pi=3.141592)
parameter(cmpc=3.085d18)
parameter(cmkpc=1000.*cmpc)
parameter(yr=3.156d7)
parameter(kbol=1.38d-16)
parameter(mu=0.61)
parameter(mp=1.67d-24)
parameter(kev_T=1.16d7)
END MODULE DATA

PROGRAM ZEUS
USE DATA
!!IMPLICIT NONE
real*8 :: xa(N), xb(N), xmax, xmin, deltax, dxa(N), dxb(N), T(N)
real*8 :: d(N), e(N), v(N), P(N), s(N), Temp(N) !DENSITA', ENERGIAINTERNA, VELOCITA', PRESSIONE, MOMENTO
real*8 :: q(N) !VISCOSITA' ARTIFICIALE
real*8 :: g2a(N), g2b(N), g31a(N), g31b(N), dvl1a(N), dvl1b(N)
real*8 :: F1(N), F2(N), F3(N), M(N),  dstar(N),  e_dstar(N), vstar(N) , e_d(N)
real*8 :: divV(N), pmax(1)
real*8 :: time(KK), rs_sedov(KK)
real*8 :: Ecin, Eter,  EterIN
real*8 :: dtmin, tmax, tt, c2, gam, cv, k, t1, t2, t3, tcontR, tcontT, LumX, cfl, T0, rho0, dTem
real*8 :: SedLaw !Costante per la legge di Sedov
integer :: sdr, Num, stamp, ncicli, tmax_list(6), ma
character(len=72) :: result_list
real*8, EXTERNAL :: Cool


!CREAZIONE DOPPIA GRIGLIA (xa e xb)

xmin=0. !!change initial values
xmax=70.*cmpc

!GRIGLIA "a"
do i=2,N
    xa(i)= xmin+(xmax-xmin)*(i-2.)/(N-1.)
end do
xa(1)=-xa(3)
xa(2)=0.

deltax=xa(3)-xa(2)

!GRIGLIA "b"
do i=2, N-1
    xb(i)=0.5*(xa(i)+xa(i+1))
end do
xb(N)=xb(N-1)+(xb(N-1)-xb(N-2))   !! add the last calculated Delta_xb to xb(N-1)

do i=2, N-1
    dxa(i)=xa(i+1)-xa(i)
    dxb(i)=xb(i)-xb(i-1)
end do

dxa(1)=xa(2)-xa(1)
dxa(N)=dxa(N-1)
dxb(1)=dxb(2)
dxb(N)=xb(N)-xb(N-1)

open(20,file='grid.dat')
do i=1,N
   write(20,1001)xa(i),xb(i),dxa(i),dxb(i)
enddo
close(20)
1001 format(4(1pe12.4))

!DEFINIZIONE FATTORI DI SCALA METRICI

sdr=1    !! this parameter selects the type of coordinates: 0 = Cartesian, 1=Spherical

if (sdr==0) then  !! Cartesian !!

    do i=1, N
    g2a(i)=1.
    g2b(i)=1.
    g31a(i)=1.
    g31b(i)=1.

    end do
    do i=1, N-1
    dvl1a(i)=xa(i+1)-xa(i)   !! Note that is centered in xb(i)
    end do
    dvl1a(N)=dvl1a(N-1)
    do i=2, N
    dvl1b(i)=xb(i)-xb(i-1)  !! Note that it is centered in xa(i)
    end do
    dvl1b(1)=dvl1b(2)



else if (sdr==1) then   !! spherical !!
    do i=1, N
    g2a(i)=xa(i)
    g31a(i)=xa(i)
    g2b(i)=xb(i)
    g31b(i)=xb(i)
    end do

    do i=1, N-1
    dvl1a(i)=(xa(i+1)**3-xa(i)**3)/3.
    end do
    dvl1a(N)=dvl1a(N-1)
    do i=2, N
    dvl1b(i)=(xb(i)**3-xb(i-1)**3)/3.
    end do
    dvl1b(1)=dvl1b(2)

end if

!IMPLEMENTAZIONE CONDIZIONI INIZIALI

 gam=1.6667   !! gam=5/3
 cv=1.99d8    !! warning: this is right for gam = 5/3 !!
 tt=0.
 c2=3.
 cfl=0.5
 !T0=1.d4      !! Need to be changed if considering radiative loss
 T0=1.d4
 rho0=2.d-24
 dTem=T0/N


tmax_list=(/ 1,2,4,6,8,10 /) !! list of tmax
result_list='result_1.datresult_2.datresult_4.datresult_6.datresult_8.datresult_0.dat'
cc=1

open(21,file='shock_radius.dat')

do x=1, size(tmax_list)
tmax=tmax_list(x)*1.d4*yr
print *,'-----------tmax=', tmax_list(x),'------------'
tt=0.
T(i)=T0
do i=1, N
          d(i)=rho0
          e(i)=cv*d(i)*T(i)
          if (T(i) .le. 1.d4) then
          T(i)=1.d4
          end if
          P(i)=(gam-1.)*e(i)
      v(i)=0.
end do

!!SN Energy Injection
e(3)=1.d51/(1.3333*pi*xa(4)**3)
e(2)=e(3)
T(3)=e(2)/(cv*rho0)
T(2)=T(3)
P(3)=e(2)*(gam-1.)
P(2)=P(3)

e(1)=e(2)
T(1)=T(2)
P(1)=P(2)

!!Shock of

open(11, file='initial.dat')
    write(11,*)'#  xa-------velocity-----pression-----temperature-----energy----density  '
        do i=1, N
            write(11,1000)xa(i)/cmpc, v(i), p(i), T(i), e(i), d(i), d(i)/(1.4*mp)
        enddo
 close(11)
1000 format(8(1pe12.4))



        ncicli=0

do while (tt<tmax)      !!!! HERE STARTS THE TIME INTEGRATION !!!!!
        ncicli=ncicli+1
!!        if(ncicli.gt.20000) goto 1111

    !!do i=1, N        !! not needed for the shock tube test !!
        !!P(i)=(gam-1.)*e(i)
    !!end do

!CALCOLO DTMIN

        dtmin=1.d30   !! any very large value !!
        cfl=0.01
        p=(gam-1.)*e !!??????????????
    do i=2, N-1
         dtmin=min(dtmin,(xb(i)-xb(i-1))/(abs(v(i))+sqrt(gam*P(i)/d(i))))
    end do
    cfl=cfl*1.1
        if (cfl.ge.0.5) then
        cfl=0.5
        end if
        dtmin=cfl*dtmin
        tt=tt+dtmin
        print*,'ncicli, dtmin = ',ncicli, real(dtmin),real(tt/yr)


!=========SOURCE STEP=======================================
!SUBSTEP I: AGGIORNAMENTO DELLA VELOCITÀ PER GRADIENTE DI P

    do i=2, N-1
        v(i)=v(i)-dtmin*2.*(P(i)-P(i-1))/((d(i)+d(i-1))*dxb(i))
    end do
    CALL BCa(v)


!CALCOLO Q
    do i=2, N-1
        if ((v(i+1)-v(i))<0.) then
            q(i)=C2*d(i)*(v(i+1)-v(i))**2
        else
            q(i)=0.
        end if
    end do
    CALL BCb(q)

!SUBSTEP II: AGGIORNAMENTO PER VISCOSITÀ ARTIFICIALE

    do i=2, N-1
        v(i)=v(i)-dtmin*2.*(q(i)-q(i-1))/((d(i)+d(i-1))*dxb(i))
    end do
    CALL BCa(v)

    do i=2, N-1
        e(i)=e(i)-dtmin*q(i)*(v(i+1)-v(i))/dxa(i)

    end do
    CALL BCb(e)

!SUBSTEP III: AGGIORNAMENTO PER RISCALDAMENTO DA COMPRESSIONE
    do i=2,N-1
        divV(i)=(g2a(i+1)*g31a(i+1)*v(i+1)-g2a(i)*g31a(i)*v(i))/dvl1a(i)
    end do
    CALL BCa(divV)

    do i=2, N-1
        e(i)=e(i)*(1.-0.5*dtmin*(gam-1.)*divV(i))/(1.+0.5*dtmin*(gam-1.)*divV(i))
    end do
    CALL BCb(e)

!!  Here update T when needed (not needed for the shock tube)
!==========Temperature with additional radiative loss term======================================
    do i=2, N-1
        e(i)=e(i)-dtmin*Cool(T(i))*(d(i)/2.17d-24)**2
    end do
    CALL Bcb(e)

    do i=1, N
    T(i)=e(i)/(cv*d(i)) !! new T(i)
    if (T(i) .le. 1.d4) then
    T(i) = 1.d4
    end if
    end do

    do i=2, N-1
    e(i)=cv*d(i)*T(i)
    end do
    CALL Bcb(e)

!!!!!!TRANSPORT STEP (use Upwind first order only)

    do i=2, N-1       !! here define the momentum density
        s(i)=0.5*(d(i)+d(i-1))*v(i)  !! this is at "i" !!
    end do

    CALL BCa(s)

!AGGIORNAMENTO DENSITÀ

    do i=2, N-1       !! here select the value of the density at the interface "i"
        if (v(i)>0.) then
            dstar(i)=d(i-1)     !! at i !!
        else
            dstar(i)=d(i)
        end if
    end do
    dstar(N)=dstar(N-1)
    dstar(1)=dstar(3)

    do i=2, N
        F1(i)=dstar(i)*v(i)*g2a(i)*g31a(i)    !! at i !!
    end do

!AGGIORNAMENTO ENERGIA

    do i=2, N-1
        M(i)=dstar(i)*v(i)
    end do
    CALL BCa(M)


    do i=2, N-1
        if (v(i)>0.) then
            e_dstar(i)=e(i-1)/d(i-1)   !! at i !!
        else
            e_dstar(i)=e(i)/d(i)
        end if
    end do
    e_dstar(N)=e_dstar(N-1)
    e_dstar(1)=e_dstar(3)

    !ORA AGGIORNO LA DENSITÀ
    do i=2, N-1
        d(i)=d(i)-dtmin*(F1(i+1)-F1(i))/dvl1a(i)
    end do
    CALL BCb(d)


    do i=2, N
        F2(i)=e_dstar(i)*M(i)*g2a(i)*g31a(i)
    end do
    CALL BCa(F2)

    do i=2, N-1
        e(i)=e(i)-dtmin*(F2(i+1)-F2(i))/dvl1a(i)
    end do

    CALL BCb(e)


!AGGIORNAMENTO MOMENTO

    do i=2, N-1
        if ((v(i-1)+v(i))*0.5>0) then
            vstar(i)=v(i-1)       !! at i-1/2  !!
        else
            vstar(i)=v(i)
        end if
    end do

    CALL BCb (vstar)

    do i=1, N-1
        F3(i)=vstar(i+1)*0.5*(M(i)+M(i+1))*g2b(i)*g31b(i)   !! questo e' a i+1/2, occhio !!
    end do

    do i=2, N-1
        s(i)=s(i)-dtmin/dvl1b(i)*(F3(i)-F3(i-1))
    end do

    CALL BCa(s)

    do i=2, N-1
        v(i)=2.*s(i)/(d(i)+d(i-1))
    end do

    CALL BCa(v)

enddo       !! here the "do while" ends !!

do i=1,N
    time(i)=(i-1)*dtmin
end do

open(20,file= result_list(cc:cc+11))


do i=1,N  !! write the results in the file "results.dat"
    write (20,1000) log10(xa(i)/cmpc), log10(xb(i)/cmpc), log10(d(i)/(1.4*mp)), v(i)/1.d5, log10(e(i)/d(i)), log10(p(i))
    pmax = maxloc(p)
    if (p(i) == pmax(1)) then
        rs_sedov(i)=(((2*1.d51)/(rho0))**0.2)*(time(i)**0.4)
        write(21,1000) log10(time(i)/yr), log10(xa(i)/cmpc), log10(xb(i)/cmpc), log10(rs_sedov(i)/cmpc)
    end if
end do
7000 format(3(1pe12.4))
!!1000 format(6(1pe12.4))

if (tmax==10*1.d4*yr) then
open(22,file='sedov_radius.dat')
do i=1, KK
    rs_sedov(i)=(((2*1.d51)/(rho0))**0.2)*(time(i)**0.4)
    write(22, 7000) log10(time(i)/yr), log10(rs_sedov(i)/cmpc)
end do
close(22)
end if

close(20)

cc=cc+12

end do !! here the tmax_list loop ends !!
close(21)
END PROGRAM ZEUS


SUBROUTINE BCa(z1) !corrette BC per velocità e momento (riflessione)
USE DATA
IMPLICIT NONE
real*8, dimension (N) :: z1

!z1(2)=0.
!z1(1)=-z1(3)
!z1(N)=z1(N-1)
z1(1)=z1(2)       !! ouflow !!
z1(N)=z1(N-1)

END SUBROUTINE BCa

SUBROUTINE BCb(z2) ! BC di outflow tradizionali
USE DATA
IMPLICIT NONE
real*8, dimension (N) :: z2
z2(1)=z2(2)
z2(N)=z2(N-1)
END SUBROUTINE BCb

!! Cooling Functions
Real*8 FUNCTION Cool(Temp1)
USE DATA
IMPLICIT NONE
Real*8:: Temp1
        if (Temp1>(1.7235D-3*kev_T) .AND. Temp1<(0.02*kev_T)) then
            Cool= 6.72d-22*((Temp1*kev_T)/0.02)**0.6
        else if (Temp1>(0.02*kev_T)) then
            Cool= 1d-22*(8.6*1d-3*(Temp1*kev_T)**(-1.7)+0.058*(Temp1*kev_T)**0.5+0.063)
        else
            Cool= 1.544d-22*((Temp1*kev_T)/0.0017235)**6
        end if

END FUNCTION Cool
