program numerik

real(16) :: c,hbar,mn,G,Ms,R_0,p,h,denklemdosya,orjinal,birimsiz,x,pcek,xkatsayicek,denkcek,xkatsayicek2,denkcek2,pcek2
real(16):: meski,atm,r,m,dPdr,dmdr,u,drhodr,uust,ualt,n0,me,e,jmev,nm2mevm3,psinir,dEdr,p2u,u2p,Energy
real(16)::Ydenge,simetri,fonk,deneme,deneme2
Character :: xlabel,ylabel,filen,alt,ust,isim,map
integer ::n,l,a,n2,l2,a2,sabit,sabit2
logical :: dosyakontrol
c=299792458
map="l"
hbar=1.0545718q-34
pi=4.0*atan(1.0)
mn=1.674927471q-27
me=9.10938q-31

n0=0.16q+45
Ms=1.989q+30

R_0=G*Ms/c**2
u=1.0q-02
sabit=1
open(50,file='basincfonk')
open(52,file='enerjifonk')
do while (u<=4.0)
write(50,*)'crust(',sabit,',1)=',u
write(50,*)'crust(',sabit,',2)=',u2p(u)/(1.60218q+32)
write(52,*)'ener(',sabit,',1)=',u
write(52,*)'ener(',sabit,',2)=',Energy(u)
sabit=sabit+1
u=u+u/50.0
end do

close(50)
close(52)
print*,"bitti"







close(6)
end program numerik


function Energy(u)
real(16)::Energy,u,mn,Ydenge,simetri,hbar,pi,ener(283,2),x
integer :: i
mn=939.56563
pi=4*atan(1.0)
x=u/3.0 - 1.0/3.0
Energy=mn-16.0+0.5*230.0*x**2.0 + (1.0/6.0)*300.0*x**3.0 - (1.0/24.0)*500.0*x**4.0+simetri(u)*(1-Ydenge(u))**2.0

return
end function

function dmdr(r,u)
real(16) :: dmdr,r,c,hbar,mn,G,Ms,gama,K,R_0,p,kacandira,u,n0,rho,Energy
c=299792458
pi=4*atan(1.0)
mn=1.674927471q-27
n0=0.16q+45
rho=(n0)*(u)*Energy(u)*(1.6022q-13)/(c*c)
dmdr=(4*pi)*(r**2.0)*rho

return
end function

function dPdr(r,u,m,p)
real(16) :: dPdr,p,m,c,hbar,mn,G,Ms,gama,K1,R_0,kacandira,rho,r,hbme,me,n0,u,Energy

c=299792458

pi=4.0*atan(1.0)
mn=1.674927471q-27
n0=0.16q+45
G=(6.67408q-11)
Ms=1.989q+30
rho=(n0)*(u)*Energy(u)*(1.6022q-13)/(c*c)
n0=0.16q+45
dPdr=-1*(G/(r**2.0))*(rho+(p)/(c**2.0))*(m+(4*pi*r**3)*(p)/(c**2))/(1-2*G*m/(r*c**2.0))

return
end function

function p2u(p)
real(16)::p2u,p,lebiderya,u,u2p,x
x=0.001
do while (p-u2p(x)>1.0q-10)
x=x+0.003
end do
p2u=x
return
end function



function u2p(u)
real(16)::u2p,u,lebiderya,h,orjinal,birimsiz,a,b,crust(283,2)
integer::leblebi
h=0.01
leblebi=u/h
lebiderya=leblebi*h
u2p=(orjinal(lebiderya+h)*(u-lebiderya)-orjinal(lebiderya)*(u-lebiderya-h))/h

return
end function

function simetri(u)
real(16)::u,simetri,x
x=u/3.0 - 1.0/3.0
simetri=E-
return
end function

function fonk(u)
real(16)::u,fonk,simetri,hbarc,pi,n0
n0=0.16
c=299792458
hbarc=197.32
pi=4.0*atan(1.0)
fonk=4.0*simetri(u)/(hbarc*(3.0*(pi**2.0)*n0*u)**(1.0/3.0))
return
end function


function Ydenge(u)
real(16)::Ydenge,u,ilkbolum,fonk,menemen,paydaninbolumu,ikincibolum
Character::map
menemen=fonk(u)
paydaninbolumu=((81.0*menemen**12.0+6.0*(menemen**9.0))**0.5-9.0*(menemen**6.0))**(1.0/3.0)
ilkbolum=-1.0/(6.0**(1.0/3.0)*paydaninbolumu)
ikincibolum=paydaninbolumu/(6.0**(2.0/3.0)*menemen**3.0) + 1.0
Ydenge=(ilkbolum+ikincibolum)/2.0
return
end function


 


function orjinal(birimsiz)
real(16)::orjinal,birimsiz,Ydenge,u,simetri,sol,sag,h,turevdenge,turevsimetri,kutle,tureveos,Energy
h=1.0q-15
kutle=(1.674927471q-27)*(0.16q+45)
!turevdenge=(Ydenge(birimsiz+h)-Ydenge(birimsiz-h))/(2.0*h)
!turevsimetri=(simetri(birimsiz+h)-simetri(birimsiz-h))/(2.0*h)
tureveos=(Energy(birimsiz+h)-Energy(birimsiz-h))/(2.0*h)

!orjinal=(0.16*((44.2/3.0)*birimsiz**(5.0/3.0)-59.1*birimsiz**2.0+(65.39*2.112/3.112)*birimsiz**3.112))*(1.60218q+32)&
!&+0.16*((2.0**(2.0/3.0)-1.0)*22.1*((2.0/3.0)*birimsiz**(5.0/3.0)-birimsiz**2.0)+30.0*birimsiz**2.0)*(1.60218q+32)
!&+(1.60218q+32)*((-4.0*(1.0-2.0*Ydenge(birimsiz))*turevdenge+turevsimetri*(1.0-2.0*Ydenge(birimsiz))**2.0)*0.16*birimsiz**2.0)
!orjinal=(2.03872092q-20)*(kutle*birimsiz)**3.0
orjinal=(1.60218q+32)*tureveos*0.16*birimsiz**2.0

return
end function
