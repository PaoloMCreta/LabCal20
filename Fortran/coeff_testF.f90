      program coeff_test
!
!     Test coeffcients:
!     A = int dg/dn dl
!     B = int g dl
!     C = int dg/dtau dtau         
!
!______________________________________________________________________
      implicit none
      integer :: i,n
      real(8) :: theta,dtheta
      real(8) :: R
      real(8) :: pi
      real(8) :: xa,xb,ya,yb
      real(8) :: xs,ys
      real(8) :: hl,taux,tauy,hnx,hny
      real(8) :: taus,hns,tauA,tauB
      real(8) :: dA,dB
      real(8) :: A,B,sumA,sumB
      real(8) :: C,sumC    
      character(10) :: cmn_string
      character(43) :: cmn_string2
!______________________________________________________________________
      pi = acos(-1.d0)
!     input
      cmn_string = ' input n'
      write(*,1000) cmn_string
      read(*,*) n
      cmn_string = ' input R'
      write(*,1000) cmn_string
      read(*,*) R
      cmn_string = ' input xs'
      write(*,1000) cmn_string
      read(*,*) xs
      cmn_string = ' input ys'
      write(*,1000) cmn_string
      read(*,*) ys
!
      cmn_string = '    n = '
      write(*,1001) cmn_string,n
      cmn_string = '    R = '
      write(*,1002) cmn_string,R
      cmn_string = '   xs = '
      write(*,1002) cmn_string,xs
      cmn_string = '   ys = '
      write(*,1002) cmn_string,ys
!
1000  format(a)
1001  format(a,i4.4)
1002  format(a,g11.4)
!
      dtheta = 2.d0*pi/float(n)
      theta = 0.d0
      xa = R*cos(theta)
      ya = R*sin(theta)
      open(unit=11,file='./data/nodes_out.dat')
      open(unit=12,file='./data/tangents.dat')
      open(unit=13,file='./data/normales.dat')
      cmn_string2 = '  nodes n = '
      write(11,2000) cmn_string2,n
      cmn_string2 = '  i      theta          x            y      '
      write(11,2001) cmn_string2
      cmn_string2 = 'tangents n= '
      write(12,2000) cmn_string2,n
      cmn_string2 = '  i      theta         taux        tauy     '
      write(12,2001) cmn_string2
      cmn_string2 = ' normals n= '
      write(13,2000) cmn_string2,n
      cmn_string2 = '  i      theta          hnx         hny     '
      write(13,2001) cmn_string2
!
      write(11,2002) i,theta,xa,ya
      sumA = 0.d0
      sumB = 0.d0
      sumC = 0.d0
      do i=1,n
         theta = float(i)*dtheta
         xb = R*cos(theta)
         yb = R*sin(theta)
!
         call panel(xa,ya,xb,yb,xs,ys,taux,tauy,hnx,hny,taus,hns, &
         tauA,tauB,dA,dB)

         write(11,2002) i,theta,xb,yb
         write(12,2002) i,theta-dtheta*0.5,taux,tauy
         write(13,2002) i,theta-dtheta*0.5,hnx,hny
!
         call coeffA(tauA,tauB,taus,hns,dA,dB,A)
         call coeffB(tauA,tauB,taus,hns,dA,dB,B)
         call coeffC(taus,hns,tauA,tauB,C)
         sumA = sumA + A
         sumB = sumB + B
         sumC = sumC + C
!
         xa = xb
         ya = yb
      end do
!
      cmn_string = ' sumA = '
      write(*,1002) cmn_string,sumA
      cmn_string = ' sumB = '
      write(*,1002) cmn_string,sumB
      cmn_string = ' sumC = '
      write(*,1002) cmn_string,sumC
!
      close(11)
      close(12)
      close(13)
!
2000  format(a,i4.4)
2001  format(a)
2002  format(i4.4,3(1x,g11.4))
!______________________________________________________________________
      end
!______________________________________________________________________
!______________________________________________________________________
      subroutine panel(xa,ya,xb,yb,xs,ys,taux,tauy,hnx,hny,taus,hns, &
      tauA,tauB,dA,dB)          
!______________________________________________________________________
!
!     defines panel geometry:
!     lenght hl, tangent, taux,tauy and normal hnx,hny,
!     and also taus,hns,tauA,tauB,dA,dB
!______________________________________________________________________
      implicit none
      real(8) :: xa,xb,ya,yb,xs,ys
      real(8) :: xm,ym
      real(8) :: hl,taux,tauy,hnx,hny
      real(8) :: taus,hns,tauA,tauB
      real(8) :: dA,dB
!______________________________________________________________________
      hl   = sqrt((xb-xa)**2 + (yb-ya)**2)
      taux = (xb-xa)/hl
      tauy = (yb-ya)/hl
      hnx  = (yb-ya)/hl
      hny  = (xb-xa)/hl
      hny  = -1.d0*hny
      xm   = (xa+xb)*0.5d0
      ym   = (ya+yb)*0.5d0
      taus = (xs-xm)*taux+(ys-ym)*tauy
      hns  = (xs-xm)*hnx+(ys-ym)*hny
      tauA = (xa-xm)*taux+(ya-ym)*tauy
      tauB = (xb-xm)*taux+(yb-ym)*tauy
      dA   = sqrt((tauA-taus)**2+hns**2)
      dB   = sqrt((tauB-taus)**2+hns**2)      
!      write(*,*) 'hns =',hns
!      write(*,*) 'taus=',taus
!______________________________________________________________________
      return
      end 
!______________________________________________________________________
!______________________________________________________________________
      subroutine coeffA(tauA,tauB,taus,hns,dA,dB,A)
!______________________________________________________________________
!
!     A = int dg/dn dl
!
!______________________________________________________________________
      implicit none
      real(8) :: tauA,tauB,taus,hns
      real(8) :: dA,dB
      real(8) :: A
      real(8) :: sn,cs,atg,atan_S
      real(8) :: atgA,atgB
      real(8) :: pi
!______________________________________________________________________
      pi   = acos(-1.d0)
!
!     B - integration limit    
      sn   = hns/dB
      cs   = (tauB-taus)/dB
      atgB = atan_S(sn,cs)
!
!     A - integration limit
      sn   = hns/dA
      cs   = (tauA-taus)/dA
      atgA = atan_S(sn,cs)
!
!      write(*,*) 'atgB=',atgB
!      write(*,*) 'atgA=',atgA
      A = -1.d0*(atgB-atgA)/(2.d0*pi)
!______________________________________________________________________
      return
      end 
!______________________________________________________________________
!______________________________________________________________________
      function atan_S(sn,cs)
!______________________________________________________________________
!
!     Given sine, sn, and cosine, cs produces atan(sn/cs) with
!     values in [0, 2.*pi]
!
!______________________________________________________________________
      implicit none
      real(8) :: atan_S,sn,cs
      real(8) :: atg
      real(8) :: pi
!______________________________________________________________________
      pi = acos(-1.d0)
      atg = atan2(sn,cs)
      if(atg.lt.0.d0) atg = 2.d0*pi+atg
      atan_S = atg
!______________________________________________________________________
      return
      end
!______________________________________________________________________
!______________________________________________________________________
       subroutine coeffB(tauA,tauB,taus,hns,dA,dB,B)
!______________________________________________________________________
!
!     B = int g dl
!     Scompongo la primitiva ottenuta in 2 addendi 
!     1) parte della primitiva pari a (1/4pi)(tau-taus)*log(hns^2)
!     2) parte della primitiva pari a "chiamando x=(tau-taus)/hns
!          (hns/4pi)(xlog(1+x^2)-2x+xatan(x))
!     Entrambe le parti andranno valutate negli estremi di integrazione
!                     
!______________________________________________________________________
      implicit none
      real(8) :: tauA,tauB,hns,taus
      real(8) :: dA,dB
      real(8) :: pi
      real(8) :: B
      real(8) :: x,sn,cs
      real(8) :: ad,add,addd
      real(8) :: atan_S
!
      pi   = acos(-1.d0)
!____________________________________________________________________
!      
!     Prima parte dell'integrale
!      
      ad   = ((tauB-tauA)*log(abs(hns)))/(2.d0*pi)
!_____________________________________________________________________
!      Seconda Parte dell'integrale
!
!      Integration limit B
!      
      x   = (tauB-taus)/hns
      sn  = hns/dB
      cs  = (tauB-taus)/dB
      add = x*log(1.d0+x**2)-2.d0*(x+atan_S(sn,cs))
!
!      Integration limit A
!      
      x    = (tauA-taus)/hns
      sn   = hns/dA
      cs   = (tauA-taus)/dA
      addd = x*log(1.d0+x**2)-2.d0*(x+atan_S(sn,cs))
!_____________________________________________________________________
!
!     Somma completa
!
      B = ad+(hns/(4.d0*pi))*(add-addd)
!      
!_____________________________________________________________________
      return
      end
!_____________________________________________________________________
!
!     Calcolo del coeff C che consiste in
!         C = int dg/dtau dtau
!        dove in 2D g=log(r)/2pi
!
      subroutine coeffC(taus,hns,tauA,tauB,C)   
!      
      real(8) :: hns,taus
      real(8) :: tauA,tauB
      real(8) :: pi
      real(8) :: dtauA,dtauB
      real(8) :: C
!
      pi    = acos(-1.d0)
      dtauA = tauA-taus
      dtauB = tauB-taus
      C     = (log(hns**2+dtauA**2)-log(hns**2+dtauB**2))/(4.d0*pi)
!              
      return
      end
!_____________________________________________________________________      
