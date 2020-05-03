      program panel_method_psi
!______________________________________________________________________
!
!     Soluzione di psi/2 = int (psi dg/dn - dpsi/dn g ) dtau
!
!______________________________________________________________________
      implicit none
!
      include 'parameter.inc'
      include 'Gauss_parameters.inc'
      real(8) :: bh(m,m_sing)
      real(8) :: A(m,m),B(m,m),bc(m),tn(m)
      real(8) :: Af(m),Bf(m),Cf(m),Df(m),Ef(m)
      real(8) :: Aw(m),Bw(m),Cw(m),Dw(m),Ew(m)
      real(8) :: utau(m),cp_wall(m)
!
      integer :: i,j,n
      real(8) :: theta,dtheta
      real(8) :: R
      real(8) :: Uinf
      real(8) :: xa(m+1),xb(m),ya(m),yb(m) 
      real(8) :: xs(m),ys(m),tauc(m)
      real(8) :: hl(m),taux(m),tauy(m),hnx(m),hny(m)
      real(8) :: AA,BB,ttn
      real(8) :: xf,yf,rf,thf
      real(8) :: sumA,sumB
!______________________________________________________________________
      call input(n,R,Uinf)
!
      call geom_circ(n,R,xs,ys,tauc,xa,ya,xb,yb,hl,taux,tauy,hnx,hny)
!
      call outgeo(n,xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny)
!
      call Mat(n,xa,xb,ya,yb,xs,ys,hl,taux,tauy,hnx,hny,A,B)
!!    call check_coeff(n,A,B)
!
      call b_cond(n,Uinf,ys,bc)
!
      call t_not(n,bc,A,tn)
!
!     A . x = tn
      if(i_gauss_ctrl.eq.1) call Gauss_PP(B,tn,bh,n)
      if(i_gauss_ctrl.eq.0) call Gauss_PP0(B,tn,n)
      call outsol(n,tauc,tn,bc,utau)
!_______________________________________________________________________
      end
!______________________________________________________________________
      subroutine input(n,R,Uinf)
!______________________________________________________________________
!
!     Input per panel_method
!______________________________________________________________________
      implicit none
      include 'parameter.inc'
      include 'panel_out.inc'
      include 'Gauss_parameters.inc'
      include 'Grid_parameters.inc'
!
      integer :: n
      real(8) :: R
      real(8) :: Uinf
!
      character(15) :: cmn_string
!______________________________________________________________________
!              
      filen = './panel_method.dat'
      open(unit=1,file=filen)
        read(1,*) dirname
        read(1,*) ver
        read(1,*) n
!
        if(n.gt.m) stop 'Err1: n > m  in Input'
!
        read(1,*) R
        read(1,*) Uinf
!
        read(1,*) nxf
        read(1,*) nyf
        read(1,*) xfI
        read(1,*) yfI
        read(1,*) xfF
        read(1,*) yfF
!
        read(1,*) i_gauss_ctrl

      close(1)
!
      cmn_string = '     dirname = '
      write(*,1001) cmn_string,dirname
      cmn_string = '         ver = '
      write(*,1001) cmn_string,ver
      cmn_string = '           n = '
      write(*,1002) cmn_string,n
      cmn_string = '           R = '
      write(*,1003) cmn_string,R
      cmn_string = '        Uinf = '
      write(*,1003) cmn_string,Uinf
      cmn_string = '         nxf = '
      write(*,1002) cmn_string,nxf
      cmn_string = '         nyf = '
      write(*,1002) cmn_string,nyf
      cmn_string = '         xfI = '
      write(*,1003) cmn_string,xfI
      cmn_string = '         yfI = '
      write(*,1003) cmn_string,yfI
      cmn_string = '         xfF = '
      write(*,1003) cmn_string,xfF
      cmn_string = '         yfF = '
      write(*,1003) cmn_string,yfF
      cmn_string = 'i_gauss_ctrl = '
      write(*,1004) cmn_string,i_gauss_ctrl
!
1000  format(a)
1001  format(a,a)
1002  format(a,i4.4)
1003  format(a,g11.4)
1004  format(a,i1.1)
!______________________________________________________________________
      return
      end
!_______________________________________________________________________
!
      subroutine geom_circ(n,R,xs,ys,tauc,xa,ya,xb,yb,hl,taux,tauy,hnx,hny)
!______________________________________________________________________
!
!     Genera un geometria circolare
!______________________________________________________________________
      implicit none
!
      include 'parameter.inc'
      real(8) :: xa(m+1),xb(m),ya(m),yb(m) 
      real(8) :: xs(m),ys(m),tauc(m)
      real(8) :: hl(m),taux(m),tauy(m),hnx(m),hny(m)
!
      integer :: n,i
      real(8) :: R
      real(8) :: dtheta,theta
      real(8) :: pi
!______________________________________________________________________
      pi = acos(-1.d0)
      dtheta = 2.d0*pi/float(n)
      theta = 0.d0
      xa(1) = R*cos(theta)
      ya(1) = R*sin(theta)
      do i=1,n
         theta = float(i)*dtheta
         xb(i) = R*cos(theta)
         yb(i) = R*sin(theta)
         call panel(xa(i),ya(i),xb(i),yb(i),hl(i),taux(i),tauy(i),hnx(i),hny(i))
         xs(i) = (xa(i) + xb(i))*0.5d0
         ys(i) = (ya(i) + yb(i))*0.5d0
!
!        write(*,*) 'gem_circ: i = ',i,' xa = ',xa(i),' ya = ',ya(i)
!
         xa(i+1) = xb(i)
         ya(i+1) = yb(i)
      end do
!
      tauc(1) = hl(1)*0.5
      do i=2,n
         tauc(i) = tauc(i-1) + (hl(i)+hl(i-1))*0.5
      end do
!______________________________________________________________________
      return
      end
!_______________________________________________________________________
      subroutine panel(xa,ya,xb,yb,hl,taux,tauy,hnx,hny)
!______________________________________________________________________
!
!     defines panel geometry:
!     lenght hl, tangent, taux,tauy and normal hnx,hny
!
!______________________________________________________________________
      implicit none
      real(8) :: xa,xb,ya,yb
      real(8) :: hl,taux,tauy,hnx,hny
!______________________________________________________________________
      hl = sqrt((xb-xa)**2 + (yb-ya)**2)
      taux = (xb-xa)/hl
      tauy = (yb-ya)/hl
      hnx  = -tauy
      hny  = taux
!______________________________________________________________________
      return
      end 
!______________________________________________________________________
!______________________________________________________________________
      subroutine Mat(n,xa,xb,ya,yb,xs,ys,hl,taux,tauy,hnx,hny,A,B)
!______________________________________________________________________
!
!     Costruisce le matrici A e B per panel_method_psi
!
!     BB * x = (1/2*I - AA) *psi
!     dove per x si intende la d(psi)/dn che sarà la 
!     nostra incognita             
!               
!     AA = int dg/dn n:= normale esternra al dominio fluido/interna al corpo
!     BB = int g dl
!
!     B * x = A * psi
!     B = BB
!     A = 1/2*I -AA
!______________________________________________________________________
      include 'parameter.inc'
      real(8) :: A(m,m),B(m,m)
!
      integer :: n
      real(8) :: xa(m+1),xb(m),ya(m),yb(m) 
      real(8) :: xs(m),ys(m)
      real(8) :: hl(m),taux(m),tauy(m),hnx(m),hny(m)
!
      integer :: i,j
      real(8) :: AA,BB,ttn,sumA
!______________________________________________________________________
!!    calcolo coefficienti
      do i=1,n
         do j=1,n
            call coeffA(xs(i),ys(i),xa(j),ya(j),xb(j),yb(j),& 
                        hl(j),taux(j),tauy(j),hnx(j),hny(j),AA)
            A(i,j) = AA
            call coeffB(xs(i),ys(i),xa(j),ya(j),xb(j),yb(j),&
                        hl(j),taux(j),tauy(j),hnx(j),hny(j),BB)
            B(i,j) = BB
         end do
         sumA = 0.d0
         do j=1,n
            sumA = sumA + A(i,j)
         end do
         sumA = sumA - A(i,i)
         A(i,i) = - sumA
      end do
!
      do i=1,n
         do j=1,n
            A(i,j) = -A(i,j)
         end do
         A(i,i) = 0.5d0 + A(i,i)
      end do
!_____________________________________________________________________
      return
      end
!_____________________________________________________________________
!______________________________________________________________________
      subroutine coeffA(xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny,A)
!______________________________________________________________________
!
!     A = int dg/dn dl
!
!______________________________________________________________________
      implicit none
      real(8) :: xa,xb,ya,yb
      real(8) :: xs,ys
      real(8) :: hl,taux,tauy,hnx,hny
      real(8) :: A
      real(8) :: sn,cs,atg,atan_C
      real(8) :: xm,ym
      real(8) :: taus,hns
      real(8) :: tau,d
      real(8) :: atgA,atgB
      real(8) :: pi
!______________________________________________________________________
      pi = acos(-1.d0)
      xm = (xa+xb)*0.5d0
      ym = (ya+yb)*0.5d0
      hns  = (xm-xs)*hnx + (ym-ys)*hny
      taus = (xm-xs)*taux+ (ym-ys)*tauy
!
!     B - integration limit
      tau = hl*0.5d0 +taus
      d = sqrt(tau**2 + hns**2)
      cs = -hns/d
      sn = -tau/d
      atgB = atan_C(-hns,sn,cs)
!
!
!     A - integration limit
      tau = -hl*0.5d0 +taus
      d = sqrt(tau**2 + hns**2)
      cs = -hns/d
      sn = -tau/d
      atgA = atan_C(-hns,sn,cs)
!
!
      A = (atgB-atgA)/(2.d0*pi)
!
1000  format(' cs = ',g11.4,'   sn = ',g11.4,' angle = ',g11.4)
1001  format('hns = ',g11.4,' taus = ',g11.4)
!______________________________________________________________________
      return
      end 
!______________________________________________________________________
!
!______________________________________________________________________
      subroutine coeffB(xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny,B)
!______________________________________________________________________
!
!     B = int g dl
!
!______________________________________________________________________
      implicit none
      real(8) :: xa,xb,ya,yb
      real(8) :: xs,ys
      real(8) :: hl,taux,tauy,hnx,hny
      real(8) :: B
!
      real(8) :: xm,ym,hns,taus
      real(8) :: tau,d,sn,cs
      real(8) :: Coe_a1,Coe_a2,Coe_a3
      real(8) :: Coe_b1,Coe_b2,Coe_b3
      real(8) :: Coe_a,Coe_b
      real(8) :: atan_C
      real(8) :: atgA,atgB
      real(8) :: pi
!______________________________________________________________________
      pi = acos(-1.d0)
      xm = (xa+xb)*0.5d0
      ym = (ya+yb)*0.5d0
      hns  = (xm-xs)*hnx + (ym-ys)*hny
      taus = (xm-xs)*taux+ (ym-ys)*tauy
!
!     B - integration limit
      tau = hl*0.5d0 -taus
      d = sqrt(tau**2 + hns**2)
      cs = -hns/d
      sn = -tau/d
      atgB = atan_C(-hns,sn,cs)
!
      Coe_b1 =  tau*log(d)
      Coe_b2 =  -tau
      Coe_b3 =  hns*atgB
!
      Coe_b = Coe_b1+Coe_b2+Coe_b3
!
!
!
!     A - integration limit
      tau = -hl*0.5d0 -taus
      d = sqrt(tau**2 + hns**2)
      cs = -hns/d
      sn = -tau/d
      atgA = atan_C(-hns,sn,cs)
!
      Coe_a1 =  tau*log(d)
      Coe_a2 =  -tau
      Coe_a3 =  hns*atgA
!
      Coe_a = Coe_a1+Coe_a2+Coe_a3
!
!   
      B = (Coe_b-Coe_a)/(2.d0*pi)
!
!
!______________________________________________________________________
      return
      end
!______________________________________________________________________
!_____________________________________________________________________n
      subroutine b_cond(n,Uinf,ys,bc)
!_____________________________________________________________________n
!
!     condizioni al contorno per panel_method_psi
!_____________________________________________________________________n
      implicit none
      include 'parameter.inc'
      integer :: n
!
      real(8) :: ys(m)
      real(8) :: bc(m)
      real(8) :: Uinf
!
      integer :: i
!_____________________________________________________________________n
      do i=1,n
         bc(i) = - Uinf*ys(i) 
      end do
!_____________________________________________________________________n
      return
      end
!_____________________________________________________________________n
!_____________________________________________________________________n
      subroutine t_not(n,bc,A,tn)
!_____________________________________________________________________n
!
!     calcolo termine noto per panel_method_psi
!_____________________________________________________________________n
      implicit none
!
      include 'parameter.inc'
      integer :: n
      real(8) :: A(m,m),bc(m),tn(m)
!
      integer :: i,j
      real :: ttn
!_____________________________________________________________________n
      do i=1,n
         ttn = 0.d0
         do j=1,n
            ttn = ttn + A(i,j)*bc(j)
         end do
         tn(i) = ttn
      end do
!_____________________________________________________________________n
      return
      end
!_____________________________________________________________________n

!_______________________________________________________________________
      subroutine check_coeff(n,A,B)
!_______________________________________________________________________
!
!     effettua la somma dei coefficienti della matrice
!_______________________________________________________________________
      implicit none
!
      include 'parameter.inc'
      integer :: n
      integer :: R
      real(8) :: A(m,m),B(m,m)
!
      integer :: i,j
      real(8) :: sumA,sumB
!_______________________________________________________________________
      do i=1,n
         write(*,*) 'i = ',i,' A(i,i) = ',A(i,i)
         sumA = 0.d0
         sumB = 0.d0
         do j=1,n
            sumA = sumA + A(i,j)
            sumB = sumB + B(i,j)
         end do
         write(*,*) 'Check_coeff: i = ',i,' sumAn = ',sumA,' sumBn = ',sumB
      end do
!_______________________________________________________________________
      return
      end
!_______________________________________________________________________
!
!_______________________________________________________________________
      subroutine Gauss_PP(A,b,bh,n)
!______________________________________________________________________
!
!     Gaussian elimination with partial pivoting
!______________________________________________________________________
      implicit none
      include 'parameter.inc'
      include 'Gauss_parameters.inc'
      real(8) :: A(m,m),b(m),bh(m,m_sing)
      integer :: n
!
      integer :: k,kk
      integer :: l,s
      integer :: kmax
      integer :: is
      real(8) :: amax2,aa2
      real(8) :: av,an,bv,bn
      character(13) :: string2,string3
      character(16) :: string
!
      real(8), parameter :: epsilon = 1.d-8
!______________________________________________________________________
      write(*,*) 'Gauss_PP'
      i_singular = 0
      do k=1,n
!
!        Partial Pivoting
         kmax = k
         amax2 = a(k,k)*a(k,k)
         do kk=k,n
            aa2 = a(kk,k)*a(kk,k)
            if(aa2.gt.amax2) kmax = kk
         end do
         bv = b(k)
         bn = b(kmax)
         b(k)    = bn
         b(kmax) = bv
         do l=k,n 
            av = a(k,l)
            an = a(kmax,l)
            a(k,l)    = an
            a(kmax,l) = av
         end do
!
!        Gauss Elimination
         if(abs(a(k,k)).gt.epsilon) then
            do l=k+1,n
               b(l) = b(l) - a(l,k)*b(k)/a(k,k)
               do s = k+1,n
                  a(l,s) = a(l,s) - a(l,k)*a(k,s)/a(k,k)
               end do
            end do
         else
            i_singular = i_singular+1
         end if
      end do
!
      do k=1,n
         do is=1,m_sing
            bh(k,is) = 0.d0
         end do
      end do
!
!     Singularity check
      if(i_singular.gt.0) then
         do is = 1,i_singular
            k = n - (is-1)
            string = ' Singular matrix'
            write(*,1000) string
            string  = 'pivot k = '
            string2 = ' a(k,k) = '
            write(*,1001) string,k,string2,a(k,k)
            string  = 'pivot k = '
            string2 = '   b(k) = '
            write(*,1001) string,k,string2,b(k)
            if(abs(b(k)).gt.epsilon) then
              string = 'Incompatible rhs'
              write(*,1000) string
              string  = '  rhs k = '
              string2 = ' b(k) = '
              write(*,1001) string,k,string2,b(k)
              stop
            end if
            b(k)  = 0.d0
            bh(k,is) = 1.d0
         end do
      end if
!
      if(i_singular.gt.m_sing) stop 'GaussPP: i_singular > m_sing'
!
!     Back substitution
      do k=n-i_singular,1,-1
         do s = k+1,n
            b(k)  = b(k)  - a(k,s)*b(s)
            do is=1,i_singular
               bh(k,is) = bh(k,is) - a(k,s)*bh(s,is)
            end do
         end do
         b(k)   = b(k)/a(k,k)
         do is=1,i_singular
            bh(k,is) = bh(k,is)/a(k,k)
         end do
      end do
!
1000  format(a)
1001  format(a,i4.4,a,g11.4)
!______________________________________________________________________
      return
      end
!______________________________________________________________________
!______________________________________________________________________
      subroutine Gauss_PP0(A,b,n)
!______________________________________________________________________
!
!     Gaussian elimination with partial pivoting
!     Matrici non singolari
!______________________________________________________________________
      implicit none
      include 'parameter.inc'
      include 'Gauss_parameters.inc'
      real(8) :: A(m,m),b(m)
      integer :: n
!
      integer :: k,kk
      integer :: l,s
      integer :: kmax
      real(8) :: amax2,aa2
      real(8) :: av,an,bv,bn
      character(13) :: string2,string3
      character(16) :: string
!
      real(8), parameter :: epsilon = 1.d-8
!______________________________________________________________________
      do k=1,n
!        Partial Pivoting
         kmax = k
         amax2 = a(k,k)*a(k,k)
         do kk=k,n
            aa2 = a(kk,k)*a(kk,k)
            if(aa2.gt.amax2) then 
              kmax = kk
              amax2 = aa2
            end if
         end do
         bv = b(k)
         bn = b(kmax)
         b(k)    = bn
         b(kmax) = bv
         do l=k,n 
            av = a(k,l)
            an = a(kmax,l)
            a(k,l)    = an
            a(kmax,l) = av
         end do
!
!        Gauss Elimination
         do l=k+1,n
            b(l) = b(l) - a(l,k)*b(k)/a(k,k)
            do s = k+1,n
               a(l,s) = a(l,s) - a(l,k)*a(k,s)/a(k,k)
            end do
         end do
      end do
!
!     Back substitution
      do k=n,1,-1
         do s = k+1,n
            b(k)  = b(k)  - a(k,s)*b(s)
         end do
         b(k)   = b(k)/a(k,k)
      end do
!
1000  format(a)
1001  format(a,i4.4,a,g11.4)
!______________________________________________________________________
      return
      end
!______________________________________________________________________
!
!______________________________________________________________________
      subroutine outsol(n,tauc,tn,bc,utau)
!______________________________________________________________________
!
!     Output soluzione equazione integrale
!______________________________________________________________________
      implicit none
      include 'parameter.inc'
      include 'panel_out.inc'
!
      integer :: n
      real(8) :: tauc(m),tn(m),bc(m),utau(m)
!______________________________________________________________________
      filen_out = 'solution'
      str_intestazione = ' k        tauc          utau    '
      call out_vector(n,tauc,tn)
!
      filen_out = 'boundary_condition'
      str_intestazione = ' k        tauc          bc     '
      call out_vector(n,tauc,bc)
!______________________________________________________________________
      return
      end
!______________________________________________________________________
!______________________________________________________________________
      function atan_C(ys,sn,cs)
!______________________________________________________________________
!
!     Arco tangente in [0,2*pi] :
!          tan_C = atan2(sn,cs)
!          se (ys <0) e (sn <0) => tan_C = atan2(sn,cs) + 2*pi
!______________________________________________________________________
      implicit none
      real(8) :: atan_C,ys,sn,cs
!
      real(8) :: sign_ys,sign_sn,theta,thetaC
      real(8) :: pi
!______________________________________________________________________
      pi = acos(-1.d0)
!
      sign_ys = sign(1.d0,ys)
      sign_sn = sign(1.d0,sn)
      thetaC = 0.5*(1.d0-sign_ys)*(1.d0-sign_sn)*pi
      theta = atan2(sn,cs)
      atan_C = theta + thetaC
!______________________________________________________________________
      return
      end
!______________________________________________________________________
!______________________________________________________________________
      subroutine coeffC(xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny,C)
!______________________________________________________________________
!
!     Integrale dg/dtau
!______________________________________________________________________
      implicit none
!
      real(8) :: xa,xb,ya,yb
      real(8) :: xs,ys
      real(8) :: hl,taux,tauy,hnx,hny
      real(8) :: C
!
      real(8) :: xm,ym,hns,taus
      real(8) :: tau,Ca,Cb
      real(8) :: pi
!______________________________________________________________________
      pi = acos(-1.d0)
      xm = (xa+xb)*0.5d0
      ym = (ya+yb)*0.5d0
      hns  = (xm-xs)*hnx + (ym-ys)*hny
      taus = (xm-xs)*taux+ (ym-ys)*tauy
!
!     B - integration limit
      tau = hl*0.5d0 -taus
      Cb = 0.5*log(hns**2 + tau**2)

!     A - integration limit
      tau = -hl*0.5d0 -taus
      Ca = 0.5*log(hns**2 + tau**2)
!
      C = (Cb - Ca)/(2.d0*pi)
!______________________________________________________________________
      return
      end
!______________________________________________________________________
!______________________________________________________________________
      subroutine coeffD(xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny,D)
!______________________________________________________________________
!
!     Integrale 
!______________________________________________________________________
      implicit none
!
      real(8) :: xa,xb,ya,yb
      real(8) :: xs,ys
      real(8) :: hl,taux,tauy,hnx,hny
      real(8) :: D
!
      real(8) :: xm,ym,hns,taus
      real(8) :: tau,Da,Db,ag
      real(8) :: pi
!______________________________________________________________________
      pi = acos(-1.d0)
      xm = (xa+xb)*0.5d0
      ym = (ya+yb)*0.5d0
      hns  = (xm-xs)*hnx + (ym-ys)*hny
      taus = (xm-xs)*taux+ (ym-ys)*tauy
!
!     B - integration limit
      tau = hl*0.5d0 -taus
      ag = tau/hns
      Db = ag/(1.d0+ag**2)

!     A - integration limit
      tau = -hl*0.5d0 -taus
      ag = tau/hns
      Da = ag/(1.d0+ag**2)
!
      D = (Db - Da)/(2.d0*pi*hns)
!______________________________________________________________________
      return
      end
!______________________________________________________________________
!______________________________________________________________________
      subroutine coeffE(xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny,E)
!______________________________________________________________________
!
!     Integrale 
!______________________________________________________________________
      implicit none
!
      real(8) :: xa,xb,ya,yb
      real(8) :: xs,ys
      real(8) :: hl,taux,tauy,hnx,hny
      real(8) :: E
!
      real(8) :: xm,ym,hns,taus
      real(8) :: tau,Ea,Eb,ag
      real(8) :: pi
!______________________________________________________________________
      pi = acos(-1.d0)
      xm = (xa+xb)*0.5d0
      ym = (ya+yb)*0.5d0
      hns  = (xm-xs)*hnx + (ym-ys)*hny
      taus = (xm-xs)*taux+ (ym-ys)*tauy
!
!     B - integration limit
      tau = hl*0.5d0 -taus
      ag = tau/hns
      Eb = 1.d0/(1.d0+ag**2)

!     A - integration limit
      tau = -hl*0.5d0 -taus
      ag = tau/hns
      Ea = 1.d0/(1.d0+ag**2)
!
      E = (Eb - Ea)/(2.d0*pi*hns)
!______________________________________________________________________
      return
      end
!______________________________________________________________________
!
      subroutine Field_Mat(n,xf,yf,xa,xb,ya,yb,hl,taux,tauy,hnx,hny,Af,Bf,Cf,Df,Ef)
!_______________________________________________________________________
!
!     calcola i coefficienti per il campo
!_______________________________________________________________________
      include 'parameter.inc'
      real(8) :: Af(m),Bf(m),Cf(m),Df(m),Ef(m)
!
      integer :: n
      real(8) :: xf,yf
      real(8) :: xa(m+1),xb(m),ya(m),yb(m) 
      real(8) :: hl(m),taux(m),tauy(m),hnx(m),hny(m)
!
      integer :: i,j
      real(8) :: AA,BB,CC,DD,EE
!_______________________________________________________________________
      do j=1,n
         call coeffA(xf,yf,xa(j),ya(j),xb(j),yb(j),& 
                     hl(j),taux(j),tauy(j),hnx(j),hny(j),AA)
         Af(j) = AA
         call coeffB(xf,yf,xa(j),ya(j),xb(j),yb(j),&
                     hl(j),taux(j),tauy(j),hnx(j),hny(j),BB)
         Bf(j) = BB
         call coeffC(xf,yf,xa(j),ya(j),xb(j),yb(j),&
                     hl(j),taux(j),tauy(j),hnx(j),hny(j),CC)
         Cf(j) = CC
         call coeffD(xf,yf,xa(j),ya(j),xb(j),yb(j),&
                     hl(j),taux(j),tauy(j),hnx(j),hny(j),DD)
         Df(j) = DD
         call coeffE(xf,yf,xa(j),ya(j),xb(j),yb(j),&
                     hl(j),taux(j),tauy(j),hnx(j),hny(j),EE)
         Ef(j) = EE
      end do
!_______________________________________________________________________
      return
      end
!______________________________________________________________________
!______________________________________________________________________
!_____________________________________________________________________n
      subroutine outgeo(n,xs,ys,xa,ya,xb,yb,hl,taux,tauy,hnx,hny)
!_____________________________________________________________________n
!
!     output geometria
!_____________________________________________________________________n
      implicit none
      include 'parameter.inc'
      include 'panel_out.inc'
!
      integer :: i,j,n
      real(8) :: xa(m+1),xb(m),ya(m),yb(m) 
      real(8) :: xs(m),ys(m)
      real(8) :: hl(m),taux(m),tauy(m),hnx(m),hny(m)
!_____________________________________________________________________n
      filen_out = 'node_a_panel_geom'
      str_intestazione = ' k        xa          ya     '
      call out_vector(n+1,xa,ya)
!
      filen_out = 'node_b_panel_geom'
      str_intestazione = ' k        xa          ya     '
      call out_vector(n,xb,yb)
!
      filen_out = 'center_panel_geom'
      str_intestazione = ' k        xs          ys     '
      call out_vector(n,xs,ys)
!
      filen_out = 'tangents_panel_geom'
      str_intestazione = ' k       taux        tauy    '
      call out_vector(n,taux,tauy)
!
      filen_out = 'normals_panel_geom'
      str_intestazione = ' k        hnx         hny    '
      call out_vector(n,hnx,hny)
!_____________________________________________________________________n
      return
      end
!_____________________________________________________________________n
!_____________________________________________________________________n
      subroutine out_vector(n,x,y)
!_______________________________________________________________________
      implicit none
      include 'parameter.inc'
      include 'panel_out.inc'
! 
      integer :: n
      real(8) :: x(m),y(m)
!
      integer :: k
!_______________________________________________________________________
      call sver()
      call strng()
      open(unit=1,file=outfile)
         write(1,1001)
         write(1,1002)
         write(1,1001)
         write(1,1005) str_intestazione
         do k=1,n
           write(1,1004) k,x(k),y(k)
         end do
      close(unit=1)
!
!            '1234*12345678901*12345678901
1001  format('------------------------------')
1002  format('           output             ')
1004  format(i4.4,2(1x,g11.4))
1005  format(a)
!_______________________________________________________________________
      return
      end
!_______________________________________________________________________
!_______________________________________________________________________
      subroutine strng()
!_______________________________________________________________________
!     generates the output file name nome: outfile = dirname + filen
!_______________________________________________________________________
      implicit none
      include 'panel_out.inc'
!_______________________________________________________________________
!
      ilength = LEN_TRIM(ADJUSTL(dirname))
      write(outfile(1:),1000) dirname
      ilength = LEN_TRIM(ADJUSTL(outfile))
      write(outfile(ilength+1:),1000) filen
1000  format(a)
!
!!    write(*,*) 'strng: outfile = ',outfile
!_______________________________________________________________________
      return
      end
!_______________________________________________________________________
!_______________________________________________________________________
      subroutine sver()
!_______________________________________________________________________
!     generates the file name adding the version 
!_______________________________________________________________________
      implicit none
      include 'panel_out.inc'
!
!_______________________________________________________________________
!
      write(filen(1:1),1002) ver
      ilength = LEN_TRIM(ADJUSTL(filen_out))
!
      write(filen(2:),1002) filen_out
      ilength = LEN_TRIM(ADJUSTL(filen))
      write(filen(ilength+1:),1000)
1000  format('.dat')
1001  format(i6.6)
1002  format(a)
!
!     write(*,*) 'sver: filen = ',filen
!_______________________________________________________________________
      return
      end
!_______________________________________________________________________
!_______________________________________________________________________
      subroutine stim()
!_______________________________________________________________________
!     generates the file name adding the time idex: filen = filen + it
!_______________________________________________________________________
      implicit none
      include 'panel_out.inc'
!
      integer :: it
!_______________________________________________________________________
!
      ilength = LEN_TRIM(ADJUSTL(filen))
      write(filen(ilength+1:),1001) it
      ilength = LEN_TRIM(ADJUSTL(filen))
      write(filen(ilength+1:),1000)
1000  format('.dat')
1001  format(i6.6)
!
      write(*,*) 'stim: filen = ',filen
!_______________________________________________________________________
      return
      end
!_______________________________________________________________________

