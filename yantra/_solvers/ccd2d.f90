!=======================================================================================
!This File is part of Yantra: A lattice Boltzmann method based tool for multiscale/
!multiphyics simulations
!=======================================================================================
!
!Copyright (C) 2016-2019  <Author> Ravi A. Patel <Email> ravee.a.patel@gmail.com
!
!This program is free software: you can redistribute it and/or modify it under the
!terms of the GNU General Public License as published by the Free Software
!Foundation, either version 3 of the License, or any later version.
!This program is distributed in the hope that it will be useful, but WITHOUT ANY
!WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!This file contains relavent subroutines for optimized implementation
!of  2D convection-conduction equation.
!
!=======================================================================================
!
subroutine edfij(feq,T,u)
    Implicit None
    real(8), Intent(InOut):: feq(5)
    real(8),Intent(In):: T, u(2)
    real(8),parameter::w1 = 2.d0 / 6.d0, w2 = 1.d0 / 6.d0
    feq(1) = w1 * T
    feq(2) = w2 * T * (1.d0+3.d0*u(1))
    feq(3) = w2 * T * (1.d0+3.d0*u(2))
    feq(4) = w2 * T * (1.d0-3.d0*u(1))
    feq(5) = w2 * T * (1.d0-3.d0*u(2))   
end subroutine edfij
!
subroutine Tij (T,f)
  real(8), Intent (In):: f(5)
  real(8), Intent(Out):: T
  T = sum(f)
End subroutine
!
Subroutine gradTij(gradT,f,u,tau)
  real(8), Intent(out)::gradT(2)
  real(8), Intent(In)::f(5)
  real(8), Intent(In)::u(2)
  real(8), Intent(In)::tau
  real(8):: T,OneTauEs2
  real(8),parameter:: es2 = 1.d0/3.d0
  call Tij(T,f)
  OneTauEs2 = -1.d0/(tau*es2)
  gradT(1)= OneTauEs2*((f(2)-f(4))-u(1)*T)
  gradT(2)= OneTauEs2*((f(3)-f(5))-u(2)*T)
End subroutine
!  
Subroutine TotalFluxij(flux,f, u, rho, cp, kappa,  tau)
  real(8), Intent(Out):: flux(2)
  real(8), Intent(In):: f(5)
  real(8), Intent(In):: u(2)
  real(8), Intent(In)::rho,cp,kappa
  real(8), Intent(In):: tau
  real(8):: T,dT(2)
  call Tij(T,f)
  call gradTij(dT,f,u,tau)
  flux (1) = rho * cp * T  * u (1) -kappa*dT(1)
  flux (2) = rho * cp * T  * u (2) -kappa*dT(2)
End subroutine
!
subroutine total_flux(flux,f,u,rho,cp,kappa,tau,nodetype,ly,lx)
  Implicit None
  real(8),Intent(In)::f(ly,lx,5)
  real(8),Intent(Out):: flux(2,ly,lx)
  real(8),Intent(In)::u(2,ly,lx)
  real(8), Intent(In)::rho(ly,lx)
  real(8), Intent(In):: cp(ly,lx)
  real(8), Intent(In):: kappa(ly,lx)
  real(8),Intent(In)::tau(ly,lx)
  real(8),Intent(In):: nodetype(ly,lx)
  Integer, Intent (In) :: ly, lx
  Integer :: i, j
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
  !$OMP DO
  Do j = 1, lx
     Do i = 1, ly
        flux (:,i, j) = 0.d0
        If (nodetype(i, j) <= 0) Then
          call TotalFluxij(flux (:,i,j),f(i,j,:),u(:,i,j),rho(i,j),cp(i,j),kappa(i,j),tau(i,j))
        End If
     End Do
  End Do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine
!
subroutine grad_t(gradT,f,u,tau,nodetype,ly,lx)
  Implicit None
  real(8),Intent(In)::f(ly,lx,5)
  real(8),Intent(Out):: gradT(2,ly,lx)
  real(8),Intent(In)::u(2,ly,lx)
  real(8),Intent(In)::tau(ly,lx)
  real(8),Intent(In):: nodetype(ly,lx)
  Integer, Intent (In) :: ly, lx
  Integer :: i, j
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
  !$OMP DO
  Do j = 1, lx
     Do i = 1, ly
        gradT (:,i, j) = 0.d0
        If (nodetype(i, j) <= 0) Then
          call gradTij(gradT(:, i, j),f(i,j,:),u (:,i, j),tau(i, j))
        End If
     End Do
  End Do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine
!
Subroutine compute_macro_var (f, T, nodetype, ly, lx)
    !computes temprature
      Implicit None
      Real (8), Intent (In) :: f (ly, lx, 5)
      Real (8), Intent (Inout) :: T (ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: i, j
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            T (i, j) = 0.d0
            If (nodetype(i, j) <= 0) Then
               call Tij(T(i, j),f (i, j,:))
            End If
         End Do
	   End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine compute_macro_var
!
Subroutine compute_edf (f, T, u, nodetype, ly, lx)
	!computes equilibrium distribution function
      Implicit None
      Real (8), Intent (out) :: f (ly, lx, 5)
      Real (8), Intent (In) :: T (ly, lx)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: i, j
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            If (nodetype(i, j) <= 0) Then
               call edfij (f(i,j,:),T(i,j),u(:,i,j))
            Else 
                f(i,j,:)=0.d0
            End If
         End Do
      End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine compute_edf
!
Subroutine collide_srt (f, T, u, nodetype, tau, ss,kappa, rho, cp ,periodicity, ly, lx)
    !performs BGK collision
      Implicit None
      Real (8), Intent (Inout) :: f (ly, lx, 5)
      Real (8), Intent (In) :: T (ly, lx)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: tau (ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)   
      Real (8), Intent (In) :: ss (ly, lx)
      Real (8), Intent (In) :: kappa(ly,lx)
      Real (8), Intent (In) :: rho(ly,lx)
      Real (8), Intent (In) :: cp(ly,lx)
      Integer, Intent(In) :: periodicity(2)       
      Integer, Intent (In) :: ly, lx
      Integer :: i, j, k
      Real (8) :: omega, omega1, ssij,feq(5), ftemp, a(2,ly,lx), rhoCp(ly,lx)
      Real(8),parameter, dimension(5):: w=(/2.d0 / 6.d0,1.d0 / 6.d0,1.d0 / 6.d0,&
      	& 1.d0 / 6.d0, 1.d0 / 6.d0/)
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,omega,omega1,ssij,feq, ftemp)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            rhoCp(i,j) =  (rho(i,j) * cp(i,j))
            call gradTij(a(:,i,j),f(i,j,:),u(:,i,j),tau(i,j))
            a(:,i,j) = -kappa(i,j)*a(:,i,j)
         End Do
     End Do
    !$OMP End Do
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            If (nodetype(i, j) <= 0) Then
               ssij = ss (i, j)+sscorr(i,j)
               omega = (1./tau(i, j))
               omega1 = 1. - omega
               call edfij(feq,T(i,j),u(:,i,j))
               f (i, j, :) = omega1 * f (i, j, :) + omega * feq + w * ssij
             !f2 becomes f4 and f4 becomes f2
             !f3 becomes f5 and f5 becomes f3
             !swap the data
               Do k = 2, 3
                  ftemp = f (i, j, k+2)
                  f (i, j, k+2) = f (i, j, k)
                  f (i, j, k) = ftemp
               End Do
            End If
         End Do
      End Do
	!$OMP END DO
  !$OMP END PARALLEL
  contains 
  		real(8) function sscorr(i,j) 
           !second order finite difference at interior points 
           !and first order at boundary if not periodic
            Integer, Intent(In)::i,j
            real(8)::mkdTdOnebyrhoCp(2),aij,rhoCp_p,&
            & rhoCp_n,dOneByRhoCp_p,dOneByRhoCp_n
            Integer:: q, nextX,nextY
            mkdTdOnebyrhoCp = 0.d0
            !x-gradient
            nextX = j+ 1
            nextY = i 
            If (periodicity(1)>0) Then
              If (nextX < 1) nextX = lx
              If (nextX > lx) nextX = 1
            Else
              If (nextX < 1) nextX= 1
              If (nextX > lx) nextX = lx
            End If 
            If (nodetype(nextY,nextX) > 0) nextX = j
            rhoCp_p = 0.5*(rhoCp(i,j) + rhoCp(nextY,nextX))
            !TODO: correction needed for nodetype > 0
            nextX = j- 1
            nextY = i 
            If (periodicity(1)>0) Then
              If (nextX < 1) nextX = lx
              If (nextX > lx) nextX = 1
            Else
              If (nextX < 1) nextX= 1
              If (nextX > lx) nextX = lx
            End If 
            If (nodetype(nextY,nextX) > 0) nextX = j
            rhoCp_n = 0.5*(rhoCp(i,j) + rhoCp(nextY,nextX))
            dOneByRhoCp_n = 2.d0*((1.d0/rhoCp_p)-(1.d0/rhoCp(i,j)))
            dOneByRhoCp_p = 2.d0*((1.d0/rhoCp(i,j))-(1.d0/rhoCp_n))
            aij= a(1,i,j)
            mkdTdOnebyrhoCp(1)=max(aij,0.d0)*dOneByRhoCp_n + min(aij,0.d0)*dOneByRhoCp_p
            !y-gradient
            nextX = j 
            nextY = i - 1
            If (periodicity(2)>0) Then
              If (nextY < 1) nextY = ly
              If (nextY > ly) nextY = 1
            Else
              If (nextY < 1) nextY = 1
              If (nextY > ly) nextY = ly
            End If 
            If (nodetype(nextY,nextX) > 0) nextY = i
            rhoCp_p = 0.5*(rhoCp(i,j) + rhoCp(nextY,nextX))
            nextX = j 
            nextY = i + 1
            If (periodicity(2)>0) Then
              If (nextY < 1) nextY = ly
              If (nextY > ly) nextY = 1
            Else
              If (nextY < 1) nextY = 1
              If (nextY > ly) nextY = ly
            End If 
            If (nodetype(nextY,nextX) > 0) nextY = i
            rhoCp_n = 0.5*(rhoCp(i,j) + rhoCp(nextY,nextX))
            dOneByRhoCp_n = 2.d0*((1.d0/rhoCp_p)-(1.d0/rhoCp(i,j)))
            dOneByRhoCp_p = 2.d0*((1.d0/rhoCp(i,j))-(1.d0/rhoCp_n))
            aij= a(2,i,j)
            mkdTdOnebyrhoCp(2)=max(aij,0.d0)*dOneByRhoCp_n + min(aij,0.d0)*dOneByRhoCp_p
            sscorr = sum(mkdTdOnebyrhoCp)
  		end function
End Subroutine collide_srt
!
Subroutine collide_trt (f, T, u, nodetype, tau_a, MagicPara, ss, kappa, rho, cp, periodicity, ly, &
& lx)
    !performs trt collision substep
      Implicit None
      Real (8), Intent (Inout) :: f (ly, lx, 5)
      Real (8), Intent (In) :: T (ly, lx)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: tau_a (ly, lx)
      Real (8), Intent (In) :: MagicPara
      Real (8), Intent (In) :: nodetype (ly, lx)
      Real (8), Intent (In) :: ss (ly, lx)
      Real (8), Intent (In) :: kappa(ly,lx)
      Real (8), Intent (In) :: rho(ly,lx)
      Real (8), Intent (In) :: cp(ly,lx)
      Integer, Intent (In) :: periodicity(2)
      Integer, Intent (In) :: ly, lx
      Integer :: i, j, k
      Real (8) :: feq(5), fseq(5), faeq(5), fs(5), &
     & fa(5), omega_a, omega_s, tau_s, ssij, ftemp,rhoCp(ly,lx),a(2,ly,lx)
      Real(8),parameter, dimension(5):: w=(/2.d0 / 6.d0,1.d0 / 6.d0,1.d0 / 6.d0,&
      	& 1.d0 / 6.d0, 1.d0 / 6.d0/)
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,feq,fseq,faeq,fs,ftemp,fa,tau_s,omega_a,omega_s,ssij)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            rhoCp(i,j) =  (rho(i,j) * cp(i,j))
            call gradTij(a(:,i,j),f(i,j,:),u(:,i,j),tau_a(i,j))
            a(:,i,j) = -kappa(i,j)*a(:,i,j)
         End Do
     End Do
    !$OMP End Do
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            If (nodetype(i, j) <= 0) Then
               ssij = ss (i, j)+sscorr(i,j)
               tau_s = 0.5 + (MagicPara/(tau_a(i, j)-0.5))
               omega_a = 1.d0 / tau_a (i, j)
               omega_s = 1.d0 / tau_s
               call edfij(feq,T(i,j),u(:,i,j))
               fseq(1) = feq(1)
               faeq(1) = 0.d0
               fs(1) = f(i,j,1)
               fa(1) = 0.d0
               Do k = 2,3
                  fseq(k) = (feq(k)+feq(k+2)) * 0.5
                  faeq(k) = (feq(k)-feq(k+2)) * 0.5
                  fs(k) = (f(i, j, k)+f(i, j, k+2)) * 0.5
	              fa(k) = (f(i, j, k)-f(i, j, k+2)) * 0.5
			   End Do
               Do k = 4,5
                  fseq(k) = fseq(k-2)
                  faeq(k) =-faeq(k-2)
                  fs(k) = fs(k-2)
	              fa(k) =-fa(k-2)
	           End Do
               f (i, j, :) = f (i, j, :) + omega_s * (fseq-fs)  + &
              & omega_a * (faeq-fa) + w * ssij
              !f2 becomes f4 and f4 becomes f2
              !f3 becomes f5 and f5 becomes f3
               Do k = 2, 3
                !swap the data
                  ftemp = f (i, j, k+2)
                  f (i, j, k+2) = f (i, j, k)
                  f (i, j, k) = ftemp
               End Do
            End If
         End Do
      End Do
	!$OMP END DO
	!$OMP END PARALLEL
	contains 
  		real(8) function sscorr(i,j) 
           !second order finite difference at interior points 
           !and first order at boundary if not periodic
            Integer, Intent(In)::i,j
            real(8)::mkdTdOnebyrhoCp(2),aij,rhoCp_p,&
            & rhoCp_n,dOneByRhoCp_p,dOneByRhoCp_n
            Integer:: q, nextX,nextY
            mkdTdOnebyrhoCp = 0.d0
            !x-gradient
            nextX = j+ 1
            nextY = i 
            If (periodicity(1)>0) Then
              If (nextX < 1) nextX = lx
              If (nextX > lx) nextX = 1
            Else
              If (nextX < 1) nextX= 1
              If (nextX > lx) nextX = lx
            End If 
            If (nodetype(nextY,nextX) > 0) nextX = j
            rhoCp_p = 0.5*(rhoCp(i,j) + rhoCp(nextY,nextX))
            nextX = j- 1
            nextY = i 
            If (periodicity(1)>0) Then
              If (nextX < 1) nextX = lx
              If (nextX > lx) nextX = 1
            Else
              If (nextX < 1) nextX= 1
              If (nextX > lx) nextX = lx
            End If 
            If (nodetype(nextY,nextX) > 0) nextX = j
            rhoCp_n = 0.5*(rhoCp(i,j) + rhoCp(nextY,nextX))
            dOneByRhoCp_n = 2.d0*((1.d0/rhoCp_p)-(1.d0/rhoCp(i,j)))
            dOneByRhoCp_p = 2.d0*((1.d0/rhoCp(i,j))-(1.d0/rhoCp_n))
            aij= a(1,i,j)
            mkdTdOnebyrhoCp(1)=max(aij,0.d0)*dOneByRhoCp_n + min(aij,0.d0)*dOneByRhoCp_p
            !y-gradient
            nextX = j 
            nextY = i - 1
            If (periodicity(2)>0) Then
              If (nextY < 1) nextY = ly
              If (nextY > ly) nextY = 1
            Else
              If (nextY < 1) nextY = 1
              If (nextY > ly) nextY = ly
            End If 
            If (nodetype(nextY,nextX) > 0) nextY = i
            rhoCp_p = 0.5*(rhoCp(i,j) + rhoCp(nextY,nextX))
            nextX = j 
            nextY = i + 1
            If (periodicity(2)>0) Then
              If (nextY < 1) nextY = ly
              If (nextY > ly) nextY = 1
            Else
              If (nextY < 1) nextY = 1
              If (nextY > ly) nextY = ly
            End If 
            If (nodetype(nextY,nextX) > 0) nextY = i
            rhoCp_n = 0.5*(rhoCp(i,j) + rhoCp(nextY,nextX))
            dOneByRhoCp_n = 2.d0*((1.d0/rhoCp_p)-(1.d0/rhoCp(i,j)))
            dOneByRhoCp_p = 2.d0*((1.d0/rhoCp(i,j))-(1.d0/rhoCp_n))
            aij= a(2,i,j)
            mkdTdOnebyrhoCp(2)=max(aij,0.d0)*dOneByRhoCp_n + min(aij,0.d0)*dOneByRhoCp_p
            sscorr = sum(mkdTdOnebyrhoCp)
  		end function
End Subroutine collide_trt
!
Subroutine stream_and_bounce (f, nodetype, periodicity, ly, lx)
	!Uses two step swap algorithim porposed by J. Latt(2008)
    !performs streaming alonngwith bounce back condition
      Implicit None
      Real (8), Intent (Inout) :: f (ly, lx, 5)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Integer, Intent(In) :: periodicity(2) 
      Integer, Intent (In) :: ly, lx
      Integer :: i, j, k, nextX, nextY, ex (5), ey (5)
      Real (8) :: ftemp
      ex = (/ 0, 1, 0, - 1, 0 /)
      ey = (/ 0, 0, 1, 0, - 1 /)
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,nextX,nextY,ftemp)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            If (nodetype(i, j) <= 0) Then
               Do k = 2, 3
                  nextX = j + ex (k)
                  nextY = i - ey (k)
	            !Apply periodicity
                  If (periodicity(1)>0) Then
                    If (nextX < 1) nextX = lx
                    If (nextX > lx) nextX = 1
                  End If 
                  If (periodicity(2)>0) Then
                    If (nextY > ly) nextY = 1
                    If (nextY < 1) nextY = ly
                  End If                  
                  If (nextX > 0 .And. nextX < lx+1 .And. nextY > 0 &
                 & .And. nextY < ly+1) Then
                     If (nodetype(nextY, nextX) <= 0) Then
                        ftemp = f (nextY, nextX, k)
                        f (nextY, nextX, k) = f (i, j, k+2)
                        f (i, j, k+2) = ftemp
                     End If
                  End If
               End Do
            End If
         End Do
      End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine stream_and_bounce
!
Subroutine boundary_conditions (f, u, nodetype, tau,  kappa, rho, cp, interp, topbc, &
& topval, bottombc, bottomval, leftbc, leftval, rightbc, rightval, &
&  ly, lx)
    !Imposes boundary conditions on straight walls
    !Accepted inputs
    !topbc    = open,c,tot_flux,periodic,nothing
    !bottombc = open,c,tot_flux,periodic,nothing
    !leftbc   = open,c,tot_flux,periodic,nothing
    !rightbc  = open,c,tot_flux,periodic,nothing
      Implicit None
      Real (8), Intent (Inout) :: f (ly, lx, 5)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Character (Len=*), Intent (In) :: topbc
      Character (Len=*), Intent (In) :: bottombc
      Character (Len=*), Intent (In) :: leftbc
      Character (Len=*), Intent (In) :: rightbc
      Real (8), Intent (In) :: topval(lx)
      Real (8), Intent (In) :: bottomval(lx)
      Real (8), Intent (In) :: leftval(ly)
      Real (8), Intent (In) :: rightval(ly)
      Real (8), Intent (In) :: tau (ly, lx)
      Real (8), Intent (In) :: kappa(ly,lx)
      Real (8), Intent (In) :: rho(ly,lx)
      Real (8), Intent (In) :: cp(ly,lx)      
      Integer, Intent (In) :: interp
      Integer, Intent (In) :: ly, lx
      Integer :: i, j
      Real (8) :: bndval, T0, flux0(2), ftemp
    !------------------------
    !open boundary
    !------------------------
    !top boundary
      If (topbc == 'open') Then
    !$OMP PARALLEL Do DEFAULT(SHARED) PRIVATE(j)
         Do j = 1, lx
            f (1, j, 5) = f (2, j, 5)
         End Do
    !$OMP End PARALLEL Do
      End If
    !bottom boundary
      If (bottombc == 'open') Then
    !$OMP PARALLEL Do DEFAULT(SHARED) PRIVATE(j)
         Do j = 1, lx
            f (ly, j, 3) = f (ly-1, j, 3)
         End Do
    !$OMP End PARALLEL Do
      End If
    !left boundary
      If (leftbc == 'open') Then
    !$OMP PARALLEL Do DEFAULT(SHARED) PRIVATE(i)
         Do i = 1, ly
            f (i, 1, 2) = f (i, 2, 2)
         End Do
    !$OMP End PARALLEL Do
      End If
    !right boundary
      If (rightbc == 'open') Then
	!$OMP PARALLEL Do DEFAULT(SHARED) PRIVATE(i)
         Do i = 1, ly
            f (i, lx, 4) = f (i, lx-1, 4)
         End Do
    !$OMP End PARALLEL Do
      End If
    !------------------------
    !flux boundary
    !------------------------
    !top boundary
      If (topbc == 'tot_flux') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,bndval,flux0)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(1, j) <= 0.) Then
              bndval = topval(j)
              If (interp == 1) Then
                call TotalFluxij(flux0,f(2,j,:),u(:,2,j),rho(2,j),cp(2,j),kappa(2,j),tau(2,j))
                bndval = (bndval + 0.5d0*flux0(2))/1.5d0 
             End If
             f (1, j, 5) = OutFFluxBc (bndval, f(1, j, 3), 3, u(2, 1, &
             & j),kappa(1, j), rho(1, j),cp(1, j), tau(1, j))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !bottom boundary
      If (bottombc == 'tot_flux') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,bndval,flux0)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(ly, j) <= 0.) Then
               bndval = bottomval(j)
               If (interp == 1) Then
                  call TotalFluxij(flux0,f(ly-1,j,:),u(:,ly-1,j),rho(ly-1,j),cp(ly-1,j),kappa(ly-1,j),tau(ly-1,j))
                  bndval = (bndval + 0.5d0*flux0(2))/1.5d0
               End If 
               f (ly, j, 3) = OutFFluxBc (bndval, f(ly, j, 5), 5, u(2, ly, &
              & j),kappa(ly, j), rho(ly, j),cp(ly, j), tau(ly, j))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !left boundary
      If (leftbc == 'tot_flux') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,bndval,flux0)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, 1) <= 0.) Then
               bndval = leftval(i)
               If (interp == 1) Then
                call TotalFluxij(flux0,f(i,2,:),u(:,i,2),rho(i,2),cp(i,2),kappa(i,2),tau(i,2))
                bndval = (bndval + 0.5d0*flux0(1))/1.5d0
               End If
               f (i, 1, 2) = OutFFluxBc (bndval, f(i, 1, 4), 4, u(1, i, &
              & 1),kappa(i, 1), rho(i, 1),cp(i, 1), tau(i, 1))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !right boundary
      If (rightbc == 'tot_flux') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,bndval,flux0)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, lx) <= 0.) Then
                bndval = rightval(i)
                If (interp == 1) Then
                 call TotalFluxij(flux0,f(i,lx-1,:),u(:,i,lx-1),rho(i,lx-1),cp(i,lx-1),kappa(i,lx-1),tau(i,lx-1))
                 bndval = (bndval + 0.5d0*flux0(1))/1.5d0
                End If
               f (i, lx, 4) = OutFFluxBc (bndval, f(i, lx, 2), 4, u(1, &
              & i, lx),kappa(i, lx), rho(i, lx),cp(i, lx), tau(i, lx))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !------------------------
    !grad_T boundary
    !------------------------
    !top boundary
      If (topbc == 'grad_T') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,bndval,flux0)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(1, j) <= 0.) Then
              bndval = topval(j)
              If (interp == 1) Then
                call gradTij(flux0,f(2,j,:),u(:,2,j),tau(2,j))
                bndval = (bndval + 0.5d0*flux0(2))/1.5d0 
             End If
             f (1, j, 5) = OutFgradBc (bndval, f(1, j, 3), 3, u(2, 1, &
             & j), tau(1, j))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !bottom boundary
      If (bottombc == 'grad_T') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,bndval,flux0)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(ly, j) <= 0.) Then
               bndval = bottomval(j)
               If (interp == 1) Then
                  call gradTij(flux0,f(ly-1,j,:),u(:,ly-1,j),tau(ly-1,j))
                  bndval = (bndval + 0.5d0*flux0(2))/1.5d0
               End If 
               f (ly, j, 3) = OutFgradBc (bndval, f(ly, j, 5), 5, u(2, ly, &
              & j), tau(ly, j))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !left boundary
      If (leftbc == 'grad_T') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,bndval,flux0)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, 1) <= 0.) Then
               bndval = leftval(i)
               If (interp == 1) Then
                call gradTij(flux0,f(i,2,:),u(:,i,2),tau(i,2))
                bndval = (bndval + 0.5d0*flux0(1))/1.5d0
               End If
               f (i, 1, 2) = OutFgradBc (bndval, f(i, 1, 4), 4, u(1, i, &
              & 1), tau(i, 1))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !right boundary
      If (rightbc == 'grad_T') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,bndval,flux0)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, lx) <= 0.) Then
                bndval = rightval(i)
                If (interp == 1) Then
                 call gradTij(flux0,f(i,lx-1,:),u(:,i,lx-1),tau(i,lx-1))
                 bndval = (bndval + 0.5d0*flux0(1))/1.5d0
                End If
               f (i, lx, 4) = OutFgradBc (bndval, f(i, lx, 2), 4, u(1, &
              & i, lx), tau(i, lx))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !------------------------
    !T boundary
    !------------------------
    !top boundary
      If (topbc == 'T') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,bndval,T0)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(1, j) <= 0.) Then
               bndval = topval(j)
               If (interp == 1) Then
                  call Tij(T0,f(2,j,:))
                  bndval = (bndval + 0.5d0*T0)/1.5d0
               End If 
               f (1, j, 5) = OutFTempBc (bndval, f(1, j, 1), f(1, j, &
              & 2), f(1, j, 3), f(1, j, 4))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !bottom boundary
      If (bottombc == 'T') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,T0)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(ly, j) <= 0.) Then
              bndval= bottomval(j)
              If (interp == 1) Then
                  call Tij(T0,f(ly-1,j,:))
                  bndval = (bndval + 0.5d0*T0)/1.5d0
               End If
               f (ly, j, 3) = OutFTempBc (bndval, f(ly, j, 1), f(ly, j, &
              & 2), f(ly, j, 4), f(ly, j, 5))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !left boundary
      If (leftbc == 'T') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,bndval,T0)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, 1) <= 0.) Then
                bndval = leftval(i)
                If (interp == 1) Then
                  call Tij(T0,f(i,2,:))
                  bndval = (bndval + 0.5d0*T0)/1.5d0
                End If
               f (i, 1, 2) = OutFTempBc (bndval, f(i, 1, 1), f(i, 1, &
              & 3), f(i, 1, 4), f(i, 1, 5))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !right boundary
      If (rightbc == 'T') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,bndval)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, lx) <= 0.) Then
               bndval = rightval(i)
               If (interp == 1) Then
                  call Tij(T0,f(i,lx-1,:))
                  bndval = (bndval + 0.5d0*T0)/1.5d0
               End If
               f (i, lx, 4) = OutFTempBc (bndval, f(i, lx, 1), f(i, lx, &
              & 2), f(i, lx, 3), f(i, lx, 5))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
Contains
	  !function to compute boundary conditions
      Real (8) Function OutFTempBc (bndval, fin1, fin2, fin3, fin4)
         Implicit None
         Real (8), Intent (In) :: bndval, fin1, fin2, fin3, fin4
         Real (8) :: OutF
         OutF = bndval - fin1 - fin2 - fin3 - fin4
         OutFTempBc = OutF
      End Function OutFTempBc
      FUNCTION genbb(fin,a1,a2,a3,u,tau,inIdx) result (fout)
        !implements  a general bounce-back reaction of form
		! a1 C + a2 gradC  = a3
		IMPLICIT NONE
		REAL(8),Intent(IN)::fin,a1,a2,a3,u,tau
		Integer, Intent(IN):: inIdx
		REAL(8)::fout
		Real(8)::neu,den
		Real(8), PARAMETER:: es2 = 1.d0/3.d0
		if (inIdx>3) Then
			neu= (a1/es2) + (a2/(tau*es2))*(1.d0+(u/es2))
			den= (a1/es2) - (a2/(tau*es2))*(1.d0-(u/es2))
		else
			neu= (a1/es2) + (a2/(tau*es2))*(1.d0-(u/es2))
			den= (a1/es2) - (a2/(tau*es2))*(1.d0+(u/es2))
		endif 
		fout = (a3-neu*fin)/ den
	 END FUNCTION genbb
	 !
	 Real(8) Function OutFgradBc (bndval, fin, inIdx, u, tau)
         Implicit None
         Real (8), Intent (In) :: bndval, fin, u, tau
         Integer, Intent (In) :: inIdx
         Real(8):: a1,a2,a3
         a1 = 0.d0
         a2 = 1.d0
         a3 = bndval
         OutFgradBc = genbb(fin,a1,a2,a3,u,tau,inIdx)
     End Function OutFgradBc
     !     !
     Real (8) Function OutFFluxBc (bndval, fin, inIdx, u, kappa, rho, cp, tau)
         Implicit None
         Real (8), Intent (In) :: bndval, fin, u, tau, kappa, rho, cp
         Integer, Intent (In) :: inIdx
         Real(8):: a1,a2,a3
         a1 = rho * cp * u
         a2 = -kappa
         a3 = bndval
         OutFFluxBc = genbb(fin,a1,a2,a3,u,tau,inIdx)
      End Function OutFFluxBc
End Subroutine boundary_conditions
