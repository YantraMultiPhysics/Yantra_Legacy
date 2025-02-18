!=======================================================================================
!This File is part of Yantra: A lattice Boltzmann method based tool for multiscale/
!multiphyics simulations
!=======================================================================================
!
!Copyright (C) 2016-2017  <Author> Ravi A. Patel <Email> ravee.a.patel@gmail.com
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
!of pore-scale 2D Advection-diffusion equation with varying relaxation parameter
!
!=======================================================================================
!
Subroutine edfij(f,c,u)
	Real(8), Intent(Out):: f(5)
	Real(8), Intent(In):: c,u(2) 
	f (1) = w1 * c
	f (2) = w2 * c * (1.d0+3.d0*u(1))
	f (3) = w2 * c * (1.d0+3.d0*u(2))
	f (4) = w2 * c * (1.d0-3.d0*u(1))
	f (5) = w2 * c * (1.d0-3.d0*u(2))
End Subroutine
!
Subroutine cij(c,f)
	Real(8), Intent(Out):: c
	Real(8), Intent(In):: f(5)
	c = sum(f)
End Subroutine
!
Subroutine fluxij(flux,f,u,tau)
	Real(8), Intent(out)::flux(2)
	Real(8), Intent(in)::f(5),u(2),tau
	Real(8):: c
	call cij(c,f)
	flux (1) = (1./(2.*tau)) * c  * u (1) + (1-(1./(2.*tau))) * (f(2)-f(4))
    flux (2) = (1./(2.*tau)) * c  * u (2) + (1-(1./(2.*tau))) * (f(3)-f(5))
end subroutine
!
Subroutine compute_macro_var (f, c, flux, u, nodetype, tau, ly, lx)
    !computes the concentration and flux
      Implicit None
      Real (8), Intent (In) :: f (ly, lx, 5)
      Real (8), Intent (Inout) :: c (ly, lx)
      Real (8), Intent (Inout) :: flux (2, ly, lx)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: tau (ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: i, j
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            c (i, j) = 0.d0
            flux (:, i, j) = 0.d0
            If (nodetype(i, j) <= 0) Then
               c (i, j) = f (i, j, 1) + f (i, j, 2) + f (i, j, 3) + f &
              & (i, j, 4) + f (i, j, 5)
               flux (1, i, j) = (1./(2.*tau(i, j))) * c (i, j) * u (1, &
              & i, j) + (1-(1./(2.*tau(i, j)))) * (f(i, j, 2)-f(i, j, &
              & 4))
               flux (2, i, j) = (1./(2.*tau(i, j))) * c (i, j) * u (2, &
              & i, j) + (1-(1./(2.*tau(i, j)))) * (f(i, j, 3)-f(i, j, &
              & 5))
            End If
         End Do
		End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine compute_macro_var
!
Subroutine compute_edf (f, c, u, nodetype, ly, lx)
	!computes equilibrium distribution function
      Implicit None
      Real (8), Intent (out) :: f (ly, lx, 5)
      Real (8), Intent (In) :: c (ly, lx)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: i, j
      Real (8) :: ftemp, w1, w2
    !lattice related parameters
      w2 = 1.d0 / 6.d0
      w1 = 2.d0 / 6.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ftemp)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            If (nodetype(i, j) <= 0) Then
              !f1
               f (i, j, 1) = w1 * c (i, j)
             !f2
               ftemp = (1.d0+3.d0*u(1, i, j))
               f (i, j, 2) = w2 * c (i, j) * ftemp
             !f3
               ftemp = (1.d0+3.d0*u(2, i, j))
               f (i, j, 3) = w2 * c (i, j) * ftemp
             !f4
               ftemp = (1.d0-3.d0*u(1, i, j))
               f (i, j, 4) = w2 * c (i, j) * ftemp
             !f5
               ftemp = (1.d0-3.d0*u(2, i, j))
               f (i, j, 5) = w2 * c (i, j) * ftemp
            Else
               f (i, j, 1) = 0.d0
               f (i, j, 2) = 0.d0
               f (i, j, 3) = 0.d0
               f (i, j, 4) = 0.d0
               f (i, j, 5) = 0.d0
            End If
         End Do
      End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine compute_edf
!
Subroutine collide_srt (f, c, u, nodetype, tau, ss, ly, lx)
    !performs BGK collision
      Implicit None
      Real (8), Intent (Inout) :: f (ly, lx, 5)
      Real (8), Intent (In) :: c (ly, lx)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: tau (ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Real (8), Intent (In) :: ss (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: i, j, k, ex (5), ey (5)
      Real (8) :: ftemp, omega, omega1, w1, w2, ssij
    !lattice related parameters
      w2 = 1.d0 / 6.d0
      w1 = 2.d0 / 6.d0
      ex = (/ 0, 1, 0, - 1, 0 /)
      ey = (/ 0, 0, 1, 0, - 1 /)
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,ftemp,omega,omega1,ssij)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            If (nodetype(i, j) <= 0) Then
               ssij = ss (i, j)
               omega = (1./tau(i, j))
               omega1 = 1. - omega
             !compute equillibrium distribution for f1
               ftemp = w1 * c (i, j)
             !collision step for f1
               f (i, j, 1) = omega1 * f (i, j, 1) + omega * ftemp + w1 &
              & * ssij
              !compute equillibrium distribution for f2
               ftemp = (1.d0+3.d0*u(1, i, j))
               ftemp = w2 * c (i, j) * ftemp
             !collision step for f2
               f (i, j, 2) = omega1 * f (i, j, 2) + omega * ftemp + w2 &
              & * ssij
              !compute equillibrium distribution for f3
               ftemp = (1.d0+3.d0*u(2, i, j))
               ftemp = w2 * c (i, j) * ftemp
             !collision step for f3
               f (i, j, 3) = omega1 * f (i, j, 3) + omega * ftemp + w2 &
              & * ssij
              !compute equillibrium distribution for f4
               ftemp = (1.d0-3.d0*u(1, i, j))
               ftemp = w2 * c (i, j) * ftemp
             !collision step for f4
               f (i, j, 4) = omega1 * f (i, j, 4) + omega * ftemp + w2 &
              & * ssij
              !compute equillibrium distribution for f5
               ftemp = (1.d0-3.d0*u(2, i, j))
               ftemp = w2 * c (i, j) * ftemp
             !collision step for f5
               f (i, j, 5) = omega1 * f (i, j, 5) + omega * ftemp + w2 &
              & * ssij
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
End Subroutine collide_srt
!
Subroutine collide_trt (f, c, u, nodetype, tau_a, MagicPara, ss, ly, &
& lx)
    !performs BGK collision substep
      Implicit None
      Real (8), Intent (Inout) :: f (ly, lx, 5)
      Real (8), Intent (In) :: c (ly, lx)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: tau_a (ly, lx)
      Real (8), Intent (In) :: MagicPara
      Real (8), Intent (In) :: nodetype (ly, lx)
      Real (8), Intent (In) :: ss (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: i, j, k
      Real (8) :: w1, w2, f1eq, f2eq, f3eq, f4eq, f5eq, fseq, faeq, fs, &
     & fa, omega_a, omega_s, tau_s, ssij
    !lattice related parameters
      w2 = 1.d0 / 6.d0
      w1 = 2.d0 / 6.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,f1eq,f2eq,f3eq,f4eq,f5eq,fseq,faeq,fs,fa,tau_s,omega_a,omega_s,ssij)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            If (nodetype(i, j) <= 0) Then
               ssij = ss (i, j)
               tau_s = 0.5 + (MagicPara/(tau_a(i, j)-0.5))
               omega_a = 1.d0 / tau_a (i, j)
               omega_s = 1.d0 / tau_s
              !compute equilibrium distribution functions
              !f1
               f1eq = w1 * c (i, j)
              !f2
               f2eq = (1.d0+3.d0*u(1, i, j))
               f2eq = w2 * c (i, j) * f2eq
              !f3
               f3eq = (1.d0+3.d0*u(2, i, j))
               f3eq = w2 * c (i, j) * f3eq
              !f4
               f4eq = (1.d0-3.d0*u(1, i, j))
               f4eq = w2 * c (i, j) * f4eq
              !f5
               f5eq = (1.d0-3.d0*u(2, i, j))
               f5eq = w2 * c (i, j) * f5eq
              !collision step
              !f1
              fseq = (f1eq+f1eq) * 0.5
              fs = (f(i, j, 1)+f(i, j, 1)) * 0.5
              f (i, j, 1) = f (i, j, 1) + omega_s * (fseq-fs) + w1 * &
             & ssij
              !f2
               fseq = (f2eq+f4eq) * 0.5
               faeq = (f2eq-f4eq) * 0.5
               fs = (f(i, j, 2)+f(i, j, 4)) * 0.5
               fa = (f(i, j, 2)-f(i, j, 4)) * 0.5
               f (i, j, 2) = f (i, j, 2) + omega_s * (fseq-fs) + &
              & omega_a * (faeq-fa) + w2 * ssij
              !f4
               faeq = - faeq
               fa = - fa
               f (i, j, 4) = f (i, j, 4) + omega_s * (fseq-fs) + &
              & omega_a * (faeq-fa) + w2 * ssij
              !f3
               fseq = (f3eq+f5eq) * 0.5
               faeq = (f3eq-f5eq) * 0.5
               fs = (f(i, j, 3)+f(i, j, 5)) * 0.5
               fa = (f(i, j, 3)-f(i, j, 5)) * 0.5
               f (i, j, 3) = f (i, j, 3) + omega_s * (fseq-fs) + &
              & omega_a * (faeq-fa) + w2 * ssij
              !f5
               faeq = - faeq
               fa = - fa
               f (i, j, 5) = f (i, j, 5) + omega_s * (fseq-fs) + &
              & omega_a * (faeq-fa) + w2 * ssij
              !f2 becomes f4 and f4 becomes f2
              !f3 becomes f5 and f5 becomes f3
               Do k = 2, 3
                !swap the data
                  fa = f (i, j, k+2)
                  f (i, j, k+2) = f (i, j, k)
                  f (i, j, k) = fa
               End Do
            End If
         End Do
      End Do
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine collide_trt
!
Subroutine collide_diff_vel (f, c, Dr, u, u_adv, nodetype, tau, ss, ly, &
& lx)
    !performs BGK collision
      Implicit None
      Real (8), Intent (Inout) :: f (ly, lx, 5)
      Real (8), Intent (Inout) :: u (2, ly, lx)
      Real (8), Intent (In) :: c (ly, lx)
      Real (8), Intent (In) :: u_adv (2, ly, lx)
      Real (8), Intent (In) :: tau (ly, lx)
      Real (8), Intent (In) :: Dr (ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Real (8), Intent (In) :: ss (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: i, j, k
      Real (8) :: ftemp, omega, omega1, w1, w2, ssij
    !lattice related parameters
      w2 = 1.d0 / 6.d0
      w1 = 2.d0 / 6.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,ftemp,ssij)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            If (nodetype(i, j) <= 0) Then
               omega = (1./tau(i, j))
               omega1 = 1. - omega
			   ssij = ss (i, j)
              !update veolcity field using diffuse velocity formulation
               u (1, i, j) = ud (Dr(i, j), tau(i,j), c(i, j), u_adv(1, i, &
              & j), f(i, j, 2), f(i, j, 4))
               u (1, i, j) = u (1, i, j) + u_adv (1, i, j)
               u (2, i, j) = ud (Dr(i, j), tau(i,j), c(i, j), u_adv(2, i, &
              & j), f(i, j, 3), f(i, j, 5))
               u (2, i, j) = u (2, i, j) + u_adv (2, i, j)
              !compute equillibrium distribution for f1
               ftemp = w1 * c (i, j)
              !collision step for f1
               f (i, j, 1) = omega1 * f (i, j, 1) + omega * ftemp + w1 &
              & * ssij
              !compute equillibrium distribution for f2
               ftemp = (1.d0+3.d0*u(1, i, j))
               ftemp = w2 * c (i, j) * ftemp
              !collision step for f2
               f (i, j, 2) = omega1 * f (i, j, 2) + omega * ftemp + w2 &
              & * ssij
              !compute equillibrium distribution for f3
               ftemp = (1.d0+3.d0*u(2, i, j))
               ftemp = w2 * c (i, j) * ftemp
              !collision step for f3
               f (i, j, 3) = omega1 * f (i, j, 3) + omega * ftemp + w2 &
              & * ssij
              !compute equillibrium distribution for f4
               ftemp = (1.d0-3.d0*u(1, i, j))
               ftemp = w2 * c (i, j) * ftemp
              !collision step for f4
               f (i, j, 4) = omega1 * f (i, j, 4) + omega * ftemp + w2 &
              & * ssij
              !compute equillibrium distribution for f5
               ftemp = (1.d0-3.d0*u(2, i, j))
               ftemp = w2 * c (i, j) * ftemp
              !collision step for f5
               f (i, j, 5) = omega1 * f (i, j, 5) + omega * ftemp + w2 &
              & * ssij
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
Contains
      Real (8) Function ud (Dr, tau, c, ua, fi, fii)
        !function to compute diffuse velocity
         Real (8) :: Dtau, Dr, ua, fi, fii, c, tau
         Dtau = 3.d0 * Dr / tau
         Dtau = Dtau / (1.d0+Dtau)
         ud = Dtau * ((((fi-fii)/((c+1e-30)*(c+1e-30)))*c)-ua)
      End Function ud
End Subroutine collide_diff_vel
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
Subroutine boundary_conditions (f, u, nodetype, tau, interp, location, topbc, &
& topval, bottombc, bottomval, leftbc, leftval, rightbc, rightval, &
&  ly, lx)
    !Imposes boundary conditions on straight walls
    !Accepted inputs
    !topbc    = open,c,flux,periodic,nothing
    !bottombc = open,c,flux,periodic,nothing
    !leftbc   = open,c,flux,periodic,nothing
    !rightbc  = open,c,flux,periodic,nothing
    !location = midway/nodal
      Implicit None
      Real (8), Intent (Inout) :: f (ly, lx, 5)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Character (Len=*), Intent (In) :: topbc
      Character (Len=*), Intent (In) :: bottombc
      Character (Len=*), Intent (In) :: leftbc
      Character (Len=*), Intent (In) :: rightbc
      Character (Len=*), Intent (In) :: location
      Real (8), Intent (In) :: topval
      Real (8), Intent (In) :: bottomval
      Real (8), Intent (In) :: leftval
      Real (8), Intent (In) :: rightval
      Real (8), Intent (In) :: tau (ly, lx)
      Integer, Intent (In) :: interp
      Integer, Intent (In) :: ly, lx
      Integer :: i, j
      Real (8) :: bndval, c0, flux0, ftemp
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
    !nodal
      If (topbc == 'flux' .And. location == 'nodal') Then
         bndval = topval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(1, j) <= 0.) Then
               f (1, j, 5) = OutFFluxBc (bndval, f(1, j, 3), 3, u(2, 1, &
              & j), tau(1, j))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (topbc == 'flux' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,c0,flux0,bndval)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(1, j) <= 0.) Then
               If (interp == 1) Then
                  c0 = f (2, j, 1) + f (2, j, 2) + f (2, j, 3) + f (2, &
                 & j, 4) + f (2, j, 5)
                  flux0 = (1./(2.*tau(2, j))) * c0 * u (2, 2, j) + &
                 & (1.-(1./(2.*tau(2, j)))) * (f(2, j, 3)-f(2, j, 5))
                  bndval = (topval + 0.5d0*flux0)/1.5d0!2. * topval - flux0
                  f (1, j, 5) = OutFFluxBc (bndval, f(1, j, 3), 3, u(2, &
                 & 1, j), tau(1, j))
               Else
                  bndval = topval
               End If
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !bottom boundary
    !nodal
      If (bottombc == 'flux' .And. location == 'nodal') Then
         bndval = bottomval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(ly, j) <= 0.) Then
               f (ly, j, 3) = OutFFluxBc (bndval, f(ly, j, 5), 5, u(2, &
              & ly, j), tau(ly, j))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (bottombc == 'flux' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,c0,flux0,bndval)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(ly, j) <= 0.) Then
               If (interp == 1) Then
                  c0 = f (ly-1, j, 1) + f (ly-1, j, 2) + f (ly-1, j, 3) &
                 & + f (ly-1, j, 4) + f (ly-1, j, 5)
                  flux0 = (1./(2.*tau(ly-1, j))) * c0 * u (2, ly-1, j) &
                 & + (1.-(1./(2.*tau(ly-1, j)))) * (f(ly-1, j, &
                 & 3)-f(ly-1, j, 5))
                  bndval = (bottomval + 0.5d0*flux0)/1.5d0!2. * bottomval - flux0
                  f (ly, j, 3) = OutFFluxBc (bndval, f(ly, j, 5), 5, &
                 & u(2, ly, j), tau(ly, j))
               Else
                  bndval = bottomval
               End If
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !left boundary
    !nodal
      If (leftbc == 'flux' .And. location == 'nodal') Then
         bndval = leftval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, 1) <= 0.) Then
               f (i, 1, 2) = OutFFluxBc (bndval, f(i, 1, 4), 4, u(1, i, &
              & 1), tau(i, 1))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (leftbc == 'flux' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,c0,flux0,bndval)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, 1) <= 0.) Then
               If (interp == 1) Then
                  c0 = f (i, 2, 1) + f (i, 2, 2) + f (i, 2, 3) + f (i, &
                 & 2, 4) + f (i, 2, 5)
                  flux0 = (1./(2.*tau(i, 2))) * c0 * u (1, i, 2) + &
                 & (1.-(1./(2.*tau(i, 2)))) * (f(i, 2, 2)-f(i, 2, 4))
                  bndval = (leftval + 0.5d0*flux0)/1.5d0!2. * leftval - flux0
                  f (i, 1, 2) = OutFFluxBc (bndval, f(i, 1, 4), 4, u(1, &
                 & i, 1), tau(i, 1))
               Else
                  bndval = leftval
               End If
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !right boundary
    !nodal
      If (rightbc == 'flux' .And. location == 'nodal') Then
         bndval = rightval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, lx) <= 0.) Then
               f (i, lx, 4) = OutFFluxBc (bndval, f(i, lx, 2), 2, u(1, &
              & i, lx), tau(i, lx))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (rightbc == 'flux' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,c0,flux0,bndval)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, lx) <= 0.) Then
               If (interp == 1) Then
                  c0 = f (i, lx-1, 1) + f (i, lx-1, 2) + f (i, lx-1, 3) &
                 & + f (i, lx-1, 4) + f (i, lx-1, 5)
                  flux0 = (1./(2.*tau(i, lx-1))) * c0 * u (1, i, &
                 & lx-1) + (1.-(1./(2.*tau(i, lx-1)))) * (f(i, lx-1, &
                 & 2)-f(i, lx-1, 4))
                  bndval =  (rightval + 0.5d0*flux0)/1.5d0!2. * rightval - flux0
                  f (i, lx, 4) = OutFFluxBc (bndval, f(i, lx, 2), 2, &
                 & u(1, i, lx), tau(i, lx))
               Else
                  bndval = rightval
               End If
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
!
    !------------------------
    !c boundary
    !------------------------
    !top boundary
    !nodal
      If (topbc == 'c' .And. location == 'nodal') Then
         bndval = topval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(1, j) <= 0.) Then
               f (1, j, 5) = OutFConcBc (bndval, f(1, j, 1), f(1, j, &
              & 2), f(1, j, 3), f(1, j, 4))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (topbc == 'c' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,c0,bndval)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(1, j) <= 0.) Then
               If (interp == 1) Then
                  c0 = f (2, j, 1) + f (2, j, 2) + f (2, j, 3) + f (2, &
                 & j, 4) + f (2, j, 5)
                  bndval = (topval + 0.5d0*c0)/1.5d0!2.d0 * topval - c0
                  f (1, j, 5) = OutFConcBc (bndval, f(1, j, 1), f(1, j, &
                 & 2), f(1, j, 3), f(1, j, 4))
               Else
                  bndval = topval
               End If
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !bottom boundary
    !nodal
      If (bottombc == 'c' .And. location == 'nodal') Then
         bndval = bottomval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(ly, j) <= 0.) Then
               f (ly, j, 3) = OutFConcBc (bndval, f(ly, j, 1), f(ly, j, &
              & 2), f(ly, j, 4), f(ly, j, 5))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (bottombc == 'c' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,c0,bndval)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(ly, j) <= 0.) Then
               If (interp == 1) Then
                  c0 = f (ly-1, j, 1) + f (ly-1, j, 2) + f (ly-1, j, 3) &
                 & + f (ly-1, j, 4) + f (ly-1, j, 5)
                  bndval = (bottomval + 0.5d0*c0)/1.5d0!2. * bottomval - c0
                  f (ly, j, 3) = OutFConcBc (bndval, f(ly, j, 1), f(ly, &
                 & j, 2), f(ly, j, 4), f(ly, j, 5))
               Else
                  bndval = bottomval
               End If
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !left boundary
    !nodal
      If (leftbc == 'c' .And. location == 'nodal') Then
         bndval = leftval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, 1) <= 0.) Then
               f (i, 1, 2) = OutFConcBc (bndval, f(i, 1, 1), f(i, 1, &
              & 3), f(i, 1, 4), f(i, 1, 5))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (leftbc == 'c' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,c0,bndval)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, 1) <= 0.) Then
               If (interp == 1) Then
                  c0 = f (i, 2, 1) + f (i, 2, 2) + f (i, 2, 3) + f (i, &
                 & 2, 4) + f (i, 2, 5)
                  bndval = (leftval + 0.5d0*c0)/1.5d0
                  f (i, 1, 2) = OutFConcBc (bndval, f(i, 1, 1), f(i, 1, &
                 & 3), f(i, 1, 4), f(i, 1, 5))
               Else
                  bndval = leftval
               End If
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !right boundary
    !nodal
      If (rightbc == 'c' .And. location == 'nodal') Then
         bndval = rightval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, lx) <= 0.) Then
               f (i, lx, 4) = OutFConcBc (bndval, f(i, lx, 1), f(i, lx, &
              & 2), f(i, lx, 3), f(i, lx, 5))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (rightbc == 'c' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,c0,bndval)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, lx) <= 0.) Then
               If (interp == 1) Then
                  c0 = f (i, lx-1, 1) + f (i, lx-1, 2) + f (i, lx-1, 3) &
                 & + f (i, lx-1, 4) + f (i, lx-1, 5)
                  bndval =  (rightval + 0.5d0*c0)/1.5d0
                  f (i, lx, 4) = OutFConcBc (bndval, f(i, lx, 1), f(i, &
                 & lx, 2), f(i, lx, 3), f(i, lx, 5))
               Else
                  bndval = rightval
               End If
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
Contains
      !function to compute boundary conditions
      Real (8) Function OutFFluxBc (bndval, Inf, inIdx, u, tau)
         Implicit None
         Real (8), Intent (In) :: bndval, Inf, u, tau
         Integer, Intent (In) :: inIdx
         Real (8) :: OutF
         Real (8) :: neum, den
         If (inIdx <= 3) Then
            neum = ((1.d0/(2.d0*tau))*(1.d0-3.d0*u)) - 1.d0
            den = ((1.d0/(2.d0*tau))*(-1.d0-3.d0*u)) + 1.d0
            OutF = (bndval-neum*Inf) / den
         Else
            neum = ((1./(2.d0*tau))*(1.d0+3.d0*u)) - 1.d0
            den = ((1./(2.d0*tau))*(-1.d0+3.d0*u)) + 1.d0
            OutF = (bndval-neum*Inf) / den
         End If
         OutFFluxBc = OutF
      End Function OutFFluxBc
      Real (8) Function OutFConcBc (bndval, fin1, fin2, fin3, fin4)
         Implicit None
         Real (8), Intent (In) :: bndval, fin1, fin2, fin3, fin4
         Real (8) :: OutF
         OutF = bndval - fin1 - fin2 - fin3 - fin4
         OutFConcBc = OutF
      End Function OutFConcBc
End Subroutine boundary_conditions
