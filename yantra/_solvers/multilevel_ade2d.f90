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
!This file contains relavent subroutines for optimized implementation of
!of multilevel 2D advection-diffusion equation using TRT scheme.
!The relaxation parameter tau can be varying in domain
!
!=======================================================================================
Subroutine compute_macro_var (f, c, flux, u, poros, nodetype, tau_a, &
& ly, lx)
    !computes the concentration and flux
      Implicit None
      Real (8), Intent (Inout) :: f (ly, lx, 5)
      Real (8), Intent (Inout) :: c (ly, lx)
      Real (8), Intent (Inout) :: flux (2, ly, lx)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: tau_a (ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Real (8), Intent (In) :: poros (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: i, j
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            c (i, j) = 0.d0
            flux (1, i, j) = 0.d0
            flux (2, i, j) = 0.d0
            If (nodetype(i, j) <= 0) Then
               c (i, j) = f (i, j, 1) + f (i, j, 2) + f (i, j, 3) + f &
              & (i, j, 4) + f (i, j, 5)
               c (i, j) = c (i, j) / poros (i, j)
               flux (1, i, j) = (1./(2.*tau_a(i, j))) * poros (i, j) * &
              & c (i, j) * u (1, i, j) + (1-(1./(2.*tau_a(i, j)))) * &
              & (f(i, j, 2)-f(i, j, 4))
               flux (2, i, j) = (1./(2.*tau_a(i, j))) * poros (i, j) * &
              & c (i, j) * u (1, i, j) + (1-(1./(2.*tau_a(i, j)))) * &
              & (f(i, j, 3)-f(i, j, 5))
            End If
         End Do
      End Do
	!$OMP END DO
    !$OMP END PARALLEL
End Subroutine compute_macro_var
!
Subroutine compute_edf (f, c, u, cphi, poros, nodetype, ly, lx)
      Implicit None
      Real (8), Intent (out) :: f (ly, lx, 5)
      Real (8), Intent (In) :: c (ly, lx)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Real (8), Intent (In) :: poros (ly, lx)
      Real (8), Intent (In) :: cphi
      Integer, Intent (In) :: ly, lx
      Integer :: i, j
      Real (8) :: ftemp, w
    !lattice related parameters
      w = 1.d0 / 2.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,cphi, ftemp)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            If (nodetype(i, j) <= 0) Then
              !f2
               ftemp = (cphi+poros(i, j)*u(1, i, j)+poros(i, j)*u(1, i, &
              & j)*u(1, i, j))
               f (i, j, 2) = w * c (i, j) * ftemp
              !f3
               ftemp = (cphi+poros(i, j)*u(2, i, j)+poros(i, j)*u(2, i, &
              & j)*u(2, i, j))
               f (i, j, 3) = w * c (i, j) * ftemp
              !f4
               ftemp = (cphi-poros(i, j)*u(1, i, j)+poros(i, j)*u(1, i, &
              & j)*u(1, i, j))
               f (i, j, 4) = w * c (i, j) * ftemp
              !f5
               ftemp = (cphi-poros(i, j)*u(2, i, j)+poros(i, j)*u(2, i, &
              & j)*u(2, i, j))
               f (i, j, 5) = w * c (i, j) * ftemp
              !f1
               f (i, j, 1) = poros (i, j) * c (i, j) - (f(i, j, 2)+f(i, &
              & j, 3)+f(i, j, 4)+f(i, j, 5))
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
Subroutine collide (f, c, u, cphi, poros, nodetype, tau_a, MagicPara, &
& ss, ly, lx)
      Implicit None
      Real (8), Intent (Inout) :: f (ly, lx, 5)
      Real (8), Intent (In) :: c (ly, lx)
      Real (8), Intent (In) :: u (2, ly, lx)
      Real (8), Intent (In) :: tau_a (ly, lx)
      Real (8), Intent (In) :: MagicPara
      Real (8), Intent (In) :: nodetype (ly, lx)
      Real (8), Intent (In) :: cphi
      Real (8), Intent (In) :: poros (ly, lx)
      Real (8), Intent (In) :: ss (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: i, j, k
      Real (8) :: w, f1eq, f2eq, f3eq, f4eq, f5eq, fseq, faeq, fs, fa, &
     & omega_a, omega_s, tau_s, ssij, w1, w2
    !lattice related parameters
      w = 1.d0 / 2.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,f1eq,f2eq,f3eq, &
    !$OMP& f4eq,f5eq,fseq,faeq,fs,fa,tau_s,omega_a,omega_s,ssij,w1,w2)
    !$OMP DO
      Do j = 1, lx
         Do i = 1, ly
            If (nodetype(i, j) <= 0) Then
               ssij = ss (i, j)
               tau_s = 0.5 + (MagicPara/(tau_a(i, j)-0.5))
               omega_a = 1.d0 / tau_a (i, j)
               omega_s = 1.d0 / tau_s
               w2 = cphi * w
               w1 = 1.d0 - 2.d0 * cphi
              !f2
               f2eq = (cphi+poros(i, j)*u(1, i, j)+poros(i, j)*u(1, i, &
              & j)*u(1, i, j))
               f2eq = w * c (i, j) * f2eq
              !f3
               f3eq = (cphi+poros(i, j)*u(2, i, j)+poros(i, j)*u(2, i, &
              & j)*u(2, i, j))
               f3eq = w * c (i, j) * f3eq
              !f4
               f4eq = (cphi-poros(i, j)*u(1, i, j)+poros(i, j)*u(1, i, &
              & j)*u(1, i, j))
               f4eq = w * c (i, j) * f4eq
              !f5
               f5eq = (cphi-poros(i, j)*u(2, i, j)+poros(i, j)*u(2, i, &
              & j)*u(2, i, j))
               f5eq = w * c (i, j) * f5eq
!
              !f1
               f1eq = poros (i, j) * c (i, j) - (f2eq+f3eq+f4eq+f5eq)
              !collision step
              !f1
               fseq = (f1eq+f1eq) * 0.5
               fs = (f(i, j, 1)+f(i, j, 1)) * 0.5
               f (i, j, 1) = f (i, j, 1) + omega_s * (fseq-fs) + w1 * &
              & ssij
!
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
End Subroutine collide
!
Subroutine stream_and_bounce (f, nodetype, periodicity, ly, lx)
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
Subroutine boundary_conditions (f, u, nodetype, poros, cphi, tau_a, &
& interp, location, topbc, topval, bottombc, bottomval, leftbc, leftval, rightbc, &
& rightval, ly, lx)
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
      Real (8), Intent (In) :: poros (ly, lx)
      Real (8), Intent (In) :: cphi
      Character (Len=*), Intent (In) :: topbc
      Character (Len=*), Intent (In) :: bottombc
      Character (Len=*), Intent (In) :: leftbc
      Character (Len=*), Intent (In) :: rightbc
      Character (Len=*), Intent (In) :: location
      Real (8), Intent (In) :: topval
      Real (8), Intent (In) :: bottomval
      Real (8), Intent (In) :: leftval
      Real (8), Intent (In) :: rightval
      Real (8), Intent (In) :: tau_a (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer, Intent (In) :: interp
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
               f (1, j, 5) = OutFFluxBc (bndval, f(1, j, 3), 3, &
              & poros(1, j), cphi, u(2, 1, j), tau_a(1, j))
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
                  flux0 = (1./(2.*tau_a(2, j))) * c0 * poros (2, j) * u &
                 & (2, 2, j) + (1.-(1./(2.*tau_a(2, j)))) * (f(2, j, &
                 & 3)-f(2, j, 5))
                  bndval = ( topval + 0.5d0 * flux0)/1.5d0
               Else
                  bndval = topval
               End If
               f (1, j, 5) = OutFFluxBc (bndval, f(1, j, 3), 3, &
              & poros(1, j), cphi, u(2, 1, j), tau_a(1, j))
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
               f (ly, j, 3) = OutFFluxBc (bndval, f(ly, j, 5), 5, &
              & poros(ly, j), cphi, u(2, ly, j), tau_a(ly, j))
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
                  flux0 = (1./(2.*tau_a(ly-1, j))) * c0 * poros (ly-1, &
                 & j) * u (2, ly-1, j) + (1.-(1./(2.*tau_a(ly-1, j)))) &
                 & * (f(ly-1, j, 3)-f(ly-1, j, 5))
                  bndval = ( bottomval + 0.5d0*flux0)/1.5d0
               Else
                  bndval = bottomval
               End If
               f (ly, j, 3) = OutFFluxBc (bndval, f(ly, j, 5), 5, &
              & poros(ly, j), cphi, u(2, ly, j), tau_a(ly, j))
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
               f (i, 1, 2) = OutFFluxBc (bndval, f(i, 1, 4), 4, &
              & poros(i, 1), cphi, u(1, i, 1), tau_a(i, 1))
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
                  flux0 = (1./(2.*tau_a(i, 2))) * c0 * poros (i, 2) * u &
                 & (1, i, 2) + (1.-(1./(2.*tau_a(i, 2)))) * (f(i, 2, &
                 & 2)-f(i, 2, 4))
                  bndval = ( leftval + 0.5d0* flux0)/1.5d0
               Else
                  bndval = leftval
               End If
               f (i, 1, 2) = OutFFluxBc (bndval, f(i, 1, 4), 4, &
              & poros(i, 1), cphi, u(1, i, 1), tau_a(i, 1))
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
               f (i, lx, 4) = OutFFluxBc (bndval, f(i, lx, 2), 2, &
              & poros(i, lx), cphi, u(1, i, lx), tau_a(i, lx))
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
                  flux0 = (1./(2.*tau_a(i, lx-1))) * c0 * poros (i, &
                 & lx-1) * u (1, i, lx-1) + (1.-(1./(2.*tau_a(i, &
                 & lx-1)))) * (f(i, lx-1, 2)-f(i, lx-1, 4))
                  bndval = ( rightval + 0.5* flux0)/1.5d0
               Else
                  bndval = rightval
               End If
               f (i, lx, 4) = OutFFluxBc (bndval, f(i, lx, 2), 2, &
              & poros(i, lx), cphi, u(1, i, lx), tau_a(i, lx))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !------------------------
    !c boundary
    !------------------------
    !top boundary
    !nodal
      If (topbc == 'c' .And. location == 'nodal') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,bndval)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(1, j) <= 0.) Then
               bndval = topval * poros (1, j)
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
                  bndval = (topval +0.5d0*c0/ poros (2, j))/1.5d0 
                  bndval = bndval * poros (1, j)
               Else
                  bndval = topval * poros (1, j)
               End If
               f (1, j, 5) = OutFConcBc (bndval, f(1, j, 1), f(1, j, &
              & 2), f(1, j, 3), f(1, j, 4))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !bottom boundary
    !nodal
      If (bottombc == 'c' .And. location == 'nodal') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j, bndval)
      !$OMP DO
         Do j = 1, lx
            If (nodetype(ly, j) <= 0.) Then
               bndval = bottomval * poros (ly, j)
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
                  bndval = (bottomval +0.5d0*c0/ poros (ly-1, j))/1.5d0 
                  bndval = bndval * poros (ly, j)
               Else
                  bndval = bottomval * poros (ly, j)
               End If
               f (ly, j, 3) = OutFConcBc (bndval, f(ly, j, 1), f(ly, j, &
              & 2), f(ly, j, 4), f(ly, j, 5))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !left boundary
    !nodal
      If (leftbc == 'c' .And. location == 'nodal') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,bndval)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, 1) <= 0.) Then
               bndval = leftval * poros (i, 1)
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
                  bndval =  (leftval +0.5d0*c0/poros (i, 2))/1.5d0   
                  bndval = bndval * poros (i, 1)
               Else
                  bndval = leftval * poros (i, 1)
               End If
               f (i, 1, 2) = OutFConcBc (bndval, f(i, 1, 1), f(i, 1, &
              & 3), f(i, 1, 4), f(i, 1, 5))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !right boundary
    !nodal
      If (rightbc == 'c' .And. location == 'nodal') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,bndval)
      !$OMP DO
         Do i = 1, ly
            If (nodetype(i, lx) <= 0.) Then
               bndval = rightval * poros (i, lx)
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
                  bndval =  (rightval +0.5d0*c0/poros (i, lx-1))/1.5d0 
                  bndval = bndval * poros (i, lx)
               Else
                  bndval = rightval * poros (i, lx)
               End If
               f (i, lx, 4) = OutFConcBc (bndval, f(i, lx, 1), f(i, lx, &
              & 2), f(i, lx, 3), f(i, lx, 5))
            End If
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
Contains
      Real (8) Function OutFFluxBc (bndval, inF, inIdx, poros, cphi, u, &
     & tau_a)
         Implicit None
         Real (8), Intent (In) :: bndval, inF, u, tau_a, cphi, poros
         Integer, Intent (In) :: inIdx
         Real (8) :: OutF, neum, den, a
         a = poros * u / (cphi-u*u)
         If (inIdx <= 3) Then
            neum = ((1.d0/(2.d0*tau_a))*(a-1.d0)) + 1.d0
            den = ((1.d0/(2.d0*tau_a))*(a+1.d0)) - 1.d0
            OutF = (bndval-neum*inF) / den
         Else
            neum = ((1.d0/(2.d0*tau_a))*(a+1.d0)) - 1.d0
            den = ((1.d0/(2.d0*tau_a))*(a-1.d0)) + 1.d0
            OutF = (bndval-neum*inF) / den
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
