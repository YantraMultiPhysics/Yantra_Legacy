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
!of multilevel 3D diffusion equation using TRT scheme.
!The relaxation parameter tau can be varying in domain
!
!=======================================================================================
Subroutine compute_macro_var (f, c, flux, poros, nodetype, tau_a, lx, ly, lz)
    !computes the concentration and flux
      Implicit None
      Real (8), Intent (In) :: f (lz, ly, lx, 7)
      Real (8), Intent (Inout) :: c (lz, ly, lx)
      Real (8), Intent (Inout) :: flux (3, lz, ly, lx)
      Real (8), Intent (In) :: tau_a (lz, ly, lx)
      Real (8), Intent (In) :: nodetype (lz, ly, lx)
      Real (8), Intent (In) :: poros (lz, ly, lx)
      Integer, Intent (In) :: lx, ly, lz
      Integer :: i, j, k
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
    !$OMP DO
      Do k = 1, lx
         Do j = 1, ly
            Do i = 1, lz
               c (i, j, k) = 0.d0
               flux (:, i, j, k) = 0.d0
               If (nodetype(i, j, k) <= 0) Then
                  c (i, j, k) = f (i, j, k, 1) + f (i, j, k, 2) + f (i, &
                 & j, k, 3) + f (i, j, k, 4) + f (i, j, k, 5) + f (i, &
                 & j, k, 6) + f (i, j, k, 7)
                  c (i, j, k) = c (i, j, k) / poros (i, j, k)
                  flux (1, i, j, k) = (1-(1./(2.*tau_a(i, j, k)))) * &
                 & (f(i, j, k, 2)-f(i, j, k, 5))
                  flux (2, i, j, k) = (1-(1./(2.*tau_a(i, j, k)))) * &
                 & (f(i, j, k, 3)-f(i, j, k, 6))
                  flux (3, i, j, k) = (1-(1./(2.*tau_a(i, j, k)))) * &
                 & (f(i, j, k, 4)-f(i, j, k, 7))
               End If
            End Do
         End Do
      End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine compute_macro_var
!
Subroutine compute_edf (f, c, cphi, poros, nodetype, lz, ly, lx)
    !computes equillibrium distribution function
      Implicit None
      Real (8), Intent (Out) :: f (lz, ly, lx, 7)
      Real (8), Intent (In) :: c (lz, ly, lx)
      Real (8), Intent (In) :: nodetype (lz, ly, lx)
      Real (8), Intent (In) :: poros (lz, ly, lx)
      Real (8), Intent (In) :: cphi
      Integer, Intent (In) :: lz, ly, lx
      Integer :: i, j, k
      Real (8) :: w
    !lattice related parameters
      w = 1.d0 / 2.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
    !$OMP DO
      Do k = 1, lx
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, k) <= 0) Then
                  !f2
                  f (i, j, k, 2) = w * c (i, j, k) * cphi
                  !f3
                  f (i, j, k, 3) = w * c (i, j, k) * cphi
                  !f4
                  f (i, j, k, 4) = w * c (i, j, k) * cphi
                  !f5
                  f (i, j, k, 5) = w * c (i, j, k) * cphi
                  !f6
                  f (i, j, k, 6) = w * c (i, j, k) * cphi
                  !f7
                  f (i, j, k, 7) = w * c (i, j, k) * cphi
                  !f1
                  f (i, j, k, 1) = poros (i, j, k) * c (i, j, k) - f &
                 & (i, j, k, 2) - f (i, j, k, 3) - f (i, j, k, 4) - f &
                 & (i, j, k, 5) - f (i, j, k, 6) - f (i, j, k, 7)
               Else
                  f (i, j, k, :) = 0.d0
               End If
            End Do
         End Do
      End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine compute_edf
!
Subroutine collide (f, c, cphi, poros, nodetype, tau_a, MagicPara, ss, &
& lz, ly, lx)
    !performs collision step
      Implicit None
      Real (8), Intent (Inout) :: f (lz, ly, lx, 7)
      Real (8), Intent (In) :: c (lz, ly, lx)
      Real (8), Intent (In) :: tau_a (lz, ly, lx)
      Real (8), Intent (In) :: MagicPara
      Real (8), Intent (In) :: nodetype (lz, ly, lx)
      Real (8), Intent (In) :: cphi
      Real (8), Intent (In) :: poros (lz, ly, lx)
      Real (8), Intent (In) :: ss (lz, ly, lx)
      Integer, Intent (In) :: lx, ly, lz
      Integer :: i, j, k, n
      Real (8) :: f1eq, f2eq, f3eq, f4eq, f5eq, f6eq, f7eq, fs, fa, &
     & fseq, faeq, tau_s, omega_s, omega_a, w, ssijk, w1, w2
    !lattice related parameters
      w = 1.d0 / 2.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,n,f1eq, &
    !$OMP f2eq,f3eq,f4eq,f5eq,f6eq,f7eq,fs,fa,fseq,faeq,tau_s,omega_s,omega_a,ssijk,w1,w2)
    !$OMP DO
      Do k = 1, lx
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, k) <= 0) Then
                  ssijk = ss (i, j, k)
                  tau_s = 0.5 + (MagicPara/(tau_a(i, j, k)-0.5))
                  omega_s = 1.d0 / tau_s
                  omega_a = (1.d0/tau_a(i, j, k))
                  !compute equillibrium distribution for f2
                  f2eq = w * c (i, j, k) * cphi
                  !compute equillibrium distribution for f3
                  f3eq = w * c (i, j, k) * cphi
                  !compute equillibrium distribution for f4
                  f4eq = w * c (i, j, k) * cphi
                  !compute equillibrium distribution for f5
                  f5eq = w * c (i, j, k) * cphi
                  !compute equillibrium distribution for f6
                  f6eq = w * c (i, j, k) * cphi
                  !compute equillibrium distribution for f7
                  f7eq = w * c (i, j, k) * cphi
                  !compute equillibrium distribution for f1
                  f1eq = poros (i, j, k) * c (i, j, k) - f2eq - f3eq - &
                 & f4eq - f5eq - f6eq - f7eq
                  !compute weights for  source sink
                  w2 = cphi * w
                  w1 = 1.d0 - 6.d0 * w2
                  !collision step for f1
                  fseq = (f1eq+f1eq) / 2.d0
                  fs = (f(i, j, k, 1)+f(i, j, k, 1)) / 2.d0
                  f (i, j, k, 1) = f (i, j, k, 1) + omega_s * (fseq-fs) &
                 & + w1 * ssijk
                  !collision step for f2
                  fseq = (f2eq+f5eq) / 2.d0
                  faeq = (f2eq-f5eq) / 2.d0
                  fs = (f(i, j, k, 2)+f(i, j, k, 5)) / 2.d0
                  fa = (f(i, j, k, 2)-f(i, j, k, 5)) / 2.d0
                  f (i, j, k, 2) = f (i, j, k, 2) + omega_s * (fseq-fs) &
                 & + omega_a * (faeq-fa) + w2 * ssijk
                  !collision step for f5
                  faeq = - faeq
                  fa = - fa
                  f (i, j, k, 5) = f (i, j, k, 5) + omega_s * (fseq-fs) &
                 & + omega_a * (faeq-fa) + w2 * ssijk
                  !collision step for f3
                  fseq = (f3eq+f6eq) / 2.d0
                  faeq = (f3eq-f6eq) / 2.d0
                  fs = (f(i, j, k, 3)+f(i, j, k, 6)) / 2.d0
                  fa = (f(i, j, k, 3)-f(i, j, k, 6)) / 2.d0
                  f (i, j, k, 3) = f (i, j, k, 3) + omega_s * (fseq-fs) &
                 & + omega_a * (faeq-fa) + w2 * ssijk
                  !collision step for f6
                  faeq = - faeq
                  fa = - fa
                  f (i, j, k, 6) = f (i, j, k, 6) + omega_s * (fseq-fs) &
                 & + omega_a * (faeq-fa) + w2 * ssijk
                  !collision step for f4
                  fseq = (f4eq+f7eq) / 2.d0
                  faeq = (f4eq-f7eq) / 2.d0
                  fs = (f(i, j, k, 4)+f(i, j, k, 7)) / 2.d0
                  fa = (f(i, j, k, 4)-f(i, j, k, 7)) / 2.d0
                  f (i, j, k, 4) = f (i, j, k, 4) + omega_s * (fseq-fs) &
                 & + omega_a * (faeq-fa) + w2 * ssijk
                  !collision step for f7
                  faeq = - faeq
                  fa = - fa
                  f (i, j, k, 7) = f (i, j, k, 7) + omega_s * (fseq-fs) &
                 & + omega_a * (faeq-fa) + w2 * ssijk
                  !f2 becomes f5 and f5 becomes f2
                  !f3 becomes f6 and f6 becomes f3
                  !f4 becomes f7 and f7 becomes f4
                  Do n = 2, 4
                    !swap the data
                     fa = f (i, j, k, n+3)
                     f (i, j, k, n+3) = f (i, j, k, n)
                     f (i, j, k, n) = fa
                  End Do
               End If
            End Do
         End Do
      End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine collide
!
Subroutine stream_and_bounce (f, nodetype, periodicity, lx, ly, lz)
    !performs streaming alonngwith bounce back condition
      Implicit None
      Real (8), Intent (Inout) :: f (lz, ly, lx, 7)
      Real (8), Intent (In) :: nodetype (lz, ly, lx)
      Integer, Intent(In) :: periodicity(3)  
      Integer, Intent (In) :: lx, ly, lz
      Integer :: i, j, k, n, nextX, nextY, nextZ, ex (7), ey (7), ez &
     & (7)
      Real (8) :: ftemp
      ex = (/ 0, 1, 0, 0, - 1, 0, 0 /)
      ey = (/ 0, 0, 1, 0, 0, - 1, 0 /)
      ez = (/ 0, 0, 0, 1, 0, 0, - 1 /)
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,n,nextX,nextY,nextZ,ftemp)
    !$OMP DO
      Do k = 1, lx
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, k) <= 0) Then
                  Do n = 2, 4
                     nextX = k + ex (n)
                     nextY = j - ey (n)
                     nextZ = i - ez (n)
       	          !Apply periodicity
                     If (periodicity(1)>0) Then
                       If (nextX < 1) nextX = lx
                       If (nextX > lx) nextX = 1
                     End If 
                     If (periodicity(2)>0) Then
                       If (nextY > ly) nextY = 1
                       If (nextY < 1) nextY = ly
                     End If
                     If (periodicity(3)>0) Then
                       If (nextZ < 1)  nextZ = lz
                       If (nextZ > lz) nextZ = 1
                     End If 
                     If (nextX > 0 .And. nextX < lx+1 .And. nextY > 0 &
                    & .And. nextY < ly+1 .And. nextZ > 0 .And. nextZ < &
                    & lz+1) Then
                        If (nodetype(nextZ, nextY, nextX) <= 0) Then
                           ftemp = f (nextZ, nextY, nextX, n)
                           f (nextZ, nextY, nextX, n) = f (i, j, k, &
                          & n+3)
                           f (i, j, k, n+3) = ftemp
                        End If
                     End If
                  End Do
               End If
            End Do
         End Do
      End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine stream_and_bounce
!
Subroutine boundary_conditions (f, nodetype, poros, tau_a, &
& interp, location, topbc, topval, bottombc, bottomval, leftbc, leftval, rightbc, &
& rightval, frontbc, frontval, backbc, backval, lx, ly, lz)
    !Imposes boundary conditions on straight walls
    !Accepted inputs
    !topbc    = open,c,flux,periodic,nothing
    !bottombc = open,c,flux,periodic,nothing
    !leftbc   = open,c,flux,periodic,nothing
    !rightbc  = open,c,flux,periodic,nothing
    !frontbc  = open,c,flux,periodic,nothing
    !backbc   = open,c,flux,periodic,nothing
    !location = midway/nodal
      Implicit None
      Real (8), Intent (Inout) :: f (lz, ly, lx, 7)
      Real (8), Intent (In) :: poros (lz, ly, lx)
      Character (Len=*), Intent (In) :: topbc
      Character (Len=*), Intent (In) :: bottombc
      Character (Len=*), Intent (In) :: leftbc
      Character (Len=*), Intent (In) :: rightbc
      Character (Len=*), Intent (In) :: frontbc
      Character (Len=*), Intent (In) :: backbc
      Character (Len=*), Intent (In) :: location
      Real (8), Intent (In) :: topval
      Real (8), Intent (In) :: bottomval
      Real (8), Intent (In) :: leftval
      Real (8), Intent (In) :: rightval
      Real (8), Intent (In) :: frontval
      Real (8), Intent (In) :: backval
      Real (8), Intent (In) :: tau_a (lz, ly, lx)
      Real (8), Intent (In) :: nodetype (lz, ly, lx)
      Integer, Intent (In) :: lz, ly, lx
      Integer, Intent (In) :: interp
      Integer :: i, j
      Real (8) :: bndval, c0, flux0, ftemp
    !--------------------
    !open boundary
    !--------------------
    !top boundary
      If (topbc == 'open') Then
    !$OMP PARALLEL Do DEFAULT(SHARED) PRIVATE(i,j)
         Do j = 1, lx
            Do i = 1, lz
               f (i, 1, j, 6) = f (i, 2, j, 6)
            End Do
         End Do
    !$OMP End PARALLEL Do
      End If
    !bottom boundary
      If (bottombc == 'open') Then
    !$OMP PARALLEL Do DEFAULT(SHARED) PRIVATE(i,j)
         Do j = 1, lx
            Do i = 1, lz
               f (i, ly, j, 3) = f (i, ly-1, j, 3)
            End Do
         End Do
    !$OMP End PARALLEL Do
      End If
    !left boundary
      If (leftbc == 'open') Then
    !$OMP PARALLEL Do DEFAULT(SHARED) PRIVATE(i,j)
         Do j = 1, ly
            Do i = 1, lz
               f (i, j, 1, 2) = f (i, j, 2, 2)
            End Do
         End Do
    !$OMP End PARALLEL Do
      End If
    !right boundary
      If (rightbc == 'open') Then
    !$OMP PARALLEL Do DEFAULT(SHARED) PRIVATE(i,j)
         Do j = 1, ly
            Do i = 1, lz
               f (i, j, lx, 5) = f (i, j, lx-1, 5)
            End Do
         End Do
    !$OMP End PARALLEL Do
      End If
    !back boundary
      If (backbc == 'open') Then
    !$OMP PARALLEL Do DEFAULT(SHARED) PRIVATE(i,j)
         Do j = 1, lx
            Do i = 1, ly
               f (lz, i, j, 4) = f (lz-1, i, j, 4)
            End Do
         End Do
    !$OMP End PARALLEL Do
      End If
    !front boundary
      If (frontbc == 'open') Then
    !$OMP PARALLEL Do DEFAULT(SHARED) PRIVATE(i,j)
         Do j = 1, lx
            Do i = 1, ly
               f (1, i, j, 7) = f (2, i, j, 7)
            End Do
         End Do
    !$OMP End PARALLEL Do
      End If
    !----------------
    !flux boundary
    !----------------
    !top boundary
    !Nodal
      If (topbc == 'flux' .And. location == 'nodal') Then
         bndval = topval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, lz
               If (nodetype(i, 1, j) <= 0.) Then
                  f (i, 1, j, 6) = OutFFluxBc (bndval, f(i, 1, j, 3), &
                 & 3, poros(i, 1, j),  tau_a(i, 1, j))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (topbc == 'flux' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,flux0,bndval)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, lz
               If (nodetype(i, 1, j) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (i, 2, j, 1) + f (i, 2, j, 2) + f (i, 2, j, &
                    & 3) + f (i, 2, j, 4) + f (i, 2, j, 5) + f (i, 2, &
                    & j, 6) + f (i, 2, j, 7)
                     flux0 = (1.-(1./(2.*tau_a(i, 2, j)))) * (f(i, 2, &
                    & j, 3)-f(i, 2, j, 6))
                     bndval = ( topval + 0.5d0 * flux0)/1.5d0
                  Else
                     bndval = topval
                  End If
                  f (i, 1, j, 6) = OutFFluxBc (bndval, f(i, 1, j, 3), &
                 & 3, poros(i, 1, j),  tau_a(i, 1, j))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !bottom boundary
    !nodal
      If (bottombc == 'flux' .And. location == 'nodal') Then
         bndval = bottomval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, lz
               If (nodetype(i, ly, j) <= 0.) Then
                  f (i, ly, j, 3) = OutFFluxBc (bndval, f(i, ly, j, 6), &
                 & 6, poros(i, ly, j),  tau_a(i, ly, j))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (bottombc == 'flux' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,flux0,bndval)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, lz
               If (nodetype(i, ly, j) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (i, ly-1, j, 1) + f (i, ly-1, j, 2) + f (i, &
                    & ly-1, j, 3) + f (i, ly-1, j, 4) + f (i, ly-1, j, &
                    & 5) + f (i, ly-1, j, 6) + f (i, ly-1, j, 7)
                     flux0 = (1.-(1./(2.*tau_a(i, ly-1, j)))) * (f(i, &
                    & ly-1, j, 3)-f(i, ly-1, j, 6))
                     bndval = ( bottomval + 0.5d0*flux0)/1.5d0
                  Else
                     bndval = bottomval
                  End If
                  f (i, ly, j, 3) = OutFFluxBc (bndval, f(i, ly, j, 6), &
                 & 6, poros(i, ly, j),  tau_a(i, ly, j))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !left boundary
    !nodal
      If (leftbc == 'flux' .And. location == 'nodal') Then
         bndval = leftval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, 1) <= 0.) Then
                  f (i, j, 1, 2) = OutFFluxBc (bndval, f(i, j, 1, 5), &
                 & 5, poros(i, j, 1),  tau_a(i, j, 1))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (leftbc == 'flux' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,flux0,bndval)
      !$OMP DO
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, 1) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (i, j, 2, 1) + f (i, j, 2, 2) + f (i, j, 2, &
                    & 3) + f (i, j, 2, 4) + f (i, j, 2, 5) + f (i, j, &
                    & 2, 6) + f (i, j, 2, 7)
                     flux0 = (1.-(1./(2.*tau_a(i, j, 2)))) * (f(i, j, &
                    & 2, 2)-f(i, j, 2, 5))
                     bndval = ( leftval + 0.5d0* flux0)/1.5d0
                  Else
                     bndval = leftval
                  End If
                  f (i, j, 1, 2) = OutFFluxBc (bndval, f(i, j, 1, 5), &
                 & 5, poros(i, j, 1),  tau_a(i, j, 1))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !right boundary
    !nodal
      If (rightbc == 'flux' .And. location == 'nodal') Then
         bndval = rightval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, lx) <= 0.) Then
                  f (i, j, lx, 5) = OutFFluxBc (bndval, f(i, j, lx, 2), &
                 & 2, poros(i, j, lx),  tau_a(i, j, lx))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (rightbc == 'flux' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,flux0,bndval)
      !$OMP DO
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, lx) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (i, j, lx-1, 1) + f (i, j, lx-1, 2) + f (i, &
                    & j, lx-1, 3) + f (i, j, lx-1, 4) + f (i, j, lx-1, &
                    & 5) + f (i, j, lx-1, 6) + f (i, j, lx-1, 7)
                     flux0 = (1.-(1./(2.*tau_a(i, j, lx-1)))) * (f(i, &
                    & j, lx-1, 2)-f(i, j, lx-1, 5))
                     bndval = ( rightval + 0.5* flux0)/1.5d0
                  Else
                     bndval = rightval
                  End If
                  f (i, j, lx, 5) = OutFFluxBc (bndval, f(i, j, lx, 2), &
                 & 2, poros(i, j, lx),  tau_a(i, j, lx))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !back boundary
    !nodal
      If (backbc == 'flux' .And. location == 'nodal') Then
         bndval = backval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, ly
               If (nodetype(lz, i, j) <= 0.) Then
                  f (lz, i, j, 4) = OutFFluxBc (bndval, f(lz, i, j, 7), &
                 & 7, poros(lz, i, j),  tau_a(lz, i, j))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (backbc == 'flux' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,flux0,bndval)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, ly
               If (nodetype(lz, i, j) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (lz-1, i, j, 1) + f (lz-1, i, j, 2) + f &
                    & (lz-1, i, j, 3) + f (lz-1, i, j, 4) + f (lz-1, i, &
                    & j, 5) + f (lz-1, i, j, 6) + f (lz-1, i, j, 7)
                     flux0 = (1.-(1./(2.*tau_a(lz-1, i, j)))) * &
                    & (f(lz-1, i, j, 4)-f(lz-1, i, j, 7))
                     bndval = ( backval + 0.5* flux0)/1.5d0
                  Else
                     bndval = backval
                  End If
                  f (lz, i, j, 4) = OutFFluxBc (bndval, f(lz, i, j, 7), &
                 & 7, poros(lz, i, j),  tau_a(lz, i, j))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !front boundary
    !nodal
      If (frontbc == 'flux' .And. location == 'nodal') Then
         bndval = frontval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, ly
               If (nodetype(1, i, j) <= 0.) Then
                  f (1, i, j, 7) = OutFFluxBc (bndval, f(1, i, j, 4), &
                 & 4, poros(1, i, j),  tau_a(1, i, j))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (frontbc == 'flux' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,bndval,flux0)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, ly
               If (nodetype(1, i, j) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (2, i, j, 1) + f (2, i, j, 2) + f (2, i, j, &
                    & 3) + f (2, i, j, 4) + f (2, i, j, 5) + f (2, i, &
                    & j, 6) + f (2, i, j, 7)
                     flux0 = (1.-(1./(2.*tau_a(2, i, j)))) * (f(2, i, &
                    & j, 4)-f(2, i, j, 7))
                     bndval = ( frontval + 0.5* flux0)/1.5d0
                  Else
                     bndval = frontval
                  End If
                  f (1, i, j, 7) = OutFFluxBc (bndval, f(1, i, j, 4), &
                 & 4, poros(1, i, j),  tau_a(1, i, j))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !------------------------
    !concentration boundary
    !------------------------
    !top boundary
    !nodal
      If (topbc == 'c' .And. location == 'nodal') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,bndval)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, lz
               If (nodetype(i, 1, j) <= 0.) Then
                  bndval = topval * poros (i, 1, j)
                  f (i, 1, j, 6) = OutFConcBc (bndval, f(i, 1, j, 1), &
                 & f(i, 1, j, 2), f(i, 1, j, 3), f(i, 1, j, 4), f(i, 1, &
                 & j, 5), f(i, 1, j, 7))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (topbc == 'c' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,bndval)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, lz
               If (nodetype(i, 1, j) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (i, 2, j, 1) + f (i, 2, j, 2) + f (i, 2, j, &
                    & 3) + f (i, 2, j, 4) + f (i, 2, j, 5) + f (i, 2, &
                    & j, 6) + f (i, 2, j, 7)
                     bndval = (topval +0.5d0*c0/ poros (i, 2, j))/1.5d0 
                     bndval = bndval * poros (i, 1, j)
                  Else
                     bndval = topval * poros (i, 1, j)
                  End If
                  f (i, 1, j, 6) = OutFConcBc (bndval, f(i, 1, j, 1), &
                 & f(i, 1, j, 2), f(i, 1, j, 3), f(i, 1, j, 4), f(i, 1, &
                 & j, 5), f(i, 1, j, 7))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !bottom boundary
    !nodal
      If (bottombc == 'c' .And. location == 'nodal') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j, bndval)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, lz
               If (nodetype(i, ly, j) <= 0.) Then
                  bndval = bottomval * poros (i, ly, j)
                  f (i, ly, j, 3) = OutFConcBc (bndval, f(i, ly, j, 1), &
                 & f(i, ly, j, 2), f(i, ly, j, 4), f(i, ly, j, 5), f(i, &
                 & ly, j, 6), f(i, ly, j, 7))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (bottombc == 'c' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,bndval)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, lz
               If (nodetype(i, ly, j) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (i, ly-1, j, 1) + f (i, ly-1, j, 2) + f (i, &
                    & ly-1, j, 3) + f (i, ly-1, j, 4) + f (i, ly-1, j, &
                    & 5) + f (i, ly-1, j, 6) + f (i, ly-1, j, 7)
                     bndval = (bottomval +0.5d0*c0/ poros (i, ly-1, j))/1.5d0 
                     bndval = bndval * poros (i, ly, j)
                  Else
                     bndval = bottomval * poros (i, ly, j)
                  End If
                  f (i, ly, j, 3) = OutFConcBc (bndval, f(i, ly, j, 1), &
                 & f(i, ly, j, 2), f(i, ly, j, 4), f(i, ly, j, 5), f(i, &
                 & ly, j, 6), f(i, ly, j, 7))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !left boundary
    !nodal
      If (leftbc == 'c' .And. location == 'nodal') Then
!
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j, bndval)
      !$OMP DO
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, 1) <= 0.) Then
                  bndval = leftval * poros (i, j, 1)
                  f (i, j, 1, 2) = OutFConcBc (bndval, f(i, j, 1, 1), &
                 & f(i, j, 1, 3), f(i, j, 1, 4), f(i, j, 1, 5), f(i, j, &
                 & 1, 6), f(i, j, 1, 7))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (leftbc == 'c' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,bndval)
      !$OMP DO
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, 1) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (i, j, 2, 1) + f (i, j, 2, 2) + f (i, j, 2, &
                    & 3) + f (i, j, 2, 4) + f (i, j, 2, 5) + f (i, j, &
                    & 2, 6) + f (i, j, 2, 7)
                     bndval = (leftval +0.5d0*c0/ poros (i, j, 2))/1.5d0 
                     bndval = bndval * poros (i, j, 1)
                  Else
                     bndval = leftval * poros (i, j, 1)
                  End If
                  f (i, j, 1, 2) = OutFConcBc (bndval, f(i, j, 1, 1), &
                 & f(i, j, 1, 3), f(i, j, 1, 4), f(i, j, 1, 5), f(i, j, &
                 & 1, 6), f(i, j, 1, 7))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !right boundary
    !nodal
      If (rightbc == 'c' .And. location == 'nodal') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, lx) <= 0.) Then
                  bndval = rightval * poros (i, j, lx)
                  f (i, j, lx, 5) = OutFConcBc (bndval, f(i, j, lx, 1), &
                 & f(i, j, lx, 2), f(i, j, lx, 3), f(i, j, lx, 4), f(i, &
                 & j, lx, 6), f(i, j, lx, 7))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (rightbc == 'c' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,bndval)
      !$OMP DO
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, lx) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (i, j, lx-1, 1) + f (i, j, lx-1, 2) + f (i, &
                    & j, lx-1, 3) + f (i, j, lx-1, 4) + f (i, j, lx-1, &
                    & 5) + f (i, j, lx-1, 6) + f (i, j, lx-1, 7)
                     bndval = (rightval +0.5d0*c0/ poros (i, j, lx-1))/1.5d0 
                     bndval = bndval * poros (i, j, lx)
                  Else
                     bndval = rightval * poros (i, j, lx)
                  End If
                  f (i, j, lx, 5) = OutFConcBc (bndval, f(i, j, lx, 1), &
                 & f(i, j, lx, 2), f(i, j, lx, 3), f(i, j, lx, 4), f(i, &
                 & j, lx, 6), f(i, j, lx, 7))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !back boundary
    !nodal
      If (backbc == 'c' .And. location == 'nodal') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, ly
               If (nodetype(lz, i, j) <= 0.) Then
                  bndval = backval * poros (lz, i, j)
                  f (lz, i, j, 4) = OutFConcBc (bndval, f(lz, i, j, 1), &
                 & f(lz, i, j, 2), f(lz, i, j, 3), f(lz, i, j, 7), &
                 & f(lz, i, j, 5), f(lz, i, j, 6))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (backbc == 'c' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,c0,bndval)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, ly
               If (nodetype(lz, i, j) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (lz-1, i, j, 1) + f (lz-1, i, j, 2) + f &
                    & (lz-1, i, j, 3) + f (lz-1, i, j, 4) + f (lz-1, i, &
                    & j, 5) + f (lz-1, i, j, 6) + f (lz-1, i, j, 7)
                     bndval = (backval +0.5d0*c0/ poros (lz-1, i, j))/1.5d0 
                     bndval = bndval * poros (lz, i, j)
                  Else
                     bndval = backval * poros (lz, i, j)
                  End If
                  f (lz, i, j, 4) = OutFConcBc (bndval, f(lz, i, j, 1), &
                 & f(lz, i, j, 2), f(lz, i, j, 3), f(lz, i, j, 7), &
                 & f(lz, i, j, 5), f(lz, i, j, 6))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !front boundary
    !nodal
      If (frontbc == 'c' .And. location == 'nodal') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, ly
               If (nodetype(1, i, j) <= 0.) Then
                  bndval = frontval * poros (1, i, j)
                  f (1, i, j, 7) = OutFConcBc (bndval, f(1, i, j, 1), &
                 & f(1, i, j, 2), f(1, i, j, 3), f(1, i, j, 5), f(1, i, &
                 & j, 6), f(1, i, j, 4))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
    !midway
      If (frontbc == 'c' .And. location == 'midway') Then
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,bndval)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, ly
               If (nodetype(1, i, j) <= 0.) Then
                  If (interp == 1) Then
                     c0 = f (2, i, j, 1) + f (2, i, j, 2) + f (2, i, j, &
                    & 3) + f (2, i, j, 4) + f (2, i, j, 5) + f (2, i, &
                    & j, 6) + f (2, i, j, 7)
                     bndval = (frontval +0.5d0*c0/ poros (2, i, j))/1.5d0 
                     bndval = bndval * poros (1, i, j)
                  Else
                     bndval = frontval * poros (1, i, j)
                  End If
                  f (1, i, j, 7) = OutFConcBc (bndval, f(1, i, j, 1), &
                 & f(1, i, j, 2), f(1, i, j, 3), f(1, i, j, 5), f(1, i, &
                 & j, 6), f(1, i, j, 4))
               End If
            End Do
         End Do
      !$OMP END DO
      !$OMP END PARALLEL
      End If
Contains
      Real (8) Function OutFFluxBc (bndval, inF, inIdx, poros,  &
     & tau_a)
         Implicit None
         Real (8), Intent (In) :: bndval, inF, tau_a,  poros
         Integer, Intent (In) :: inIdx
         Real (8) :: OutF, neum, den, a
         a = 0.d0
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
      Real (8) Function OutFConcBc (bndval, fin1, fin2, fin3, fin4, &
     & fin5, fin6)
         Implicit None
         Real (8), Intent (In) :: bndval, fin1, fin2, fin3, fin4, fin5, &
        & fin6
         OutFConcBc = bndval - fin1 - fin2 - fin3 - fin4 - fin5 - fin6
      End Function OutFConcBc
End Subroutine boundary_conditions
!
