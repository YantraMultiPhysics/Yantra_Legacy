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
!of pore-scale 3D  diffusion equation with varying relaxation parameters
!
!=======================================================================================
!
Subroutine compute_macro_var (f, c, flux, nodetype, tau, lx, ly, lz)
    !computes the concentration and flux
      Implicit None
      Real (8), Intent (In) :: f (lz, ly, lx, 7)
      Real (8), Intent (Inout) :: c (lz, ly, lx)
      Real (8), Intent (Inout) :: flux (3, lz, ly, lx)
      Real (8), Intent (In) :: tau (lz, ly, lx)
      Real (8), Intent (In) :: nodetype (lz, ly, lx)
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
                  flux (1, i, j, k) = (1-(1./(2.*tau(i, j, k)))) * &
                 & (f(i, j, k, 2)-f(i, j, k, 5))
                  flux (2, i, j, k) = (1-(1./(2.*tau(i, j, k)))) * &
                 & (f(i, j, k, 3)-f(i, j, k, 6))
                  flux (3, i, j, k) = (1-(1./(2.*tau(i, j, k)))) * &
                 & (f(i, j, k, 4)-f(i, j, k, 7))
               End If
            End Do
         End Do
      End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine compute_macro_var
!
Subroutine compute_edf (f, c, nodetype, lx, ly, lz)
    !computes equillibrium distribution function
      Implicit None
      Real (8), Intent (out) :: f (lz, ly, lx, 7)
      Real (8), Intent (In) :: c (lz, ly, lx)
      Real (8), Intent (In) :: nodetype (lz, ly, lx)
      Integer, Intent (In) :: lz, ly, lx
      Integer :: i, j, k
      Real (8) :: ftemp, omega, omega1, w1, w2
    !lattice related parameters
      w2 = 1.d0 / 7.d0
      w1 = 1.d0 / 7.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,ftemp)
    !$OMP DO
      Do k = 1, lx
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, k) <= 0) Then
                  f (i, j, k, 1) = w1 * c (i, j, k)
                  f (i, j, k, 2) = w2 * c (i, j, k)
                  f (i, j, k, 3) = w2 * c (i, j, k)
                  f (i, j, k, 4) = w2 * c (i, j, k)
                  f (i, j, k, 5) = w2 * c (i, j, k)
                  f (i, j, k, 6) = w2 * c (i, j, k)
                  f (i, j, k, 7) = w2 * c (i, j, k)
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
Subroutine collide_srt (f, c, nodetype, tau, ss, lx, ly, lz)
    !performs collision step
      Implicit None
      Real (8), Intent (Inout) :: f (lz, ly, lx, 7)
      Real (8), Intent (In) :: c (lz, ly, lx)
      Real (8), Intent (In) :: tau (lz, ly, lx)
      Real (8), Intent (In) :: nodetype (lz, ly, lx)
      Real (8), Intent (In) :: ss (lz, ly, lx)
      Integer, Intent (In) :: lx, ly, lz
      Integer :: i, j, k, n
      Real (8) :: ftemp, omega, omega1, w1, w2, ssijk
      w2 = 1.d0 / 7.d0
      w1 = 1.d0 / 7.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,n,ftemp,omega,omega1)
    !$OMP DO
      Do k = 1, lx
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, k) <= 0) Then
                  ssijk = ss (i, j, k)
                  omega = (1./tau(i, j, k))
                  omega1 = 1. - omega
                  ftemp = w1 * c (i, j, k)
                  f (i, j, k, 1) = omega1 * f (i, j, k, 1) + omega * &
                 & ftemp + w1 * ssijk
                  ftemp = w2 * c (i, j, k)
                  f (i, j, k, 2) = omega1 * f (i, j, k, 2) + omega * &
                 & ftemp + w2 * ssijk
                  f (i, j, k, 3) = omega1 * f (i, j, k, 3) + omega * &
                 & ftemp + w2 * ssijk
                  f (i, j, k, 4) = omega1 * f (i, j, k, 4) + omega * &
                 & ftemp + w2 * ssijk
                  f (i, j, k, 5) = omega1 * f (i, j, k, 5) + omega * &
                 & ftemp + w2 * ssijk
                  f (i, j, k, 6) = omega1 * f (i, j, k, 6) + omega * &
                 & ftemp + w2 * ssijk
                  f (i, j, k, 7) = omega1 * f (i, j, k, 7) + omega * &
                 & ftemp + w2 * ssijk
                  Do n = 2, 4
                     ftemp = f (i, j, k, n+3)
                     f (i, j, k, n+3) = f (i, j, k, n)
                     f (i, j, k, n) = ftemp
                  End Do
               End If
            End Do
         End Do
      End Do
    !$OMP END DO
    !$OMP END PARALLEL
End Subroutine collide_srt
!
Subroutine collide_trt (f, c, nodetype, tau_a, MagicPara, ss, lx, ly, &
& lz)
      Implicit None
      Real (8), Intent (Inout) :: f (lz, ly, lx, 7)
      Real (8), Intent (In) :: c (lz, ly, lx)
      Real (8), Intent (In) :: tau_a (lz, ly, lx)
      Real (8), Intent (In) :: MagicPara
      Real (8), Intent (In) :: nodetype (lz, ly, lx)
      Real (8), Intent (In) :: ss (lz, ly, lx)
      Integer, Intent (In) :: lx, ly, lz
      Integer :: i, j, k, n
      Real (8) :: fs, fa, fseq, faeq, tau_s, omega_s, omega_a, w1, w2, &
     & ssijk
    !lattice related parameters
      w2 = 1.d0 / 7.d0
      w1 = 1.d0 / 7.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k, n, &
    !$OMP& fs,fa,fseq,faeq,tau_s,omega_s,omega_a,ssijk)
    !$OMP DO
      Do k = 1, lx
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, k) <= 0) Then
                  ssijk = 0.d0
                  ssijk = ss (i, j, k)
                  tau_s = 0.5 + (MagicPara/(tau_a(i, j, k)-0.5))
                  omega_s = 1.d0 / tau_s
                  omega_a = (1.d0/tau_a(i, j, k))
                  !compute equillibrium distribution for f1
                  fseq = w1 * c (i, j, k)!faeq=0
                  !collision step for f1
                  fs = f (i, j, k, 1)
                  f (i, j, k, 1) = f (i, j, k, 1) + omega_s * (fseq-fs) &
                 & + w1 * ssijk
                  !compute equillibrium distribution for f2-f7
                  fseq = w2 * c (i, j, k)!faeq=0
                  !collision step for f2
                  fs = (f(i, j, k, 2)+f(i, j, k, 5)) / 2.d0
                  fa = (f(i, j, k, 2)-f(i, j, k, 5)) / 2.d0
                  f (i, j, k, 2) = f (i, j, k, 2) + omega_s * (fseq-fs) &
                 & - omega_a * fa + w2 * ssijk
                  !collision step for f5
                  fa = - fa
                  f (i, j, k, 5) = f (i, j, k, 5) + omega_s * (fseq-fs) &
                 & - omega_a * fa + w2 * ssijk
                  !collision step for f3
                  fs = (f(i, j, k, 3)+f(i, j, k, 6)) / 2.d0
                  fa = (f(i, j, k, 3)-f(i, j, k, 6)) / 2.d0
                  f (i, j, k, 3) = f (i, j, k, 3) + omega_s * (fseq-fs) &
                 & - omega_a * fa + w2 * ssijk
                  !collision step for f6
                  fa = - fa
                  f (i, j, k, 6) = f (i, j, k, 6) + omega_s * (fseq-fs) &
                 & - omega_a * fa + w2 * ssijk
                  !collision step for f4
                  fs = (f(i, j, k, 4)+f(i, j, k, 7)) / 2.d0
                  fa = (f(i, j, k, 4)-f(i, j, k, 7)) / 2.d0
                  f (i, j, k, 4) = f (i, j, k, 4) + omega_s * (fseq-fs) &
                 & - omega_a * fa + w2 * ssijk
                  !collision step for f7
                  fa = - fa
                  f (i, j, k, 7) = f (i, j, k, 7) + omega_s * (fseq-fs) &
                 & - omega_a * fa + w2 * ssijk
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
End Subroutine collide_trt
!
Subroutine stream_and_bounce (f, nodetype, periodicity, lx, ly, lz)
    !performs streaming alonngwith bounce back condition
      Implicit None
      Real (8), Intent (Inout) :: f (lz, ly, lx, 7)
      Real (8), Intent (Inout) :: nodetype (lz, ly, lx)
      Integer, Intent(In) :: periodicity(3)  
      Integer, Intent (In) :: lx, ly, lz
      Integer :: i, j, k, n, nextX, nextY, nextZ, ex (7), ey (7), ez &
     & (7)
      Real (8) :: ftemp
      ex = (/ 0, 1, 0, 0, - 1, 0, 0 /)
      ey = (/ 0, 0, 1, 0, 0, - 1, 0 /)
      ez = (/ 0, 0, 0, 1, 0, 0, - 1 /)
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,n,nextX,nextY,nextZ,ftemp)
    !$OMP Do
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
Subroutine boundary_conditions (f, nodetype, tau, interp, location, topbc, &
& topval, bottombc, bottomval, leftbc, leftval, rightbc, rightval, &
& frontbc, frontval, backbc, backval, lx, ly, lz)
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
      Real (8), Intent (In) :: tau (lz, ly, lx)
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
                 & 3, tau(i, 1, j))
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
                     flux0 = (1.-(1./(2.*tau(i, 2, j)))) * (f(i, 2, j, &
                    & 3)-f(i, 2, j, 6))
                     bndval = ( topval + 0.5d0 * flux0)/1.5d0
                  Else
                     bndval = topval
                  End If
                  f (i, 1, j, 6) = OutFFluxBc (bndval, f(i, 1, j, 3), &
                 & 3, tau(i, 1, j))
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
                 & 6, tau(i, ly, j))
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
                     flux0 = (1.-(1./(2.*tau(i, ly-1, j)))) * (f(i, &
                    & ly-1, j, 3)-f(i, ly-1, j, 6))
                     bndval = ( bottomval + 0.5d0*flux0)/1.5d0
                  Else
                     bndval = bottomval
                  End If
                  f (i, ly, j, 3) = OutFFluxBc (bndval, f(i, ly, j, 6), &
                 & 6, tau(i, ly, j))
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
                 & 5, tau(i, j, 1))
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
                     flux0 = (1.-(1./(2.*tau(i, j, 2)))) * (f(i, j, 2, &
                    & 2)-f(i, j, 2, 5))
                     bndval = ( leftval + 0.5d0* flux0)/1.5d0
                  Else
                     bndval = leftval
                  End If
                  f (i, j, 1, 2) = OutFFluxBc (bndval, f(i, j, 1, 5), &
                 & 5, tau(i, j, 1))
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
                 & 2, tau(i, j, lx))
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
                     flux0 = (1.-(1./(2.*tau(i, j, lx-1)))) * (f(i, j, &
                    & lx-1, 2)-f(i, j, lx-1, 5))
                     bndval = ( rightval + 0.5* flux0)/1.5d0
                  Else
                     bndval = rightval
                  End If
                  f (i, j, lx, 5) = OutFFluxBc (bndval, f(i, j, lx, 2), &
                 & 2, tau(i, j, lx))
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
                 & 7, tau(lz, i, j))
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
                     flux0 = (1.-(1./(2.*tau(lz-1, i, j)))) * (f(lz-1, &
                    & i, j, 4)-f(lz-1, i, j, 7))
                     bndval = ( backval + 0.5* flux0)/1.5d0
                  Else
                     bndval = backval
                  End If
                  f (lz, i, j, 4) = OutFFluxBc (bndval, f(lz, i, j, 7), &
                 & 7, tau(lz, i, j))
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
                 & 4, tau(1, i, j))
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
                     flux0 = (1.-(1./(2.*tau(lz-1, i, j)))) * (f(2, i, &
                    & j, 4)-f(2, i, j, 7))
                     bndval = ( frontval + 0.5* flux0)/1.5d0
                  Else
                     bndval = frontval
                  End If
                  f (1, i, j, 7) = OutFFluxBc (bndval, f(1, i, j, 4), &
                 & 4, tau(1, i, j))
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
         bndval = topval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, lz
               If (nodetype(i, 1, j) <= 0.) Then
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
                     bndval = (topval +0.5d0*c0)/1.5d0
                  Else
                     bndval = topval
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
         bndval = bottomval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, lz
               If (nodetype(i, ly, j) <= 0.) Then
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
                     bndval = (bottomval +0.5d0*c0)/1.5d0
                  Else
                     bndval = bottomval
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
         bndval = leftval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, 1) <= 0.) Then
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
                     bndval = (leftval +0.5d0*c0)/1.5d0
                  Else
                     bndval = leftval
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
         bndval = rightval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, ly
            Do i = 1, lz
               If (nodetype(i, j, lx) <= 0.) Then
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
                     bndval = (rightval +0.5d0*c0)/1.5d0
                  Else
                     bndval = rightval
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
         bndval = backval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, ly
               If (nodetype(lz, i, j) <= 0.) Then
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
                     bndval = (backval +0.5d0*c0)/1.5d0
                  Else
                     bndval = backval
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
         bndval = frontval
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
      !$OMP DO
         Do j = 1, lx
            Do i = 1, ly
               If (nodetype(1, i, j) <= 0.) Then
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
                     bndval = (frontval +0.5d0*c0)/1.5d0
                  Else
                     bndval = frontval
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
      Real (8) Function OutFFluxBc (bndval, inF, inIdx, tau)
         Implicit None
         Real (8), Intent (In) :: bndval, inF, tau
         Integer, Intent (In) :: inIdx
         Real (8) :: neum, den
         If (inIdx <= 4) Then
            neum = ((1./(2.d0*tau))*(-1.)) + 1.
            den = ((1./(2.d0*tau))*(1.)) - 1.
            OutFFluxBc = (bndval-neum*inF) / den
         Else
            neum = ((1./(2.d0*tau))*(1.)) - 1.
            den = ((1./(2.d0*tau))*(-1.)) + 1.
            OutFFluxBc = (bndval-neum*inF) / den
         End If
      End Function OutFFluxBc
      Real (8) Function OutFConcBc (bndval, fin1, fin2, fin3, fin4, &
     & fin5, fin6)
         Implicit None
         Real (8), Intent (In) :: bndval, fin1, fin2, fin3, fin4, fin5, &
        & fin6
         OutFConcBc = bndval - fin1 - fin2 - fin3 - fin4 - fin5 - fin6
      End Function OutFConcBc
End Subroutine boundary_conditions
