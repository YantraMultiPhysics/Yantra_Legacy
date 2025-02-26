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
!Alogrithims to help with geometry update for 3D model
!
!=======================================================================================
!
Subroutine reassign_nodetype(nodetype, lz, ly, lx)
!marks fluid interface node as 0 and solid interface node
!as 1
  Implicit none
  Real (8), Intent (InOut) :: nodetype (lz, ly, lx)
  Integer, Intent (In) :: lz, ly, lx
  Integer :: ex(6),ey(6),ez(6),i,j,k,n,next_i,next_j,next_k,flag
  ex=(/1,0,-1,0,0,0/)
  ey=(/0,1,0,-1,0,0/)
  ez=(/0,0,0,0,1,-1/)
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,next_i,next_j,next_k,n,flag)
  Do k = 1, lx
    Do j = 1, ly
      Do i = 1, lz
        flag = 0
        check: Do n = 1, 6
          next_i = i+ez(n)
          next_j = j+ey(n)
          next_k = k+ex(n)
          If (next_i>0 .AND. next_i<lz+1 .AND. next_j>0 .AND. next_j<ly+1 &
&             .AND. next_k>0 .AND. next_k<lx+1) Then
            If ((nodetype(next_i,next_j,next_k) > 0) .AND. (nodetype(i,j,k) <=0)) Then
              flag = 1
              exit check
            End If
            If ((nodetype(next_i,next_j,next_k) <=0) .AND. (nodetype(i,j,k)  >0)) Then
              flag = 1
              exit check
            End If 
          End If
        End Do check
        If (nodetype(i,j,k)>0) Then
          If (flag == 1) Then
            nodetype(i,j,k)=1
          else
            nodetype(i,j,k)=2
          End If
        Else
          If (flag == 1) Then
            nodetype(i,j,k)=0
          Else 
            nodetype(i,j,k)=-1
          End If
        End If 
      End Do
    End Do
  End Do
  !$OMP END PARALLEL DO
End Subroutine reassign_nodetype
!
Subroutine update_nodetype(nodetype, vol, totvol, frac, lz, ly, lx)
!update node type to fluid if  vol/totvol < frac and 
!vol/totvol > frac update ntype to solid
  implicit none
  Real (8), Intent(Inout) :: nodetype(lz,ly,lx)
  Real (8), Intent(In):: vol(lz,ly,lx)
  Real (8), Intent(In):: totvol
  Real (8), Intent(In):: frac
  Integer, Intent(In) :: lz,ly, lx
  Integer:: i, j , k       
  Real(8) :: phi
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,phi)
  Do k = 1,lx
    Do j = 1,ly
      Do i = 1,lz
        phi = vol(i, j, k)/totvol
        If (nodetype(i, j, k) > 0) Then
          !dissolve
          If (phi < frac) nodetype(i, j, k) = 0
        Else If (nodetype(i, j, k) <= 0) Then
          !precipitate
          If (phi > frac) nodetype(i, j, k) = 1
        End If
      End Do
    End Do
  End Do
  !$OMP END PARALLEL DO
End Subroutine update_nodetype
!
Subroutine update_non_diffusive_phases(dphase, phaseqty, nodetype, vol, mvol, totvol, np, lz, ly, lx)
  implicit none
  Real (8), Intent (Inout) :: phaseqty (np, lz, ly, lx)
  Real (8), Intent (Inout):: vol (lz, ly, lx)
  Real (8), Intent (In) :: dphase(np, lz, ly, lx)
  Real (8), Intent (In) :: nodetype (lz, ly, lx)
  Real (8), Intent (In) :: mvol(np)
  Real (8), Intent (In) :: totvol
  Integer, Intent(In) :: lz,ly, lx, np
  Real (8) :: wi, rdphase, vreq, vavail, vreq_t
  Integer :: ex (6), ey (6), ez(6), i, j, k, p, nextI, nextJ, nextK, ni, itr
  ex = (/ 1, 0, 0,- 1,   0,   0 /)
  ey = (/ 0, 1, 0,  0, - 1,   0 /)
  ez = (/ 0, 0, 1,  0,   0, - 1 /)
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( wi, rdphase, vreq, vavail, vreq_t, i, j, k, p, nextI, nextJ, nextK, ni, itr)
  Do p = 1,np
    Do k = 1, lx
      Do j = 1, ly
        Do i = 1, lz
          If (nodetype(i, j, k)<=0) then
            If (dphase(p,i, j, k) < 0) Then
              !dissolve
              !get weights for neighbours
              wi = 0.d0
              Do ni = 1, 6
                nextI = ez (ni) + i
                nextJ = ey (ni) + j
                nextK = ex (ni) + k
                If (nextI > 0 .And. nextI < lz+1 .And. & 
&                   nextJ > 0 .And. nextJ < ly+1 .And. &
&                   nextK > 0 .And. nextK < lx+1) Then
                  If (nodetype(nextI, nextJ, nextK) > 0 .And. &
&                     phaseqty(p, nextI, nextJ, nextK) > 0) Then
                    wi = wi + 1.d0
                  End If
                End If
              End Do
              If (wi > 0) Then
                wi = 1.d0 / wi
              Else 
                wi = 0.d0
              End If
              !update phase and volumes
              If (phaseqty(p, i, j, k) > 0) Then
                If (phaseqty(p, i, j, k) >-1*dphase(p, i, j, k)) Then
                  !$OMP ATOMIC
                  vol (i, j, k) = vol (i, j, k) + dphase (p, i, j, k) * mvol(p)
                  !$OMP ATOMIC
                  phaseqty (p, i, j, k) = phaseqty (p, i, j, k) + dphase (p, i, j, k)
                Else
                  rdphase = dphase (p, i, j, k) + phaseqty (p, i, j, k)	
                  !$OMP ATOMIC
                  vol (i, j, k) = vol (i, j, k) - phaseqty (p, i, j, k) * mvol(p)
                  !$OMP ATOMIC
                  phaseqty (p, i, j, k) =  0.d0*phaseqty (p, i, j, k)
                  Do ni = 1, 6
                    nextI = ez (ni) + i
                    nextJ = ey (ni) + j
                    nextK = ex (ni) + k
                    If (nextI > 0 .And. nextI < lz+1 .And. &
&                       nextJ > 0 .And. nextJ < ly+1 .And. &
&                       nextK > 0 .And. nextK < lx+1) Then
                      If (nodetype(nextI, nextJ,nextK) > 0 .And. &
&                         phaseqty(p, nextI, nextJ,nextK) > 0) Then
                        !$OMP ATOMIC
                        vol (nextI, nextJ,nextK) = vol (nextI, nextJ,nextK) + &
&                         rdphase * mvol(p) * wi
                        !$OMP ATOMIC
                        phaseqty (p, nextI, nextJ, nextK) = phaseqty &
&                         (p, nextI, nextJ, nextK) + rdphase * wi
                      End If
                    End If
                  End Do
                End If
              Else
                Do ni = 1, 6
                  nextI = ez (ni) + i
                  nextJ = ey (ni) + j
                  nextK = ex (ni) + k
                  If (nextI > 0 .And. nextI < lz+1 .And. &
&                     nextJ > 0 .And. nextJ < ly+1 .And. &
&                     nextK > 0 .And. nextK < lx+1) Then
                    If (nodetype(nextI, nextJ, nextK) > 0 .And. &
&                       phaseqty(p, nextI, nextJ, nextK) > 0) Then
                      !$OMP ATOMIC
                      vol (nextI, nextJ, nextK) = vol (nextI, nextJ, nextK) + &
&                       dphase (p, i, j, k) * mvol(p) * wi
                      !$OMP ATOMIC
                      phaseqty (p, nextI, nextJ, nextK) = phaseqty &
&                       (p, nextI, nextJ, nextK) + dphase (p, i, j, k) * wi
                    End If
                  End If
                End Do
              End If
            Else If (dphase(p, i, j, k) > 0) Then
              !precipitate
              vreq = dphase (p, i, j, k) * mvol(p) !required volume to be filled
              Do itr = 1, 10
                !get weights for neighbours
                wi = 0.d0
                Do ni = 1, 6
                  nextI = ez (ni) + i
                  nextJ = ey (ni) + j
                  nextK = ex (ni) + k
                  If (nextI > 0 .And. nextI < lz+1 .And. &
&                     nextJ > 0 .And. nextJ < ly+1 .And. &
&                     nextK > 0 .And. nextK < lx+1 ) Then
                    vavail = totvol - vol (nextI, nextJ, nextK)
                    If (nodetype(nextI, nextJ, nextK) > 0 .And. vavail > 0) Then
                      wi = wi + 1.d0
                    End If
                  End If
                End Do
                If (wi > 0) Then
                  wi = 1.d0 / wi
                Else 
                  wi = 0.d0
                End If     
                vreq_t = wi * vreq                
                !get neighbours filled
                Do ni = 1, 6
                  nextI = ez (ni) + i
                  nextJ = ey (ni) + j
                  nextK = ex (ni) + k
                  If (nextI > 0 .And. nextI < lz+1 .And. &
&                     nextJ > 0 .And. nextJ < ly+1 .And. &
&                     nextK > 0 .And. nextK < lx+1) Then
                    vavail = totvol - vol (nextI, nextJ, nextK)
                    If (nodetype(nextI, nextJ, nextK) > 0 .And. &
&                       vavail > 0) Then
                      If (vavail > vreq_t) Then
                        vreq = vreq - vreq_t
                        !$OMP ATOMIC
                        vol (nextI, nextJ, nextK) = vol (nextI, &
&                         nextJ, nextK) + vreq_t
                        !$OMP ATOMIC
                        phaseqty (p, nextI, nextJ, nextK) = phaseqty &
&                         (p, nextI, nextJ, nextK) + (vreq_t/mvol(p)) 
                      Else
                        vreq = vreq - vavail
                        !$OMP ATOMIC
                        vol (nextI, nextJ, nextK) = vol(nextI, nextJ, nextK) + vavail
                        !$OMP ATOMIC
                        phaseqty (p, nextI, nextJ, nextK) = phaseqty &
&                         (p, nextI, nextJ, nextK) + (vavail/mvol(p)) 
                      End If
                    End If
                  End If
                End Do
                If (vreq <= 1e-30 .Or. wi <= 0) Exit
              End Do
              If (vreq > 0) Then
                !$OMP ATOMIC
                vol (i, j, k) = vol (i, j, k) + vreq
                !$OMP ATOMIC
                phaseqty (p, i, j, k) = phaseqty (p, i, j, k) + (vreq/mvol(p)) 
              End If
            End If
          End If
        End Do
      End Do
    End Do
  End Do
  !$OMP END PARALLEL DO 
End Subroutine update_non_diffusive_phases
!
subroutine update_diffusive_phases(dphase,phaseqty,nodetype, vol,mvol,np, lz,ly,lx)
   implicit none
   real(8), intent(in):: nodetype(lz,ly,lx)
   real(8), intent(in):: dphase(np,lz,ly,lx)
   real(8), intent(in):: mvol(np)
   real(8), intent(inout):: phaseqty(np,lz,ly,lx)
   real(8), intent(inout):: vol(lz,ly,lx) 
   integer,intent(in):: lz,ly,lx,np
   integer:: i, j, k , p
   vol = 0.d0
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,p)
   Do p = 1,np
     Do k = 1,lx
       Do j = 1,ly
         Do i = 1,lz
           If (nodetype(i,j,k) <=0 .AND. dphase(p,i,j,k) /= 0) then
             phaseqty(p,i,j,k) = phaseqty(p,i,j,k) + dphase(p,i,j,k)
             If (phaseqty(p,i,j,k) < 0) phaseqty(p,i,j,k) = 0.d0
           End If    
             !$OMP ATOMIC
             vol(i,j,k) = vol(i,j,k) + phaseqty(p,i,j,k) * mvol(p) 
         End Do   
       End Do
     End Do
   End Do
   !$OMP END PARALLEL DO
end subroutine update_diffusive_phases
!
subroutine porosity(poros,vol,totvol,nodetype, lz,ly,lx)
  real(8), Intent(Out):: poros(lz,ly,lx)
  real(8), Intent(In)::vol(lz,ly,lx)
  real(8), Intent(In)::nodetype(lz,ly,lx)
  real(8), Intent(In):: totvol
  integer, intent(in):: lz,ly,lx
  integer:: i , j, k
  poros = 0.d0
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
  !$OMP DO
  Do  k = 1,lx
    Do  j = 1,ly
      Do i = 1,lz
        If (nodetype(i, j, k) <=0) Then
          poros(i, j, k) = (1.d0-vol(i, j, k))/totvol
        End If
      End Do
    End Do
  End Do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine porosity
!
subroutine phaseqty_for_phrqc(phaseqty_phrqc, nodetype, phaseqty, np, lz, ly, lx)
  implicit none
  Real (8), Intent (Out) :: phaseqty_phrqc (np, lz,ly, lx)
  Real (8), Intent (In) :: nodetype (lz, ly, lx)
  Real (8), Intent (In) :: phaseqty (np, lz, ly, lx)
  Real (8):: wi 
  Integer, Intent (In) :: lz, ly, lx, np
  Integer :: ex(6),ey(6), ez(6), ni,nextI,nextJ, nextK, i, j, k, p
  ex = (/ 1, 0, 0,- 1,   0,   0 /)
  ey = (/ 0, 1, 0,  0, - 1,   0 /)
  ez = (/ 0, 0, 1,  0,   0, - 1 /)
  phaseqty_phrqc = 0.d0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ni,nextI,nextJ,nextK,i,j,k,p)
  Do p = 1, np
    Do k = 1, lx
      Do j = 1, ly
        Do i = 1, lz
          wi = 0.d0
          If (nodetype(i, j, k) > 0 .AND. phaseqty(p, i, j, k) > 0) Then 
            Do ni = 1, 6
              nextI = ez (ni) + i
              nextJ = ey (ni) + j
              nextK = ex (ni) + k
              If (nextI > 0 .And. nextI < lz+1 .And. &
&                 nextJ > 0 .And. nextJ < ly+1 .And. &
&                 nextK > 0 .And. nextK < lx+1) Then
                If (nodetype(nextI, nextJ, nextK) <= 0) Then
                  wi = wi + 1.d0
                End If
              End If
            End Do
            wi = 1.d0/wi
            !Distribute values
            Do ni = 1, 6
              nextI = ez (ni) + i
              nextJ = ey (ni) + j
              nextK = ex (ni) + k
              If (nextI > 0 .And. nextI < lz+1 .And. &
&                 nextJ > 0 .And. nextJ < ly+1 .And. &
&                 nextK > 0 .And. nextK < lx+1) Then
                If (nodetype(nextI, nextJ, nextK) <= 0) Then
                  phaseqty_phrqc(p,nextI,nextJ, nextK) = phaseqty_phrqc(p,nextI,nextJ,nextK) + &
                  wi * phaseqty(p, i, j, k)
                End If
              End If
            End Do
          Else If (nodetype(i,j,k) <=0) Then
            phaseqty_phrqc(p,i,j,k)=phaseqty_phrqc(p,i,j,k)+phaseqty(p,i,j,k)
          End If
        End Do
      End Do
    End Do
  End Do
  !$OMP END PARALLEL DO
End subroutine phaseqty_for_phrqc
