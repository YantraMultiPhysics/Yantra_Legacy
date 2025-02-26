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
!Alogrithims to help with geometry update for 2D model
!
!=======================================================================================
!
Subroutine reassign_nodetype(nodetype, ly, lx)
!marks fluid interface node as 0 and solid interface node
!as 1
      Implicit none
      Real (8), Intent (InOut) :: nodetype (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: ex(4),ey(4),i,j,k,next_i,next_j,flag
      ex=(/1,0,-1,0/)
      ey=(/0,1,0,-1/)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k,flag)
      Do j = 1, lx
        Do i = 1, ly
          flag = 0
          check: Do k = 1, 4
            next_i = i+ey(k)
            next_j = j+ex(k)
            If (next_i>0 .AND. next_i<ly+1 .AND. next_j>0 .AND. next_j<lx+1) Then
              If ((nodetype(next_i,next_j) > 0) .AND. (nodetype(i,j) <=0)) Then
                flag = 1
                exit check
              End If
              If ((nodetype(next_i,next_j) <=0) .AND. (nodetype(i,j)  >0)) Then
                flag = 1
                exit check
              End if
            End if 
          End Do check
          If (nodetype(i,j)>0) Then
            If (flag == 1) Then
              nodetype(i,j)=1
            else
              nodetype(i,j)=2
            End If
          Else
            If (flag == 1) Then
              nodetype(i,j)=0
            else  
              nodetype(i,j)=-1
            End If
          End If 
        End Do
      End Do
      !$OMP END PARALLEL DO
End Subroutine reassign_nodetype
!
Subroutine update_nodetype(nodetype, vol, totvol, frac, ly, lx)
!update node type to fluid if  vol/totvol < frac and 
!vol/totvol > frac update nodetype to solid
    Implicit None
    Real (8), Intent(Inout) :: nodetype(ly,lx)
    Real (8), Intent(In) :: vol(ly,lx)
    Real (8), Intent(In) :: totvol
    Real (8), Intent(In) :: frac
    Integer, Intent(In) :: ly, lx
    Integer:: i, j        
    Real(8) :: phi
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,phi)
    Do j = 1,lx
      Do i = 1,ly
        phi = vol(i, j)/totvol
        If (nodetype(i, j) > 0) Then
          !dissolve
          If (phi < frac) nodetype(i, j) = 0
        Else If (nodetype(i, j) <= 0) Then
          !precipitate
          If (phi > frac) nodetype(i, j) = 1
        End If
      End Do
    End Do
    !$OMP END PARALLEL DO
End Subroutine update_nodetype
!
Subroutine update_non_diffusive_phases(dphase, phaseqty, nodetype, vol, mvol, &
& totvol, np, ly, lx)
  implicit none
  Real (8), Intent (Inout) :: phaseqty (np, ly, lx)
  Real (8), Intent (Inout) :: vol (ly, lx)
  Real (8), Intent (In) :: dphase (np, ly, lx)
  Real (8), Intent (In) :: nodetype (ly, lx)
  Real (8), Intent (In) :: mvol(np) 
  Real(8), Intent (In) :: totvol
  Integer, Intent(In) :: ly, lx, np
  Real (8) :: wi, rdphase, vreq, vavail, vreq_t
  Integer :: ex (4), ey (4), i, j, p, nextI, nextJ, ni, itr
  ex = (/ 1, 0, -1, 0 /)
  ey = (/ 0, 1, 0, -1 /)
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(wi, rdphase, vreq, vavail, vreq_t, i, j, p, nextI, nextJ, ni, itr )
  Do p = 1,np
    Do j = 1, lx
      Do  i = 1, ly
        If (nodetype(i, j)<=0) then
          If (dphase(p, i, j) < 0) Then
            !dissolve
            !get weights for neighbours
            wi = 0.d0
            Do ni = 1, 4
              nextI = ey (ni) + i
              nextJ = ex (ni) + j
              If (nextI > 0 .And. nextI < ly + 1 .And. nextJ > 0 .And. nextJ &
&                 < lx + 1) Then
                If (nodetype(nextI, nextJ) > 0 .And. phaseqty(p, nextI, &
&                   nextJ) > 0) Then
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
            If (phaseqty(p, i, j) > 0) Then
              If (phaseqty(p, i, j) >-1.*dphase(p, i, j)) Then
                !$OMP ATOMIC
                vol (i, j) = vol (i, j) + dphase (p, i, j) * mvol(p)
                !$OMP ATOMIC
                phaseqty (p, i, j) = phaseqty (p, i, j) + dphase (p, i, j)
              Else
                rdphase = dphase (p, i, j) + phaseqty (p, i, j)
                !$OMP ATOMIC
                vol (i, j) = vol (i, j) - phaseqty (p, i, j) * mvol(p)
                !$OMP ATOMIC
                phaseqty (p, i, j) = 0.d0 * phaseqty (p, i, j)
                Do ni = 1, 4
                  nextI = ey (ni) + i
                  nextJ = ex (ni) + j
                  If (nextI > 0 .And. nextI < ly +1 .And. nextJ > 0 .And. &
&                     nextJ < lx+1) Then
                    If (nodetype(nextI, nextJ) > 0 .And. &
&                       phaseqty(p,nextI, nextJ) > 0) Then
                      !$OMP ATOMIC
                      vol (nextI, nextJ) = vol (nextI, nextJ) + &
&                       rdphase * mvol(p) * wi
                      !$OMP ATOMIC
                      phaseqty (p, nextI, nextJ) = phaseqty &
&                       (p, nextI, nextJ) + rdphase * wi 
                    End If
                  End If
                End Do
              End If
            Else
              Do ni = 1, 4
                nextI = ey (ni) + i
                nextJ = ex (ni) + j
                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1) Then
                  If (nodetype(nextI, nextJ) > 0 .And. &
&                     phaseqty(p, nextI, nextJ) > 0) Then
                    !$OMP ATOMIC
                    vol (nextI, nextJ) = vol (nextI, nextJ) + &
&                     dphase (p, i, j) * mvol(p) * wi
                    !$OMP ATOMIC
                    phaseqty (p, nextI, nextJ) = phaseqty &
&                     (p, nextI, nextJ) + dphase (p, i, j) * wi
                  End If
                End If
              End Do
            End If
          Else If (dphase(p, i, j) > 0) Then
            !precipitate
            vreq = dphase (p, i, j) * mvol(p) !required volume to be filled
            Do itr = 1, 10
              !get weights for neighbours
              wi = 0.d0
              Do ni = 1, 4
                nextI = ey (ni) + i
                nextJ = ex (ni) + j
                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1) Then
                  vavail = totvol - vol (nextI, nextJ)
                  If (nodetype(nextI, nextJ) > 0 .And. vavail > 0) Then
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
              Do ni = 1, 4
                nextI = ey (ni) + i
                nextJ = ex (ni) + j
                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1) Then
                  vavail = totvol - vol (nextI, nextJ)
                  If (nodetype(nextI, nextJ) > 0 .And. vavail >  0) Then
                    If (vavail > vreq_t) Then
                      vreq = vreq - vreq_t
                      !$OMP ATOMIC
                      vol (nextI, nextJ) = vol (nextI, &
&                       nextJ) + vreq_t
                      !$OMP ATOMIC
                      phaseqty (p, nextI, nextJ) = phaseqty &
&                       (p, nextI, nextJ) + (vreq_t/mvol(p)) 
                    Else
                      vreq = vreq - vavail
                      !$OMP ATOMIC
                      vol (nextI, nextJ) = vol(nextI, nextJ) + vavail
                      !$OMP ATOMIC
                      phaseqty (p, nextI, nextJ) = phaseqty &
&                       (p, nextI, nextJ) + (vavail/mvol(p)) 
                    End If
                  End If
                End If
              End Do
              If (vreq <= 1e-30 .Or. wi <= 0) Exit
            End Do
            If (vreq > 0) Then
              !$OMP ATOMIC
              vol (i, j) = vol (i, j) + vreq
              !$OMP ATOMIC
              phaseqty (p, i, j) = phaseqty (p, i, j) + (vreq/mvol(p)) 
            End If
          End If
        End If
      End Do
    End Do
  End Do
  !$OMP END PARALLEL DO 
End Subroutine
!
subroutine update_diffusive_phases(dphase,phaseqty,nodetype,vol,mvol,np,ly,lx)
   Implicit None
   real(8),intent(in):: nodetype(ly,lx)
   real(8),intent(in):: dphase(np,ly,lx)
   real(8),intent(in):: mvol(np)
   real(8),intent(inout):: phaseqty(np,ly,lx)
   real(8),intent(inout):: vol(ly,lx)
   integer,intent(in):: np,ly,lx
   integer:: i, j, p
   vol = 0.d0
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,p)
   Do p = 1,np
     Do j = 1,lx
       Do i = 1,ly
         If (nodetype(i,j) <=0 .AND. dphase(p,i,j) /= 0) then
           phaseqty(p,i,j) = phaseqty(p,i,j) + dphase(p,i,j)
           If (phaseqty(p,i,j) < 0) phaseqty(p,i,j)=0. 
         End If       
           !$OMP ATOMIC
           vol(i,j) = vol(i,j) + phaseqty(p,i,j) * mvol(p)
       End Do
     End Do
   End Do
   !$OMP END PARALLEL DO
end subroutine update_diffusive_phases
!
subroutine porosity(poros,vol,totvol,nodetype,ly,lx)
   Implicit None
   Real(8), Intent(Out):: poros(ly,lx)
   Real(8), Intent(In):: vol(ly,lx)
   Real(8), Intent(In):: nodetype(ly,lx)
   Real(8), Intent(In):: totvol
   Integer, Intent(In):: ly,lx
   Integer:: i,j
   poros = 0.d0
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
   Do j = 1,lx
     Do i = 1,ly 
        If (nodetype(i,j) <=0) Then
            poros(i, j) = (1.-vol(i,j))/totvol
        End If
     End Do
   End Do
   !$OMP END PARALLEL DO
end subroutine porosity
!
subroutine phaseqty_for_phrqc(phaseqty_phrqc, nodetype, phaseqty,np, ly, lx)
         Implicit None
         Real (8), Intent (Out) :: phaseqty_phrqc (np, ly, lx)
         Real (8), Intent (In) :: nodetype (ly, lx)
         Real (8), Intent (In) :: phaseqty (np, ly, lx)
         Integer, Intent (In) :: np,ly, lx
         Real (8):: wi 
         Integer ::ex(4),ey(4), ni,nextI,nextJ,i,j,p
         ex = (/ 1, 0, -1, 0 /)
         ey = (/ 0, 1, 0, -1 /)
         phaseqty_phrqc  = 0.d0
         !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(wi,ni,nextI,nextJ,i,j,p)
         Do p = 1,np
           Do j = 1, lx
             Do i = 1, ly
               wi = 0.d0
               If (nodetype(i, j) > 0 .AND. phaseqty(p,i,j) > 0) Then
                 !Get weights
                 Do ni = 1, 4
                   nextI = ey (ni) + i
                   nextJ = ex (ni) + j
                   If (nextI > 0 .And. nextI < ly + 1 .And. nextJ > 0 .And. nextJ &
&                      < lx + 1) Then
                     If (nodetype(nextI, nextJ) <= 0) Then
                       wi = wi + 1.d0
                     End If
                   End If
                 End Do
                 If (wi > 0) Then
                   wi = 1.d0/wi 
                 Else 
                   wi = 0.d0
                 End If
                 !Distribute values
                 Do ni = 1, 4
                   nextI = ey (ni) + i
                   nextJ = ex (ni) + j
                   If (nextI > 0 .And. nextI < ly + 1 .And. nextJ > 0 .And. nextJ &
&                      < lx + 1) Then
                     If (nodetype(nextI, nextJ) <= 0) Then
                       !$OMP ATOMIC
                       phaseqty_phrqc(p,nextI,nextJ) = phaseqty_phrqc(p,nextI,nextJ) + wi * phaseqty(p, i, j)
                     End If
                   End If
                 End Do
               Else If (nodetype(i,j) <=0) Then
                 !$OMP ATOMIC
                 phaseqty_phrqc(p,i,j)=phaseqty_phrqc(p,i,j)+phaseqty(p,i,j)
               End If
             End Do
           End Do
         End Do
         !$OMP END PARALLEL DO
End subroutine
