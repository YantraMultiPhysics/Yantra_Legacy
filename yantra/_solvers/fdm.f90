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
!This file contains relavent subroutines for finite-difference computations used in yantra
!
!=======================================================================================
!
subroutine grad2d(f,df,dx,periodic,ly,lx)
!gradient on uniform grid with spacing equal to dx
!reference for weights at boundary: Fornberg, Bengt. "Generation of finite 
!difference formulas on arbitrarily spaced grids." Mathematics of computation
! 51.184 (1988): 699-706.
    Implicit None
    Real (8), Intent(In):: f(ly,lx),dx
    Real (8), Intent(Out):: df(2,ly,lx)
    Integer, Intent(In):: periodic(2)
    Integer, Intent(In):: ly,lx
    Integer:: i,j 
    !interior point
    !$OMP Parallel default(shared) private(i,j)
    !$OMP Do
    Do i = 1,ly
        Do j = 1,lx
            if ((j > 1) .AND. (j<lx)) then   
                df(1,i,j) = 0.5/dx*(f(i,j+1)-f(i,j-1))
            end if 
            if ((i > 1) .AND. (i<ly)) then
                df(2,i,j) = -1.d0*0.5/dx*(f(i+1,j)-f(i-1,j))
            end if
        End Do
    End Do            
    !$OMP end Do 
    !x-boundary
    if (periodic(1)>0) then
        !$OMP  Do 
        Do i = 1,ly
            df(1,i,1)  =0.5/dx*(f(i,2) - f(i,lx))
            df(1,i,lx) =0.5/dx*(f(i,1) - f(i,lx-1))
        End Do
        !$OMP end Do
    else
        !$OMP  Do 
        Do i = 1,ly
            df(1,i,1)  =0.5/dx*(4*f(i,2)-3*f(i,1)- f(i,3))
            df(1,i,lx) =0.5/dx*(3*f(i,lx)+ f(i,lx-2)-4*f(i,lx-1))
        End do
        !$OMP end Do
    end if
    !y-boundary
    if (periodic(2)>0) then
        !$OMP Do
        Do j =  1,lx
            df(2,1,j)  =-1.d0*0.5/dx*(f(2,j) - f(ly,j))
            df(2,ly,j) =-1.d0*0.5/dx*(f(1,j) - f(ly-1,j))
        end Do
        !$OMP End Do
    else
        !$OMP Do
        Do j = 1,lx
            df(2,1,j)  =-1.d0*0.5/dx*(4*f(2,j)-3*f(1,j)- f(3,j))
            df(2,ly,j) =-1.d0*0.5/dx*(3*f(ly,j)+ f(ly-2,j)-4*f(ly-1,j))
        end Do
        !$OMP end Do
    end if 
    !$OMP end Parallel
end subroutine
!
subroutine div2d(v,div,dx,periodic,ly,lx)
!divergence on uniform grid with spacing equal to dx
!reference for weights at boundary: Fornberg, Bengt. "Generation of finite 
!difference formulas on arbitrarily spaced grids." Mathematics of computation
! 51.184 (1988): 699-706.
    Implicit None
    Real (8), Intent(In):: v(2,ly,lx),dx
    Real (8), Intent(Out):: div(ly,lx)
    Integer, Intent(In):: periodic(2)
    Integer, Intent(In):: ly,lx
    Integer:: i,j 
    Real (8):: dv(2,ly,lx)
    !interior point
    !$OMP Parallel default(shared) private(i,j)
    !$OMP Do
    Do i = 1,ly
        Do j = 1,lx
            if ((j > 1) .AND. (j<lx)) then   
                dv(1,i,j) = 0.5/dx*(v(1,i,j+1)-v(1,i,j-1))
            end if 
            if ((i > 1) .AND. (i<ly)) then
                dv(2,i,j) = -1.d0*0.5/dx*(v(2,i+1,j)-v(2,i-1,j))
            end if
        End Do
    End Do            
    !$OMP end Do 
    !x-boundary
    if (periodic(1)>0) then
        !$OMP  Do 
        Do i = 1,ly
            dv(1,i,1)  =0.5/dx*(v(1,i,2) - v(1,i,lx))
            dv(1,i,lx) =0.5/dx*(v(1,i,1) - v(1,i,lx-1))
        End Do
        !$OMP end Do
    else
        !$OMP  Do 
        Do i = 1,ly
            dv(1,i,1)  =0.5/dx*(4*v(1,i,2)-3*v(1,i,1)- v(1,i,3))
            dv(1,i,lx) =0.5/dx*(3*v(1,i,lx)+ v(1,i,lx-2)-4*v(1,i,lx-1))
        End do
        !$OMP end Do
    end if
    !y-boundary
    if (periodic(2)>0) then
        !$OMP Do
        Do j =  1,lx
            dv(2,1,j)  =-1.d0*0.5/dx*(v(2,2,j) - v(2,ly,j))
            dv(2,ly,j) =-1.d0*0.5/dx*(v(2,1,j) - v(2,ly-1,j))
        end Do
        !$OMP End Do
    else
        !$OMP Do
        Do j = 1,lx
            dv(2,1,j)  =-1.d0*0.5/dx*(4*v(2,2,j)-3*v(2,1,j)- v(2,3,j))
            dv(2,ly,j) =-1.d0*0.5/dx*(3*v(2,ly,j)+ v(2,ly-2,j)-4*v(2,ly-1,j))
        end Do
        !$OMP end Do
    end if 
    !$OMP Do
    Do i = 1,ly
        Do j = 1,lx
            div(i,j) = sum(dv(:,i,j))
        End Do
    End Do            
    !$OMP end Do 
    !$OMP end Parallel
end subroutine
!
subroutine div3d(v,div,dx,periodic,lz,ly,lx)
!reference for weights at boundary: Fornberg, Bengt. "Generation of finite 
!difference formulas on arbitrarily spaced grids." Mathematics of computation
! 51.184 (1988): 699-706.
    Implicit None
    Real (8), Intent(In):: v(3,lz,ly,lx),dx
    Real (8), Intent(Out)::div(lz,ly,lx)
    Integer, Intent(In):: lz,ly,lx
    Integer, Intent(In):: periodic(3)
    Integer:: i,j,k
    Real(8):: dv(3,lz,ly,lx)
    !interior point
    !$OMP Parallel default(shared) private(i,j,k)
    !$OMP Do 
    Do k = 1,lx
        Do j = 1,ly
            Do i = 1,lz
                if ((k > 1) .AND. (k < lx)) then   
                    dv(1,i,j,k) = 0.5/dx*(v(1,i,j,k+1)-v(1,i,j,k-1))
                end if 
                if ((j > 1) .AND. (j < ly)) then
                    dv(2,i,j,k) =-1.d0*0.5/dx*(v(2,i,j+1,k)-v(2,i,j-1,k))
                end if 
                if ((i > 1) .AND. (i < lz)) then
                    dv(3,i,j,k) =-1.d0*0.5/dx*(v(3,i+1,j,k)-v(3,i-1,j,k))
                end if
            End Do
        End Do
    End Do       
    !$OMP end  Do      
    !x-boundary
    if (periodic(1)>0) then
        !$OMP Do
        Do j = 1,ly
            Do i = 1,lz
                dv(1,i,j,1)  =0.5/dx*(v(1,i,j,2) - v(1,i,j,lx))
                dv(1,i,j,lx) =0.5/dx*(v(1,i,j,1) - v(1,i,j,lx-1))
            End Do
        End Do
        !$OMP end Do
    else
       !$OMP Do
        Do j = 1,ly
            Do i = 1,lx
                dv(1,i,j,1)  =0.5/dx*(4*v(1,i,j,2)-3*v(1,i,j,1)- v(1,i,j,3))
                dv(1,i,j,lx) =0.5/dx*(3*v(1,i,j,lx)+ v(1,i,j,lx-2)-4*v(1,i,j,lx-1))
            End Do
        End Do
        !$OMP end Do
    end if
    !y-boundary
    if (periodic(2)>0) then
        !$OMP DO
        Do k  = 1,lx
            Do i = 1,lz
                dv(2,i,1,k)  =-1.d0*0.5/dx*(v(2,i,2,k) - v(2,i,ly,k))
                dv(2,i,ly,k) =-1.d0*0.5/dx*(v(2,i,1,k) - v(2,i,ly-1,k))
            End Do
        End Do
        !$OMP End Do
    else
        !$OMP Do
        Do k = 1,lx
            Do i = 1,lz
                dv(2,i,1,k)  =-1.d0*0.5/dx*(4*v(2,i,2,k)-3*v(2,i,1,k)- v(2,i,3,k))
                dv(2,i,ly,k) =-1.d0*0.5/dx*(3*v(2,i,ly,k)+ v(2,i,ly-2,k)-4*v(2,i,ly-1,k))
            End Do
        End Do
        !$OMP end Do
    end if 
    !z-boundary
    if (periodic(3)>0) then
        !$OMP  do
        Do k = 1,lx
            Do j = 1,ly
                dv(3,1,j,k)  =-1.d0*0.5/dx*(v(3,2,j,k) - v(3,lz,j,k))
                dv(3,lz,j,k) =-1.d0*0.5/dx*(v(3,1,j,k) - v(3,lz-1,j,k))
            End Do
        End Do
        !$OMP End do
    else
        !$OMP Do
        Do k = 1,lx
            Do j = 1,ly
                dv(3,1,j,k)  =-1.d0*0.5/dx*(4*v(3,2,j,k)-3*v(3,1,j,k)- v(3,3,j,k))
                dv(3,lz,j,k) =-1.d0*0.5/dx*(3*v(3,lz,j,k)+ v(3,lz-2,j,k)-4*v(3,lz-1,j,k))
            End Do
        End Do
        !$OMP end do
    end if 
    !$OMP Do 
    Do k = 1,lx
        Do j = 1,ly
            Do i = 1,lz
                div(i,j,k) = sum(dv(:,i,j,k))
            End Do
        End Do
    End Do       
    !$OMP end  Do 
    !$OMP end parallel
end subroutine
!
subroutine grad3d(f,df,dx,periodic,lz,ly,lx)
!reference for weights at boundary: Fornberg, Bengt. "Generation of finite 
!difference formulas on arbitrarily spaced grids." Mathematics of computation
! 51.184 (1988): 699-706.
    Implicit None
    Real (8), Intent(In):: f(lz,ly,lx),dx
    Real (8), Intent(Out):: df(3,lz,ly,lx)
    Integer, Intent(In):: lz,ly,lx
    Integer, Intent(In):: periodic(3)
    Integer:: i,j,k
    !interior point
    !$OMP Parallel default(shared) private(i,j,k)
    !$OMP Do 
    Do k = 1,lx
        Do j = 1,ly
            Do i = 1,lz
                if ((k > 1) .AND. (k < lx)) then   
                    df(1,i,j,k) = 0.5/dx*(f(i,j,k+1)-f(i,j,k-1))
                end if 
                if ((j > 1) .AND. (j < ly)) then
                    df(2,i,j,k) =-1.d0*0.5/dx*(f(i,j+1,k)-f(i,j-1,k))
                end if 
                if ((i > 1) .AND. (i < lz)) then
                    df(3,i,j,k) =-1.d0*0.5/dx*(f(i+1,j,k)-f(i-1,j,k))
                end if
            End Do
        End Do
    End Do       
    !$OMP end  Do      
    !x-boundary
    if (periodic(1)>0) then
        !$OMP Do
        Do j = 1,ly
            Do i = 1,lz
                df(1,i,j,1)  =0.5/dx*(f(i,j,2) - f(i,j,lx))
                df(1,i,j,lx) =0.5/dx*(f(i,j,1) - f(i,j,lx-1))
            End Do
        End Do
        !$OMP end Do
    else
       !$OMP Do
        Do j = 1,ly
            Do i = 1,lx
                df(1,i,j,1)  =0.5/dx*(4*f(i,j,2)-3*f(i,j,1)- f(i,j,3))
                df(1,i,j,lx) =0.5/dx*(3*f(i,j,lx)+ f(i,j,lx-2)-4*f(i,j,lx-1))
            End Do
        End Do
        !$OMP end Do
    end if
    !y-boundary
    if (periodic(2)>0) then
        !$OMP DO
        Do k  = 1,lx
            Do i = 1,lz
                df(2,i,1,k)  =-1.d0*0.5/dx*(f(i,2,k) - f(i,ly,k))
                df(2,i,ly,k) =-1.d0*0.5/dx*(f(i,1,k) - f(i,ly-1,k))
            End Do
        End Do
        !$OMP End Do
    else
        !$OMP Do
        Do k = 1,lx
            Do i = 1,lz
                df(2,i,1,k)  =-1.d0*0.5/dx*(4*f(i,2,k)-3*f(i,1,k)- f(i,3,k))
                df(2,i,ly,k) =-1.d0*0.5/dx*(3*f(i,ly,k)+ f(i,ly-2,k)-4*f(i,ly-1,k))
            End Do
        End Do
        !$OMP end Do
    end if 
    !z-boundary
    if (periodic(3)>0) then
        !$OMP  do
        Do k = 1,lx
            Do j = 1,ly
                df(3,1,j,k)  =-1.d0*0.5/dx*(f(2,j,k) - f(lz,j,k))
                df(3,lz,j,k) =-1.d0*0.5/dx*(f(1,j,k) - f(lz-1,j,k))
            End Do
        End Do
        !$OMP End do
    else
        !$OMP Do
        Do k = 1,lx
            Do j = 1,ly
                df(3,1,j,k)  =-1.d0*0.5/dx*(4*f(2,j,k)-3*f(1,j,k)- f(3,j,k))
                df(3,lz,j,k) =-1.d0*0.5/dx*(3*f(lz,j,k)+ f(lz-2,j,k)-4*f(lz-1,j,k))
            End Do
        End Do
        !$OMP end do
    end if 
    !$OMP end parallel
end subroutine
