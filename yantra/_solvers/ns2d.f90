Subroutine  compute_macro_var(f, rho0, rho, p, u, nodetype,ly, lx)
!computation of macroscopic variables from distribution function 
!for D2Q9 lattice
	Implicit None
	Real(kind = 8), Intent(In) :: f(ly,lx,9)	!Distribution function
	Real(kind = 8), Intent(Inout) :: rho(ly,lx)	!density
	Real(kind = 8), Intent(Inout) :: p(ly,lx)	!pressure
	Real(kind = 8), Intent(Inout) :: u(2,ly,lx)	!velocity
	Real(kind = 8), Intent(In) :: nodetype(ly,lx)	!nodetype
	Real(kind = 8), Intent(In) :: rho0	!density of fluid
	Integer, Intent(In) :: ly,lx	!length in y and x direction respectively
	Integer :: i,j
	Real(kind= 8) :: M1x,M1y
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,M1x,M1y)
	!$OMP DO  
	Do j= 1, lx
	  Do i= 1, ly
	    If (nodetype(i,j) <= 0) Then
	      rho(i, j)= sum(f(i,j,:)) 
	      M1x= f(i,j,2) + f(i,j,4) + f(i,j,5) - f(i,j,6) - f(i,j,8) - f(i,j,9)
	      M1y= f(i,j,3) + f(i,j,4) - f(i,j,5) - f(i,j,7) - f(i,j,8) + f(i,j,9)
	      u(1,i, j)= (1.d0/rho(i,j))*M1x
	      u(2,i, j)= (1.d0/rho(i,j))*M1y
	      p(i, j)=  1.0d0/3.0d0*(rho(i,j) - rho0)
	    Else
	      p(i,j)=0.d0; u(1,i,j)=0.d0; u(2,i,j)=0.d0; rho(i,j)=0.d0
	    End If
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine compute_macro_var
!
Subroutine  get_feq(feq, nodetype, rho, u, ly, lx)
!subroutine to compute the equilibrium distribution function
!for D2Q9 lattice with order 2 terms of u
	Implicit None
	Real(kind = 8), Intent(Out) :: feq(ly,lx,9)	!Distribution function
	Real(kind = 8), Intent(InOut) :: nodetype(ly,lx)	!nodetype
	Real(kind = 8), Intent(In) :: rho(ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(In) :: u(2,ly,lx)	!velocity 
	Integer, Intent(In) :: ly,lx	!length in y and x direction respectively
	Integer :: i, j
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
	!$OMP DO  
	Do j= 1, lx
	  Do i= 1, ly
	   If (nodetype(i,j) <= 0) Then
	    feq(i, j, 1)= rho(i, j)*(-2.0d0/3.0d0*u(1,i, j)**2 - 2.0d0/3.0d0*u(2,i, j)**2 + 4.0d0/9.0d0)
	    feq(i, j, 2)= rho(i, j)*((1.0d0/3.0d0)*u(1,i, j)**2 + (1.0d0/3.0d0)*u(1,i, j) - 1.0d0/6.0d0*u(2,i, j)**2 + &
     	    1.0d0/9.0d0)
	    feq(i, j, 3)= rho(i, j)*(-1.0d0/6.0d0*u(1,i, j)**2 + (1.0d0/3.0d0)*u(2,i, j)**2 + (1.0d0/3.0d0)*u(2,i, j) + &
            1.0d0/9.0d0)
	    feq(i, j, 4)= rho(i, j)*(-1.0d0/24.0d0*u(1,i, j)**2 + (1.0d0/12.0d0)*u(1,i, j) - 1.0d0/24.0d0*u(2,i, j)**2 + ( &
            1.0d0/12.0d0)*u(2,i, j) + (1.0d0/8.0d0)*(u(1,i, j) + u(2,i, j))**2 + 1.0d0/36.0d0)
	    feq(i, j, 5)= rho(i, j)*(-1.0d0/24.0d0*u(1,i, j)**2 + (1.0d0/12.0d0)*u(1,i, j) - 1.0d0/24.0d0*u(2,i, j)**2 - &
            1.0d0/12.0d0*u(2,i, j) + (1.0d0/8.0d0)*(u(1,i, j) - u(2,i, j))**2 + 1.0d0/36.0d0)
	    feq(i, j, 6)= rho(i, j)*((1.0d0/3.0d0)*u(1,i, j)**2 - 1.0d0/3.0d0*u(1,i, j) - 1.0d0/6.0d0*u(2,i, j)**2 + 1.0d0 &
            /9.0d0)
	    feq(i, j, 7)= rho(i, j)*(-1.0d0/6.0d0*u(1,i, j)**2 + (1.0d0/3.0d0)*u(2,i, j)**2 - 1.0d0/3.0d0*u(2,i, j) + &
            1.0d0/9.0d0)
	    feq(i, j, 8)= rho(i, j)*(-1.0d0/24.0d0*u(1,i, j)**2 - 1.0d0/12.0d0*u(1,i, j) - 1.0d0/24.0d0*u(2,i, j)**2 - &
            1.0d0/12.0d0*u(2,i, j) + (1.0d0/8.0d0)*(-u(1,i, j) - u(2,i, j))**2 + 1.0d0/36.0d0)
	    feq(i, j, 9)= rho(i, j)*(-1.0d0/24.0d0*u(1,i, j)**2 - 1.0d0/12.0d0*u(1,i, j) - 1.0d0/24.0d0*u(2,i, j)**2 + ( &
            1.0d0/12.0d0)*u(2,i, j) + (1.0d0/8.0d0)*(-u(1,i, j) + u(2,i, j))**2 + 1.0d0/36.0d0)
          End If
	 End Do
       End Do
       !$OMP END DO
       !$OMP END PARALLEL
End Subroutine get_feq
!
Subroutine  collide_srt(f, rho, u, Fv, nodetype, tau,forcing_model, ly, lx)
!subroutine to compute collision step for single relaxation time scheme
!for D2Q9 lattice with order 2 terms of u considered in feq
	Implicit None
	Real(kind = 8), Intent(Inout) :: f(ly,lx,9)	!Distribution function
	Real(kind = 8), Intent(In) :: rho(ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(InOut) :: u(2,ly,lx)	!velocity in y direction
	Real(kind = 8), Intent(In) :: Fv(2,ly,lx)	!forcing
	Real(kind = 8), Intent(In) :: nodetype(ly,lx) !nodetype
	Real(kind = 8), Intent(In) :: tau(ly,lx)	!tau
  Character(len = *), Intent(In):: forcing_model !forcing model to be used available models 'sc', 'guo', 'mlga'
	Integer, Intent(In) :: ly,lx	!length in y and x direction respectively
	Integer :: i, j, n
	Real(kind = 8):: omega, omega1, feq(9), SS(9), fswap, utemp(2)
  SS=0.0d0
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, n, feq, fswap, SS,omega,omega1, utemp)
	!$OMP DO 
	Do j= 1, lx
	  Do i= 1, ly
            If (nodetype(i,j) <= 0) Then
              omega = 1.0d0/tau(i,j)
              omega1 = 1.0d0 - omega
              utemp = 0.d0
              !modify velocity based for forcing models  before computation of equilibrium distribution function  
              SELECT CASE(forcing_model)
                case ('guo')
                  u(1,i,j) = u(1,i,j) + Fv(1,i,j)/(2.0d0*rho(i,j))
                  u(2,i,j) = u(2,i,j) + Fv(2,i,j)/(2.0d0*rho(i,j))
                case ('sc')
                  utemp= u(:,i,j)
                  u(1,i,j) =u(1,i,j) + (tau(i,j)*Fv(1,i,j))/rho(i,j)
                  u(2,i,j) =u(2,i,j) + (tau(i,j)*Fv(2,i,j))/rho(i,j)
              End select
	      feq(1)= rho(i,j)*(-2.0d0/3.0d0*u(1,i,j)**2 - 2.0d0/3.0d0*u(2,i,j)**2 + 4.0d0/9.0d0)
	      feq(2)= rho(i,j)*((1.0d0/3.0d0)*u(1,i,j)**2 + (1.0d0/3.0d0)*u(1,i,j) - 1.0d0/6.0d0*u(2,i,j)**2 + &
              1.0d0/9.0d0)
	      feq(3)= rho(i,j)*(-1.0d0/6.0d0*u(1,i,j)**2 + (1.0d0/3.0d0)*u(2,i,j)**2 + (1.0d0/3.0d0)*u(2,i,j) + &
              1.0d0/9.0d0)
	      feq(4)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 + (1.0d0/12.0d0)*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 + ( &
              1.0d0/12.0d0)*u(2,i,j) + (1.0d0/8.0d0)*(u(1,i,j) + u(2,i,j))**2 + 1.0d0/36.0d0)
	      feq(5)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 + (1.0d0/12.0d0)*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 - &
              1.0d0/12.0d0*u(2,i,j) + (1.0d0/8.0d0)*(u(1,i,j) - u(2,i,j))**2 + 1.0d0/36.0d0)
	      feq(6)= rho(i,j)*((1.0d0/3.0d0)*u(1,i,j)**2 - 1.0d0/3.0d0*u(1,i,j) - 1.0d0/6.0d0*u(2,i,j)**2 + 1.0d0 &
              /9.0d0)
	      feq(7)= rho(i,j)*(-1.0d0/6.0d0*u(1,i,j)**2 + (1.0d0/3.0d0)*u(2,i,j)**2 - 1.0d0/3.0d0*u(2,i,j) + &
              1.0d0/9.0d0)
	      feq(8)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 - 1.0d0/12.0d0*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 - &
              1.0d0/12.0d0*u(2,i,j) + (1.0d0/8.0d0)*(-u(1,i,j) - u(2,i,j))**2 + 1.0d0/36.0d0)
	      feq(9)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 - 1.0d0/12.0d0*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 + ( &
              1.0d0/12.0d0)*u(2,i,j) + (1.0d0/8.0d0)*(-u(1,i,j) + u(2,i,j))**2 + 1.0d0/36.0d0)
              SELECT CASE(forcing_model)
              case('guo')
               Call get_Guo_forcing(SS,Fv(1,i,j),Fv(2,i,j),u(1,i,j),u(2,i,j),omega)
              case('mlga')
               Call get_mLGA_forcing(SS,Fv(1,i,j),Fv(2,i,j))
               u(1,i,j) = u(1,i,j) + Fv(1,i,j)/(2.0d0*rho(i,j))
               u(2,i,j) = u(2,i,j) + Fv(2,i,j)/(2.0d0*rho(i,j))
              case('sc')
               !only correct the velocity 
               u(1,i,j) = utemp(1) + Fv(1,i,j)/(2.0d0*rho(i,j))
               u(2,i,j) = utemp(2) + Fv(2,i,j)/(2.0d0*rho(i,j))
              End SELECT
              f(i,j,:)= omega1*f(i,j,:) + omega*feq + SS
	      Do n=2,5
	        fswap = f(i,j,n)
	        f(i,j,n)=f(i,j,n+4)
	        f(i,j,n+4)=fswap
	      End Do
	    End If
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
     Contains
       Subroutine get_mLGA_forcing(SS,Fx,Fy)
         Implicit None
         Real(kind = 8), Intent(InOut):: SS(9)
         Real(kind = 8), Intent(In):: Fx, Fy
         SS(1) = 0
         SS(2) = (1.0d0/3.0d0)*Fx
         SS(3) = (1.0d0/3.0d0)*Fy
         SS(4) = (1.0d0/12.0d0)*Fx + (1.0d0/12.0d0)*Fy
         SS(5) = (1.0d0/12.0d0)*Fx - 1.0d0/12.0d0*Fy
         SS(6) = -1.0d0/3.0d0*Fx
         SS(7) = -1.0d0/3.0d0*Fy
         SS(8) = -1.0d0/12.0d0*Fx - 1.0d0/12.0d0*Fy
         SS(9) = -1.0d0/12.0d0*Fx + (1.0d0/12.0d0)*Fy
       End Subroutine
       Subroutine get_Guo_forcing(SS,Fx,Fy,ux,uy,omega)
        Implicit none
        Real(kind = 8), Intent(InOut):: SS(9)
        Real(kind = 8), Intent(In)::Fx,Fy,ux,uy,omega
        Real(kind = 8)::t1
        t1 = -0.5d0*omega + 1.0d0
        SS(1) = t1*(-4.0d0/3.0d0*Fx*ux - 4.0d0/3.0d0*Fy*uy)
        SS(2) = t1*((2.0d0/3.0d0)*Fx*ux + (1.0d0/3.0d0)*Fx - 1.0d0/3.0d0*Fy*uy)
        SS(3) = t1*(-1.0d0/3.0d0*Fx*ux + (2.0d0/3.0d0)*Fy*uy + (1.0d0/3.0d0)*Fy)
        SS(4) = t1*(-1.0d0/12.0d0*Fx*ux + (1.0d0/12.0d0)*Fx - 1.0d0/12.0d0*Fy*uy + ( &
                     1.0d0/12.0d0)*Fy + (1.0d0/4.0d0)*(Fx + Fy)*(ux + uy))
        SS(5) = t1*(-1.0d0/12.0d0*Fx*ux + (1.0d0/12.0d0)*Fx - 1.0d0/12.0d0*Fy*uy - 1.0d0 &
                    /12.0d0*Fy + (1.0d0/4.0d0)*(Fx - Fy)*(ux - uy))
        SS(6) = t1*((2.0d0/3.0d0)*Fx*ux - 1.0d0/3.0d0*Fx - 1.0d0/3.0d0*Fy*uy)
        SS(7) = t1*(-1.0d0/3.0d0*Fx*ux + (2.0d0/3.0d0)*Fy*uy - 1.0d0/3.0d0*Fy)
        SS(8) = t1*(-1.0d0/12.0d0*Fx*ux - 1.0d0/12.0d0*Fx - 1.0d0/12.0d0*Fy*uy - 1.0d0/ &
                     12.0d0*Fy + (1.0d0/4.0d0)*(-Fx - Fy)*(-ux - uy))
        SS(9) = t1*(-1.0d0/12.0d0*Fx*ux - 1.0d0/12.0d0*Fx - 1.0d0/12.0d0*Fy*uy + (1.0d0/ &
                     12.0d0)*Fy + (1.0d0/4.0d0)*(-Fx + Fy)*(-ux + uy))
       End subroutine get_Guo_forcing
End Subroutine collide_srt
!
Subroutine  collide_trt(f, rho, u, Fv, nodetype, tau_s, magic_para,forcing_model, ly, lx)
!subroutine to compute collision step for single relaxation time scheme
!for D2Q9 lattice with order 2 terms of u considered in feq
	Implicit None
	Real(kind = 8), Intent(Inout) :: f(ly,lx,9)	!Distribution function
	Real(kind = 8), Intent(In) :: rho(ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(InOut) :: u(2,ly,lx)	!velocity in y direction
	Real(kind = 8), Intent(In) :: Fv(2,ly,lx)	!forcing
	Real(kind = 8), Intent(In) :: nodetype(ly,lx) !nodetype
  Real(kind = 8), Intent(In) :: tau_s(ly,lx)	!tau for symmetric  moments
  Real(kind = 8), Intent(In):: magic_para !magic_paramter
  Character(len = *), Intent(In):: forcing_model !forcing model to be used available models 'guo', 'mlga'
	Integer, Intent(In) :: ly,lx	!length in y and x direction respectively
	Integer :: i, j, n
	Real(kind = 8):: omega_a, omega_s, tau_a, feq(9),feq_s(9),feq_a(9),f_a(9),f_s(9), SS(9), fswap
  SS=0.0d0
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, n, feq,feq_s,feq_a,f_a,f_s, fswap, SS,omega_a,omega_s, tau_a)
	!$OMP DO 
	Do j= 1, lx
	  Do i= 1, ly
	    If (nodetype(i,j) <= 0) Then
              tau_a = 0.5 + (magic_para/(tau_s(i,j)-0.5)) 
              omega_s = 1.0d0/tau_s(i,j)
              omega_a = 1.0d0/tau_a
              !modify velocity based for forcing models  before computation of equilibrium distribution function  
              SELECT CASE(forcing_model)
                case ('guo')
                  u(1,i,j) = u(1,i,j) + Fv(1,i,j)/(2.0d0*rho(i,j))
                  u(2,i,j) = u(2,i,j) + Fv(2,i,j)/(2.0d0*rho(i,j))
              End select
              feq(1)= rho(i,j)*(-2.0d0/3.0d0*u(1,i,j)**2 - 2.0d0/3.0d0*u(2,i,j)**2 + 4.0d0/9.0d0)
	            feq(2)= rho(i,j)*((1.0d0/3.0d0)*u(1,i,j)**2 + (1.0d0/3.0d0)*u(1,i,j) - 1.0d0/6.0d0*u(2,i,j)**2 + &
              1.0d0/9.0d0)
	            feq(3)= rho(i,j)*(-1.0d0/6.0d0*u(1,i,j)**2 + (1.0d0/3.0d0)*u(2,i,j)**2 + (1.0d0/3.0d0)*u(2,i,j) + &
              1.0d0/9.0d0)
	            feq(4)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 + (1.0d0/12.0d0)*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 + ( &
              1.0d0/12.0d0)*u(2,i,j) + (1.0d0/8.0d0)*(u(1,i,j) + u(2,i,j))**2 + 1.0d0/36.0d0)
	            feq(5)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 + (1.0d0/12.0d0)*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 - &
              1.0d0/12.0d0*u(2,i,j) + (1.0d0/8.0d0)*(u(1,i,j) - u(2,i,j))**2 + 1.0d0/36.0d0)
	            feq(6)= rho(i,j)*((1.0d0/3.0d0)*u(1,i,j)**2 - 1.0d0/3.0d0*u(1,i,j) - 1.0d0/6.0d0*u(2,i,j)**2 + 1.0d0 &
              /9.0d0)
	            feq(7)= rho(i,j)*(-1.0d0/6.0d0*u(1,i,j)**2 + (1.0d0/3.0d0)*u(2,i,j)**2 - 1.0d0/3.0d0*u(2,i,j) + &
              1.0d0/9.0d0)
	            feq(8)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 - 1.0d0/12.0d0*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 - &
              1.0d0/12.0d0*u(2,i,j) + (1.0d0/8.0d0)*(-u(1,i,j) - u(2,i,j))**2 + 1.0d0/36.0d0)
	            feq(9)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 - 1.0d0/12.0d0*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 + ( &
              1.0d0/12.0d0)*u(2,i,j) + (1.0d0/8.0d0)*(-u(1,i,j) + u(2,i,j))**2 + 1.0d0/36.0d0)
!              feq(1) = rho(i,j)-sum(feq(2:9))
              !compute f_s, f_a, feq_s,feq_a
              feq_s = 0.d0
              feq_a = 0.d0
              f_s = 0.d0
              f_a = 0.d0
              feq_s(1)=feq(1)
              f_s(1)=f(i,j,1)
              feq_a(1)=0.d0
              feq_a(1)=0.d0
              Do n=2,5
                feq_s(n)= (feq(n) + feq(n+4))/2.d0
                feq_a(n)= (feq(n) - feq(n+4))/2.d0
                f_s(n) =  (f(i,j,n) + f(i,j,n+4))/2.d0
                f_a(n) =  (f(i,j,n) - f(i,j,n+4))/2.d0
              End Do
              Do n=6,9
                feq_s(n) = feq_s(n-4)
                feq_a(n) = -feq_a(n-4)
                f_s(n) = f_s(n-4)
                f_a(n) = -f_a(n-4)
              End Do
              SELECT CASE(forcing_model)
              CASE ('mlga')
                Call get_mLGA_forcing(SS,Fv(1,i,j),Fv(2,i,j))
                u(1,i,j) = u(1,i,j) + Fv(1,i,j)/(2.0d0*rho(i,j))
                u(2,i,j) = u(2,i,j) + Fv(2,i,j)/(2.0d0*rho(i,j))
              CASE ('guo')
                Call get_Guo_forcing(SS,Fv(1,i,j),Fv(2,i,j),u(1,i,j),u(2,i,j),omega_s,omega_a)
              END SELECT
 	            f(i,j,:) = f(i,j,:) - omega_s*(f_s-feq_s) - omega_a*(f_a-feq_a) + SS
        Do n=2,5
	        fswap = f(i,j,n)
	        f(i,j,n)=f(i,j,n+4)
	        f(i,j,n+4)=fswap
	      End Do
	    End If
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
     Contains
       Subroutine get_mLGA_forcing(SS,Fx,Fy)
         !for details refer to guo et al. 2015 book
         Implicit None
         Real(kind = 8), Intent(InOut):: SS(9)
         Real(kind = 8), Intent(In):: Fx, Fy
         SS(1) = 0
         SS(2) = (1.0d0/3.0d0)*Fx
         SS(3) = (1.0d0/3.0d0)*Fy
         SS(4) = (1.0d0/12.0d0)*Fx + (1.0d0/12.0d0)*Fy
         SS(5) = (1.0d0/12.0d0)*Fx - 1.0d0/12.0d0*Fy
         SS(6) = -1.0d0/3.0d0*Fx
         SS(7) = -1.0d0/3.0d0*Fy
         SS(8) = -1.0d0/12.0d0*Fx - 1.0d0/12.0d0*Fy
         SS(9) = -1.0d0/12.0d0*Fx + (1.0d0/12.0d0)*Fy
       End Subroutine
       Subroutine get_Guo_forcing(SS,Fx,Fy,ux,uy,omega_s,omega_a)
         Implicit None
         Real(kind = 8), Intent(InOut):: SS(9)
         Real(kind = 8), Intent(In):: Fx, Fy, ux, uy, omega_s, omega_a
         Real(kind = 8):: SS_a(9),SS_s(9)
         Integer::n
        SS(1) = (-4.0d0/3.0d0*Fx*ux - 4.0d0/3.0d0*Fy*uy)
        SS(2) = ((2.0d0/3.0d0)*Fx*ux + (1.0d0/3.0d0)*Fx - 1.0d0/3.0d0*Fy*uy)
        SS(3) = (-1.0d0/3.0d0*Fx*ux + (2.0d0/3.0d0)*Fy*uy + (1.0d0/3.0d0)*Fy)
        SS(4) = (-1.0d0/12.0d0*Fx*ux + (1.0d0/12.0d0)*Fx - 1.0d0/12.0d0*Fy*uy + ( &
                     1.0d0/12.0d0)*Fy + (1.0d0/4.0d0)*(Fx + Fy)*(ux + uy))
        SS(5) = (-1.0d0/12.0d0*Fx*ux + (1.0d0/12.0d0)*Fx - 1.0d0/12.0d0*Fy*uy - 1.0d0 &
                    /12.0d0*Fy + (1.0d0/4.0d0)*(Fx - Fy)*(ux - uy))
        SS(6) = ((2.0d0/3.0d0)*Fx*ux - 1.0d0/3.0d0*Fx - 1.0d0/3.0d0*Fy*uy)
        SS(7) = (-1.0d0/3.0d0*Fx*ux + (2.0d0/3.0d0)*Fy*uy - 1.0d0/3.0d0*Fy)
        SS(8) = (-1.0d0/12.0d0*Fx*ux - 1.0d0/12.0d0*Fx - 1.0d0/12.0d0*Fy*uy - 1.0d0/ &
                     12.0d0*Fy + (1.0d0/4.0d0)*(-Fx - Fy)*(-ux - uy))
        SS(9) = (-1.0d0/12.0d0*Fx*ux - 1.0d0/12.0d0*Fx - 1.0d0/12.0d0*Fy*uy + (1.0d0/ &
                     12.0d0)*Fy + (1.0d0/4.0d0)*(-Fx + Fy)*(-ux + uy))
        SS_s = 0.0d0
        SS_a = 0.0d0
        SS_s(1) =  SS(1)
        Do n = 2,5
          SS_s(n) = 0.5*(SS(n) + SS(n+4))
          SS_a(n) = 0.5*(SS(n) - SS(n+4))
        End Do
        Do n = 6,9
          SS_s(n) = SS_s(n-4)
          SS_a(n) = -SS_a(n-4)
        End Do
        SS = (1.0d0-0.5d0*omega_s)*SS_s + (1.0d0-0.5d0*omega_a)*SS_a  
       End Subroutine get_Guo_forcing
End Subroutine collide_trt
!
Subroutine  stream_and_bounce(f, nodetype, periodicity, ly, lx)
	Implicit None
	Real(kind = 8), Intent(Inout) :: f(ly,lx,9)	!Distribution function
	Real(kind = 8), Intent(In) :: nodetype(ly,lx)	!nodetype
     Integer, Intent(In) :: periodicity(2) !flag for periodicity 0 dimension for x, 1 dimension for y 
	Integer, Intent(In) :: ly,lx	!length in y and x direction respectively
	Integer :: ex(9)	!lattice velocity in x direction
	Integer :: ey(9)	!lattice velocity in y direction
	Integer :: i,j, n, half, nextX, nextY
	Real(kind = 8) :: ftemp
	ex = (/0, 1, 0, 1, 1, -1, 0, -1, -1/)
	ey = (/0, 0, 1, 1, -1, 0, -1, -1, 1/)
	half= 4
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, n, nextX, nextY, ftemp)
	!$OMP DO 
	Do j= 1, lx
	  Do i= 1, ly
	    If (nodetype(i,j) <= 0) Then
              Do  n= 2,5
                nextX = j + ex (n)
                nextY = i - ey (n)
	          !Apply periodicity
                If (periodicity(1)>0) Then
                  If (nextX < 1) nextX = lx
                  If (nextX > lx) nextX = 1
                End If 
                If (periodicity(2)>0) Then
                  If (nextY > ly) nextY = 1
                  If (nextY < 1) nextY = ly
                End If
                If (nextX > 0 .And. nextX < lx+1 &
&                   .And. nextY > 0  .And. nextY < ly+1) Then
                  If (nodetype(nextY, nextX) <= 0 .And. &
&                     nodetype(nextY,j)<=0 .OR. nodetype(i,nextX)<=0) Then
                    ftemp = f (nextY, nextX, n)
	            f(nextY, nextX, n) = f (i, j, n+half)
                    f(i, j, n+half) = ftemp
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
Subroutine apply_bc(f,rho0,nodetype,interp,left_type,left_val, &
&     right_type,right_val,top_type,top_val,bottom_type,bottom_val,ly,lx)
    !Subroutine to apply boundary condition for navier stokes equation
    !for d2q9 lattice. Pressure boundary implemented as extrapolated scheme and 
    !velocity boundary implemented as anti-bounceback scheme
    Implicit none
    Real(kind = 8), Intent(InOut):: f(ly,lx,9) ! distribution function
    Real(kind = 8), Intent(In):: rho0 ! density of fluid
    Real(kind = 8), Intent(In):: nodetype(ly,lx) ! nodetype indicator
    Character(Len=*), Intent(In):: left_type ! left boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    Character(Len=*), Intent(In):: right_type ! right boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    Character(Len=*), Intent(In):: top_type ! top boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    character(Len=*), Intent(In):: bottom_type ! bottom boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    Real(kind = 8), Intent(In):: left_val(2) ! value for left boundary
    Real(kind = 8), Intent(In):: right_val(2) ! value for right boundary
    Real(kind = 8), Intent(In):: top_val(2) ! value for top boundary
    Real(kind = 8), Intent(In):: bottom_val(2) ! value for bottom boundary
    Integer, Intent(In):: interp ! 1 to allow interpolation to set boundary condtion 0 to not interpolate at the boundary node
    Integer, Intent(In):: ly ! number of nodes in y direction
    Integer, Intent(In):: lx ! number of nodes in x direction
    !other variables...
    Integer:: i, inext
    Real(kind = 8):: bc_val1,bc_val2,es2,rhonbh,ubx,uby,l,ftemp
    es2= 0.333333333333
    !$OMP PARALLEL DEFAULT (SHARED) PRIVATE(i,inext,bc_val1,bc_val2,rhonbh,ubx,uby,l, ftemp)    
    select case (left_type)
     case ('p')
       !$OMP DO 
        Do i= 1,ly
          If (nodetype(i,1)<=0) Then
            !convert pressure to density
            bc_val1 = left_val(1)/es2 + rho0
            If (interp==1) Then
              rhonbh = sum(f(i,2,:)) 
              bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
            End If
            call pressure_extrapol_bc(f(i,1,:),f(i,2,:),bc_val1)            
          End If
        End Do
        !$OMP END DO
     case ('rho')
       !$OMP DO 
        Do i= 1,ly
          If (nodetype(i,1)<=0) Then
            !convert pressure to density
            bc_val1 = left_val(1) 
            If (interp==1) Then
              rhonbh = sum(f(i,2,:)) 
              bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
            End If
            call pressure_extrapol_bc(f(i,1,:),f(i,2,:),bc_val1)            
          End If
        End Do
        !$OMP END DO
      case ('u')
        !$OMP DO 
        Do i= 1, ly
          If (nodetype(i,1)<=0) Then
            bc_val1 = left_val(1)
            bc_val2 = left_val(2)
            If (interp==1) Then
              call left_vel_bb_bc(f(i,1,:),rho0,bc_val1,bc_val2)
            Else 
              call vel_extrapol_bc(f(i,1,:),f(i,2,:),bc_val1,bc_val2)
            End If
          End If 
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = left_val(1)
        bc_val2 = left_val(2)
        if ((nint(ly*bc_val2)-ly*bc_val2)<0) Then
          l =nint(ly*bc_val2) + 1
        else
          l=nint(ly*bc_val2)
        endif
        !$OMP DO 
        Do i= 1, ly
          If (nodetype(i,1)<=0) Then
            If (interp==1) Then
              ubx= (4.0d0*bc_val1)/(l-2.0d0)**2*(i-1.5d0)*(l-i-0.5d0)
              uby= 0.0d0
              call left_vel_bb_bc(f(i,1,:),rho0,ubx,uby)
            Else
              ubx= (4.0d0*bc_val1)/(l-1.0d0)**2*(i-1.0d0)*(l-i)
              uby= 0.0d0
              call vel_extrapol_bc(f(i,1,:),f(i,2,:),ubx,uby)
            End If         
          End If
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO 
        Do i= 1, ly
          If (nodetype(i,1)<=0) Then
            f(i,1,:) = f(i,2,:)
          End If
        End Do
        !$OMP END DO
    End select
    select case (right_type)
     case ('p')
       !$OMP DO 
        Do i= 1,ly
          If (nodetype(i,lx)<=0) Then
            !convert pressure to density
            bc_val1 = right_val(1)/es2 + rho0
            If (interp==1) Then
              rhonbh = sum(f(i,lx-1,:)) 
              bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
            End If
            call pressure_extrapol_bc(f(i,lx,:),f(i,lx-1,:),bc_val1)            
          End If
        End Do
        !$OMP END DO
     case ('rho')
       !$OMP DO 
        Do i= 1,ly
          If (nodetype(i,lx)<=0) Then
            !convert pressure to density
            bc_val1 = right_val(1)
            If (interp==1) Then
              rhonbh = sum(f(i,lx-1,:)) 
              bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
            End If
            call pressure_extrapol_bc(f(i,lx,:),f(i,lx-1,:),bc_val1)            
          End If
        End Do
        !$OMP END DO
      case ('u')
        !$OMP DO 
        Do i= 1, ly
          If (nodetype(i,lx)<=0) Then
            bc_val1 = right_val(1)
            bc_val2 = right_val(2)
            If (interp==1) Then
              call right_vel_bb_bc(f(i,lx,:),rho0,bc_val1,bc_val2)
            Else 
              call vel_extrapol_bc(f(i,lx,:),f(i,lx-1,:),bc_val1,bc_val2)
            End If
          End If 
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = right_val(1)
        bc_val2 = right_val(2)
        if ((nint(ly*bc_val2)-ly*bc_val2)<0) Then
          l =nint(ly*bc_val2) + 1
        else
          l=nint(ly*bc_val2)
        endif
        !$OMP DO 
        Do i= 1, ly
          If (nodetype(i,lx)<=0) Then
            If (interp==1) Then
              ubx= (4.0d0*bc_val1)/(l-2.0d0)**2*(i-1.5d0)*(l-i-0.5d0)
              uby= 0.0d0
              call right_vel_bb_bc(f(i,lx,:),rho0,ubx,uby)
            Else
              ubx= (4.0d0*bc_val1)/(l-1.0d0)**2*(i-1.0d0)*(l-i)
              uby= 0.0d0
              call vel_extrapol_bc(f(i,lx,:),f(i,lx-1,:),ubx,uby)
            End If         
          End If
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO 
        Do i= 1, ly
          If (nodetype(i,lx)<=0) Then
            f(i,lx,:) = f(i,lx-1,:)
          End If
        End Do
        !$OMP END DO
    End select
    select case (top_type)
     case ('p')
       !$OMP DO 
        Do i= 1,lx
          If (nodetype(1,i)<=0) Then
            !convert pressure to density
            bc_val1 = top_val(1)/es2 + rho0
            If (interp==1) Then
              rhonbh = sum(f(2,i,:)) 
              bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
            End If
            call pressure_extrapol_bc(f(1,i,:),f(2,i,:),bc_val1)            
          End If
        End Do
        !$OMP END DO
     case ('rho')
       !$OMP DO 
        Do i= 1,lx
          If (nodetype(1,i)<=0) Then
            !convert pressure to density
            bc_val1 = top_val(1) 
            If (interp==1) Then
              rhonbh = sum(f(2,i,:)) 
              bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
            End If
            call pressure_extrapol_bc(f(1,i,:),f(2,i,:),bc_val1)            
          End If
        End Do
        !$OMP END DO
      case ('u')
        !$OMP DO 
        Do i= 1, lx
          If (nodetype(1,i)<=0) Then
            bc_val1 = top_val(1)
            bc_val2 = top_val(2)
            If (interp==1) Then
              call top_vel_bb_bc(f(1,i,:),rho0,bc_val1,bc_val2)
            Else 
              call vel_extrapol_bc(f(1,i,:),f(2,i,:),bc_val1,bc_val2)
            End If
          End If 
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = top_val(1)
        bc_val2 = top_val(2)
        if ((nint(lx*bc_val2)-lx*bc_val2)<0) Then
          l =nint(lx*bc_val2) + 1
        else
          l=nint(lx*bc_val2)
        endif
        !$OMP DO 
        Do i= 1, lx
          If (nodetype(1,i)<=0) Then
            If (interp==1) Then
              ubx= 0.0d0
              uby= (4.0d0*bc_val1)/(l-2.0d0)**2*(i-1.5d0)*(l-i-0.5d0)
              call top_vel_bb_bc(f(1,i,:),rho0,ubx,uby)
            Else
              ubx= 0.0d0
              uby= (4.0d0*bc_val1)/(l-1.0d0)**2*(i-1.0d0)*(l-i)
              call vel_extrapol_bc(f(1,i,:),f(2,i,:),ubx,uby)
            End If         
          End If
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO 
        Do i= 1, lx
          If (nodetype(1,i)<=0) Then
            f(1,i,:) = f(2,i,:)
          End If
        End Do
        !$OMP END DO
    End select
    select case (bottom_type)
     case ('p')
       !$OMP DO 
        Do i= 1,lx
          If (nodetype(ly,i)<=0) Then
            !convert pressure to density
            bc_val1 = bottom_val(1)/es2 + rho0
            If (interp==1) Then
              rhonbh = sum(f(ly-1,i,:)) 
              bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
            End If
            call pressure_extrapol_bc(f(ly,i,:),f(ly-1,i,:),bc_val1)            
          End If
        End Do
        !$OMP END DO
     case ('rho')
       !$OMP DO 
        Do i= 1,lx
          If (nodetype(ly,i)<=0) Then
            !convert pressure to density
            bc_val1 = bottom_val(1)
            If (interp==1) Then
              rhonbh = sum(f(ly-1,i,:)) 
              bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
            End If
            call pressure_extrapol_bc(f(ly,i,:),f(ly-1,i,:),bc_val1)            
          End If
        End Do
        !$OMP END DO
      case ('u')
        !$OMP DO 
        Do i= 1, lx
          If (nodetype(ly,i)<=0) Then
            bc_val1 = bottom_val(1)
            bc_val2 = bottom_val(2)
            If (interp==1) Then
              call bottom_vel_bb_bc(f(ly,i,:),rho0,bc_val1,bc_val2)
            Else 
              call vel_extrapol_bc(f(ly,i,:),f(ly-1,i,:),bc_val1,bc_val2)
            End If
          End If 
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = bottom_val(1)
        bc_val2 = bottom_val(2)
        if ((nint(lx*bc_val2)-lx*bc_val2)<0) Then
          l =nint(lx*bc_val2) + 1
        else
          l=nint(lx*bc_val2)
        endif
        !$OMP DO 
        Do i= 1, lx
          If (nodetype(ly,i)<=0) Then
            If (interp==1) Then
              ubx= 0.0d0
              uby= (4.0d0*bc_val1)/(l-2.0d0)**2*(i-1.5d0)*(l-i-0.5d0)
              call bottom_vel_bb_bc(f(ly,i,:),rho0,ubx,uby)
            Else
              ubx= 0.0d0
              uby= (4.0d0*bc_val1)/(l-1.0d0)**2*(i-1.0d0)*(l-i)
              call vel_extrapol_bc(f(ly,i,:),f(ly-1,i,:),ubx,uby)
            End If         
          End If
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO 
        Do i= 1, lx
          If (nodetype(ly,i)<=0) Then
            f(ly,i,:) = f(ly-1,i,:)
          End If
        End Do
        !$OMP END DO
    End select
    !$OMP END PARALLEL
    Contains
    Subroutine pressure_extrapol_bc(f,fnbh,rho)
      Implicit None
      Real(kind = 8), Intent(In):: fnbh(9),rho
      Real(kind = 8), Intent(InOut):: f(9)
      Real(kind = 8) :: feq(9),feqnbh(9),rhonbh,uxnbh,uynbh
      rhonbh = sum(fnbh)
      call get_vel(uxnbh,uynbh,fnbh,rhonbh)
      call get_feq(feq,rho,uxnbh,uynbh)
      call get_feq(feqnbh,rhonbh,uxnbh,uynbh)
      f= feq + (fnbh-feqnbh)
    End Subroutine pressure_extrapol_bc
    Subroutine vel_extrapol_bc(f,fnbh,ux,uy)
      Implicit None
      Real(kind = 8), Intent(In):: fnbh(9),ux,uy
      Real(kind = 8), Intent(InOut):: f(9)
      Real(kind = 8) :: feq(9),feqnbh(9),rhonbh,uxnbh,uynbh
      rhonbh = sum(fnbh)
      call get_vel(uxnbh,uynbh,fnbh,rhonbh)
      call get_feq(feq,rhonbh,ux,uy)
      call get_feq(feqnbh,rhonbh,uxnbh,uynbh)
      f= feq + (fnbh-feqnbh)
    End Subroutine vel_extrapol_bc
    Subroutine left_vel_bb_bc(f,m0,ux,uy)
      Real(kind = 8), Intent(In):: m0,ux,uy
      Real(kind = 8), Intent(InOut):: f(9)
      f(2)= f(2)-(-2.0d0/3.0d0*m0*ux)
      f(4)= f(4)-((1.0d0/6.0d0)*m0*(-ux - uy))
      f(5)= f(5)-((1.0d0/6.0d0)*m0*(-ux + uy))
    End Subroutine left_vel_bb_bc
    Subroutine right_vel_bb_bc(f,m0,ux,uy)
      Real(kind = 8), Intent(In):: m0,ux,uy
      Real(kind = 8), Intent(InOut):: f(9)
      f(6)= f(6)-((2.0d0/3.0d0)*m0*ux)
      f(8)= f(8)-((1.0d0/6.0d0)*m0*(ux + uy))
      f(9)= f(9)-((1.0d0/6.0d0)*m0*(ux - uy))
    End Subroutine right_vel_bb_bc
    Subroutine top_vel_bb_bc(f,m0,ux,uy)
      Real(kind = 8), Intent(In):: m0,ux,uy
      Real(kind = 8), Intent(InOut):: f(9)
      f(7)= f(7)-((2.0d0/3.0d0)*m0*uy)
      f(8)= f(8)-((1.0d0/6.0d0)*m0*(ux + uy))
      f(5)= f(5)-((1.0d0/6.0d0)*m0*(-ux + uy))
    End Subroutine top_vel_bb_bc
    Subroutine bottom_vel_bb_bc(f,m0,ux,uy)
      Real(kind = 8), Intent(In):: m0,ux,uy
      Real(kind = 8), Intent(InOut):: f(9)
      f(9)= f(9)-((1.0d0/6.0d0)*m0*(ux - uy))
      f(3)= f(3)-(-2.0d0/3.0d0*m0*uy)
      f(4)= f(4)-((1.0d0/6.0d0)*m0*(-ux - uy))
    End Subroutine bottom_vel_bb_bc
    Subroutine get_vel(ux,uy,f,rho)
      Implicit None
      Real(kind = 8), Intent(In)::  f(9),rho
      Real(kind = 8), Intent(InOut)::ux,uy
      Real(kind = 8):: m1x,m1y
      m1x= f(2) + f(4) + f(5) - f(6) - f(8) - f(9)
      m1y= f(3) + f(4) - f(5) - f(7) - f(8) + f(9)
      ux= m1x/rho
      uy= m1y/rho
    End Subroutine get_vel
    Subroutine get_feq(feq,rho,ux,uy)
      Implicit None
      Real(kind = 8), Intent(In)::  rho,ux,uy
      Real(kind = 8), Intent(InOut)::feq(9)
      feq(1)=rho*(-0.666666666666667d0*ux**2 - 0.666666666666667d0*uy**2 + &
      0.444444444444444d0)
      feq(2)=rho*(0.333333333333333d0*ux**2 + 0.333333333333333d0*ux - &
      0.166666666666667d0*uy**2 + 0.111111111111111d0)
      feq(3)=rho*(-0.166666666666667d0*ux**2 + 0.333333333333333d0*uy**2 + &
      0.333333333333333d0*uy + 0.111111111111111d0)
      feq(4)=rho*(-0.0416666666666667d0*ux**2 + 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 + 0.0833333333333333d0*uy + 0.125d0*( &
      ux + uy)**2 + 0.0277777777777778d0)
      feq(5)=rho*(-0.0416666666666667d0*ux**2 + 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 - 0.0833333333333333d0*uy + 0.125d0*( &
      ux - uy)**2 + 0.0277777777777778d0)
      feq(6)=rho*(0.333333333333333d0*ux**2 - 0.333333333333333d0*ux - &
      0.166666666666667d0*uy**2 + 0.111111111111111d0)
      feq(7)=rho*(-0.166666666666667d0*ux**2 + 0.333333333333333d0*uy**2 - &
      0.333333333333333d0*uy + 0.111111111111111d0)
      feq(8)=rho*(-0.0416666666666667d0*ux**2 - 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 - 0.0833333333333333d0*uy + 0.125d0*( &
      -ux - uy)**2 + 0.0277777777777778d0)
      feq(9)=rho*(-0.0416666666666667d0*ux**2 - 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 + 0.0833333333333333d0*uy + 0.125d0*( &
      -ux + uy)**2 + 0.0277777777777778d0)
    End Subroutine get_feq
End Subroutine
































