Subroutine  compute_macro_var(f, rho0, rho, poros, p, u, nodetype, lz, ly, lx)
!computation of macroscopic variables from distribution function 
!for D3Q19 lattice
	Implicit None
	Real(kind = 8), Intent(In) :: f(lz,ly,lx,19)	!Distribution function
	Real(kind = 8), Intent(Inout) :: rho(lz,ly,lx)	!density
	Real(kind = 8), Intent(Inout) :: poros(lz,ly,lx)	!porosity
	Real(kind = 8), Intent(Inout) :: p(lz,ly,lx)	!pressure
	Real(kind = 8), Intent(Inout) :: u(3,lz,ly,lx)	!velocity
	Real(kind = 8), Intent(In) :: nodetype(lz,ly,lx)	!nodetype
	Real(kind = 8), Intent(In) :: rho0	!density of fluid
	Integer, Intent(In) :: lz, ly, lx	!length in z, y and x direction respectively
	Integer :: i, j, k
	Real(kind = 8) :: M1x, M1y, M1z
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,M1x,M1y,M1z)
	!$OMP DO 
	Do k= 1, lx
	  Do j= 1, ly
	    Do i= 1, lz
	      If (nodetype(i,j,k) <= 0) Then
	        M1x= 0.d0; M1y= 0.d0; M1z= 0.d0
	        rho(i, j, k)= sum(f(i,j,k,:))/poros(i,j,k)
                M1x= f(i,j,k,2) + f(i,j,k,5) + f(i,j,k,6) +f(i,j,k,7) + f(i,j,k,8)-f(i,j,k,11) &
&                -f(i,j,k,14) - f(i,j,k,15) - f(i,j,k,16) - f(i,j,k,17)  
                M1y= f(i,j,k,3) + f(i,j,k,7) - f(i,j,k,8) + f(i,j,k,9) - f(i,j,k,10) - &
&                 f(i,j,k,12) - f(i,j,k,16) + f(i,j,k,17) - f(i,j,k,18) + f(i,j,k,19)   
                M1z= f(i,j,k,4) + f(i,j,k,5) - f(i,j,k,6) + f(i,j,k,9) + f(i,j,k,10) - &
&                 f(i,j,k,13) - f(i,j,k,14) + f(i,j,k,15) - f(i,j,k,18) - f(i,j,k,19)   
	        u(1,i, j, k)= (1.d0/(rho(i,j,k)*poros(i,j,k)))*M1x
	        u(2,i, j, k)= (1.d0/(rho(i,j,k)*poros(i,j,k)))*M1y
	        u(3,i, j, k)= (1.d0/(rho(i,j,k)*poros(i,j,k)))*M1z
	        p(i, j, k)= 1/3d0*(rho(i,j,k) - rho0)
	      Else
	        rho(i,j,k)=0.d0; p(i, j, k)=0.d0; u(1,i, j, k)=0.d0; u(2,i, j, k)=0.d0; u(3,i, j, k)=0.d0
	      End If
	    End Do
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine compute_macro_var
!
Subroutine  get_feq(feq, nodetype, rho,poros, u, lz, ly, lx)
!subroutine to compute the equilibrium distribution function
!for D3Q19 lattice with order 2 terms of u
	Implicit None
	Real(kind = 8), Intent(out) :: feq(lz,ly,lx,19) !Distribution function
	Real(kind = 8), Intent(In) :: rho(lz,ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(In) :: nodetype(lz,ly,lx)!nodetype
	Real(kind = 8), Intent(In):: poros(lz,ly,lx)!porosity	
	Real(kind = 8), Intent(In) :: u(3,lz,ly,lx)	!velocity 
	Integer, Intent(In) :: lz,ly,lx	!length in z,y, and x direction respectively
	Integer :: i, j, k
	Real(kind=8)::pRho
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j, k, pRho)
	!$OMP DO  
	Do k= 1, lx
	  Do j= 1, ly
	    Do i= 1, lz
           If (nodetype(i,j,k) <= 0) Then
           	pRho = rho(i,j,k) * poros(i,j,k)
	      feq(i, j, k, 1)= pRho*(-1.0/2.0*u(1,i, j, k)**2 - 1.0/2.0*u(2,i, j, k)**2 - 1.0/2.0*u(3,i, j, k)**2 + &
      1.0/3.0)
	      feq(i, j, k, 2)= pRho*((1.0/6.0)*u(1,i, j, k)**2 + (1.0/6.0)*u(1,i, j, k) - 1.0/12.0*u(2,i, j, k)**2 - &
      1.0/12.0*u(3,i, j, k)**2 + 1.0/18.0)
	      feq(i, j, k, 3)= pRho*(-1.0/12.0*u(1,i, j, k)**2 + (1.0/6.0)*u(2,i, j, k)**2 + (1.0/6.0)*u(2,i, j, k) - &
      1.0/12.0*u(3,i, j, k)**2 + 1.0/18.0)
	      feq(i, j, k, 4)= pRho*(-1.0/12.0*u(1,i, j, k)**2 - 1.0/12.0*u(2,i, j, k)**2 + (1.0/6.0)*u(3,i, j, k)**2 + &
      (1.0/6.0)*u(3,i, j, k) + 1.0/18.0)
	      feq(i, j, k, 5)= pRho*(-1.0/24.0*u(1,i, j, k)**2 + (1.0/12.0)*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/24.0*u(3,i, j, k)**2 + (1.0/12.0)*u(3,i, j, k) + (1.0/8.0)*(u(1,i, j, k) + &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 6)= pRho*(-1.0/24.0*u(1,i, j, k)**2 + (1.0/12.0)*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/24.0*u(3,i, j, k)**2 - 1.0/12.0*u(3,i, j, k) + (1.0/8.0)*(u(1,i, j, k) - u(3,i, j, k) &
      )**2 + 1.0/36.0)
	      feq(i, j, k, 7)= pRho*(-1.0/24.0*u(1,i, j, k)**2 + (1.0/12.0)*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 + ( &
      1.0/12.0)*u(2,i, j, k) - 1.0/24.0*u(3,i, j, k)**2 + (1.0/8.0)*(u(1,i, j, k) + &
      u(2,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 8)= pRho*(-1.0/24.0*u(1,i, j, k)**2 + (1.0/12.0)*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/12.0*u(2,i, j, k) - 1.0/24.0*u(3,i, j, k)**2 + (1.0/8.0)*(u(1,i, j, k) - u(2,i, j, k) &
      )**2 + 1.0/36.0)
	      feq(i, j, k, 9)= pRho*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/24.0*u(2,i, j, k)**2 + (1.0/12.0)*u(2,i, j, k) - &
      1.0/24.0*u(3,i, j, k)**2 + (1.0/12.0)*u(3,i, j, k) + (1.0/8.0)*(u(2,i, j, k) + &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 10)= pRho*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/24.0*u(2,i, j, k)**2 - 1.0/12.0*u(2,i, j, k) - &
      1.0/24.0*u(3,i, j, k)**2 + (1.0/12.0)*u(3,i, j, k) + (1.0/8.0)*(-u(2,i, j, k) + &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 11)= pRho*((1.0/6.0)*u(1,i, j, k)**2 - 1.0/6.0*u(1,i, j, k) - 1.0/12.0*u(2,i, j, k)**2 - &
      1.0/12.0*u(3,i, j, k)**2 + 1.0/18.0)
	      feq(i, j, k, 12)= pRho*(-1.0/12.0*u(1,i, j, k)**2 + (1.0/6.0)*u(2,i, j, k)**2 - 1.0/6.0*u(2,i, j, k) - &
      1.0/12.0*u(3,i, j, k)**2 + 1.0/18.0)
	      feq(i, j, k, 13)= pRho*(-1.0/12.0*u(1,i, j, k)**2 - 1.0/12.0*u(2,i, j, k)**2 + (1.0/6.0)*u(3,i, j, k)**2 - &
      1.0/6.0*u(3,i, j, k) + 1.0/18.0)
	      feq(i, j, k, 14)= pRho*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/12.0*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/24.0*u(3,i, j, k)**2 - 1.0/12.0*u(3,i, j, k) + (1.0/8.0)*(-u(1,i, j, k) - &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 15)= pRho*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/12.0*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/24.0*u(3,i, j, k)**2 + (1.0/12.0)*u(3,i, j, k) + (1.0/8.0)*(-u(1,i, j, k) + &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 16)= pRho*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/12.0*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/12.0*u(2,i, j, k) - 1.0/24.0*u(3,i, j, k)**2 + (1.0/8.0)*(-u(1,i, j, k) - &
      u(2,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 17)= pRho*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/12.0*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 + ( &
      1.0/12.0)*u(2,i, j, k) - 1.0/24.0*u(3,i, j, k)**2 + (1.0/8.0)*(-u(1,i, j, k) + &
      u(2,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 18)= pRho*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/24.0*u(2,i, j, k)**2 - 1.0/12.0*u(2,i, j, k) - &
      1.0/24.0*u(3,i, j, k)**2 - 1.0/12.0*u(3,i, j, k) + (1.0/8.0)*(-u(2,i, j, k) - &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 19)= pRho*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/24.0*u(2,i, j, k)**2 + (1.0/12.0)*u(2,i, j, k) - &
      1.0/24.0*u(3,i, j, k)**2 - 1.0/12.0*u(3,i, j, k) + (1.0/8.0)*(u(2,i, j, k) - u(3,i, j, k) &
      )**2 + 1.0/36.0)
          End If
	    End Do
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine get_feq
!
Subroutine  collide_srt(f, rho, poros, u, Fv, nodetype, tau, forcing_model,periodicity, lz, ly, lx)
!subroutine to compute collision step for single relaxation time scheme
!for D3Q19 lattice with order 2 terms of u considered in feq
	Implicit None
	Real(kind = 8), Intent(Inout) :: f(lz, ly,lx,19)	!Distribution function
	Real(kind = 8), Intent(In) :: rho(lz,ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(In) :: poros(lz,ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(InOut) :: u(3,lz,ly,lx)	!velocity in z direction
    Real(kind = 8), Intent(In) :: Fv(3,lz,ly,lx) !Volumetric Forcing 
	Real(kind = 8), Intent(In) :: nodetype(lz,ly,lx)	!nodetype
	Integer, Intent(In):: periodicity(3) !periodicity deifinition 
	Real(kind = 8), Intent(In) :: tau(lz,ly,lx)	!nodetype
	character(len = *), Intent(In):: forcing_model !forcing model for LB. Available models 'sc', 'guo', 'mlga'
	Integer, Intent(In) :: lz,ly,lx	!length in z, y and x direction respectively
	Integer :: i, j, k, n
	Real(kind = 8):: omega, omega1, pRho, feq(19), fswap, ss(19), utemp(3),Fpcorr(19),Fvlocal(3)
    ss = 0.d0
    Fpcorr = 0.d0
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, n, feq, fswap, ss, omega, omega1, utemp, pRho,Fpcorr,Fvlocal)
	!$OMP DO 
	Do k= 1, lx
	   Do j= 1, ly
	    Do i= 1, lz
	      If (nodetype(i,j,k) <= 0) Then
	        omega = 1.d0/tau(i,j,k)
        	omega1 = 1.d0 - omega
        	pRho = poros(i,j,k)*rho(i,j,k)
        	!compute pressure corrrection force 
        	!modify velocity to calculate equilibrium distribution function for some forcing models  
        	Fvlocal = Fv(:,i,j,k)
        	SELECT CASE(forcing_model)
        	Case ('guo')
        		u(1,i,j,k) = u(1,i,j,k) + Fvlocal(1)/(2.0*pRho)
             	u(2,i,j,k) = u(2,i,j,k) + Fvlocal(2)/(2.0*pRho)
             	u(3,i,j,k) = u(3,i,j,k) + Fvlocal(3)/(2.0*pRho)
            Case ('sc')
            	utemp= u(:,i,j,k)
                u(1,i,j,k) =u(1,i,j,k) + (tau(i,j,k)*Fvlocal(1))/pRho
                u(2,i,j,k) =u(2,i,j,k) + (tau(i,j,k)*Fvlocal(2))/pRho
                u(3,i,j,k) =u(3,i,j,k) + (tau(i,j,k)*Fvlocal(3))/pRho
            END SELECT
  		    feq(1)=  pRho*(-1.0/2.0*u(1,i,j,k)**2 - 1.0/2.0*u(2,i,j,k)**2 - 1.0/2.0*u(3,i,j,k)**2 + &
      1.0/3.0)
	        feq(2)=  pRho*((1.0/6.0)*u(1,i,j,k)**2 + (1.0/6.0)*u(1,i,j,k) - 1.0/12.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(3)=  pRho*(-1.0/12.0*u(1,i,j,k)**2 + (1.0/6.0)*u(2,i,j,k)**2 + (1.0/6.0)*u(2,i,j,k) - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(4)=  pRho*(-1.0/12.0*u(1,i,j,k)**2 - 1.0/12.0*u(2,i,j,k)**2 + (1.0/6.0)*u(3,i,j,k)**2 + &
      (1.0/6.0)*u(3,i,j,k) + 1.0/18.0)
	        feq(5)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(u(1,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(6)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(u(1,i,j,k) - u(3,i,j,k) &
      )**2 + 1.0/36.0)
	        feq(7)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 + ( &
      1.0/12.0)*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(u(1,i,j,k) + &
      u(2,i,j,k))**2 + 1.0/36.0)
	        feq(8)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(u(1,i,j,k) - u(2,i,j,k) &
      )**2 + 1.0/36.0)
	        feq(9)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 + (1.0/12.0)*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(u(2,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(10)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 - 1.0/12.0*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(-u(2,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(11)=  pRho*((1.0/6.0)*u(1,i,j,k)**2 - 1.0/6.0*u(1,i,j,k) - 1.0/12.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(12)=  pRho*(-1.0/12.0*u(1,i,j,k)**2 + (1.0/6.0)*u(2,i,j,k)**2 - 1.0/6.0*u(2,i,j,k) - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(13)=  pRho*(-1.0/12.0*u(1,i,j,k)**2 - 1.0/12.0*u(2,i,j,k)**2 + (1.0/6.0)*u(3,i,j,k)**2 - &
      1.0/6.0*u(3,i,j,k) + 1.0/18.0)
	        feq(14)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(-u(1,i,j,k) - &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(15)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(-u(1,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(16)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(-u(1,i,j,k) - &
      u(2,i,j,k))**2 + 1.0/36.0)
	        feq(17)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 + ( &
      1.0/12.0)*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(-u(1,i,j,k) + &
      u(2,i,j,k))**2 + 1.0/36.0)
	        feq(18)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 - 1.0/12.0*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(-u(2,i,j,k) - &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(19)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 + (1.0/12.0)*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(u(2,i,j,k) - u(3,i,j,k) &
      )**2 + 1.0/36.0)
	        SELECT CASE(forcing_model)
	        	CASE('guo')
	            	Call get_Guo_forcing(SS,Fvlocal(1),Fvlocal(2),Fvlocal(3),u(1,i,j,k),u(2,i,j,k),u(3,i,j,k),omega)
              	CASE('mlga')
                	Call get_mLGA_forcing(SS,Fvlocal(1),Fvlocal(2),Fvlocal(3))
                	u(1,i,j,k) = u(1,i,j,k) + Fvlocal(1)/(2.0d0*pRho)
                	u(2,i,j,k) = u(2,i,j,k) + Fvlocal(2)/(2.0d0*pRho)
                	u(3,i,j,k) = u(3,i,j,k) + Fvlocal(3)/(2.0d0*pRho)
              	CASE('sc')
                	!only correct the velocity 
                	u(1,i,j,k) = utemp(1) + Fvlocal(1)/(2.0d0*pRho)
                	u(2,i,j,k) = utemp(2) + Fvlocal(2)/(2.0d0*pRho)
                	u(3,i,j,k) = utemp(3) + Fvlocal(3)/(2.0d0*pRho)
	        End SELECT
             f(i,j,k,:)= omega1*f(i,j,k,:) + omega*feq + ss
              call get_Fpcorr(Fpcorr,i,j,k,rho(i,j,k),poros, periodicity,lz,ly,lx)
             f(i,j,k,:)=f(i,j,k,:)+Fpcorr
	        Do n=2,10
	          fswap = f(i,j,k,n)
	          f(i,j,k,n)=f(i,j,k,n+9)
	          f(i,j,k,n+9)=fswap
	        End Do
	      End If
	    End Do
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
     Contains
     	Subroutine get_Fpcorr(Fpcorr,i,j,k,rho,poros, periodicity,lz,ly,lx)
     		Implicit None
     		Real(kind = 8), Intent(Inout)::Fpcorr(19)
        Real(kind = 8), Intent(In):: rho
     		Real(kind = 8), Intent(In):: poros(lz,ly,lx)
        Integer,Intent(In)::i,j,k,lz,ly,lx,periodicity(3)
     		Integer,parameter, dimension(19)::	ex = (/0, 1, 0, 0, 1, 1, 1, 1, 0, 0, -1, 0, 0, -1, -1, -1, -1, 0, 0/)
			  Integer,parameter, dimension(19)::	ey = (/0, 0, 1, 0, 0, 0, 1, -1, 1, -1, 0, -1, 0, 0, 0, -1, 1, -1, 1/)
			  Integer,parameter, dimension(19)::	ez = (/0, 0, 0, 1, 1, -1, 0, 0, 1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1/)
			  Real(8),parameter, dimension(19)::	w = (/0.333333333333333d0, 0.0555555555555556d0, 0.0555555555555556d0, &
			                                            0.0555555555555556d0, 0.0277777777777778d0, 0.0277777777777778d0, &
			                                             0.0277777777777778d0, 0.0277777777777778d0, 0.0277777777777778d0, &
			                                            0.0277777777777778d0, 0.0555555555555556d0, 0.0555555555555556d0, &
			                                            0.0555555555555556d0, 0.0277777777777778d0, 0.0277777777777778d0, &
			                                            0.0277777777777778d0, 0.0277777777777778d0, 0.0277777777777778d0, &
			                                            0.0277777777777778d0/)
			 Integer::n, nextX,nextY,nextZ
  			Fpcorr= 0.d0
     		Do n  = 2, 19
     			nextZ = i-ez(n)
     			nextY = j-ey(n)
     			nextX = k+ex(n)
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
                If (nextX > 0 .And. nextX < lx+1 &
 	          & .And. nextY > 0  .And. nextY < ly+1 &
	          & .And. nextZ > 0  .And. nextZ < lz+1) Then
     				Fpcorr(n)=  w(n)* (poros(nextZ,nextY,nextX)-poros(i,j,k)) *rho
     			End If
     		End Do
     	End Subroutine
        Subroutine get_mLGA_forcing(SS,Fx,Fy,Fz)
          Implicit None
          Real(kind = 8), Intent(InOut):: SS(19)
          Real(kind = 8), Intent(In):: Fx, Fy, Fz
          SS(1) = 0
          SS(2) = (1.0/6.0)*Fx
          SS(3) = (1.0/6.0)*Fy
          SS(4) = (1.0/6.0)*Fz
          SS(5) = (1.0/12.0)*Fx + (1.0/12.0)*Fz
          SS(6) = (1.0/12.0)*Fx - 1.0/12.0*Fz
          SS(7) = (1.0/12.0)*Fx + (1.0/12.0)*Fy
          SS(8) = (1.0/12.0)*Fx - 1.0/12.0*Fy
          SS(9) = (1.0/12.0)*Fy + (1.0/12.0)*Fz
          SS(10) = -1.0/12.0*Fy + (1.0/12.0)*Fz
          SS(11) = -1.0/6.0*Fx
          SS(12) = -1.0/6.0*Fy
          SS(13) = -1.0/6.0*Fz
          SS(14) = -1.0/12.0*Fx - 1.0/12.0*Fz
          SS(15) = -1.0/12.0*Fx + (1.0/12.0)*Fz
          SS(16) = -1.0/12.0*Fx - 1.0/12.0*Fy
          SS(17) = -1.0/12.0*Fx + (1.0/12.0)*Fy
          SS(18) = -1.0/12.0*Fy - 1.0/12.0*Fz
          SS(19) = (1.0/12.0)*Fy - 1.0/12.0*Fz
        End Subroutine
        Subroutine get_Guo_forcing(SS,Fx,Fy,Fz,ux,uy,uz,omega)
        	Real(8), Intent(InOut):: SS(19)
        	Real(8), Intent(In)::Fx, Fy, Fz, ux, uy, uz, omega
        	Real(8):: t1
        	t1 = 1.0d0 - 0.5d0*omega
 			SS(1) = t1*(-1.0d0*Fx*ux - 1.0d0*Fy*uy - 1.0d0*Fz*uz)
 			SS(2) = t1*(0.333333333333333d0*Fx*ux + 0.166666666666667d0*Fx - &
      		0.166666666666667d0*Fy*uy - 0.166666666666667d0*Fz*uz)
			 SS(3) = t1*(-0.166666666666667d0*Fx*ux + 0.333333333333333d0*Fy*uy + &
			      0.166666666666667d0*Fy - 0.166666666666667d0*Fz*uz)
			 SS(4) = t1*(-0.166666666666667d0*Fx*ux - 0.166666666666667d0*Fy*uy + &
			      0.333333333333333d0*Fz*uz + 0.166666666666667d0*Fz)
			 SS(5) = t1*(-0.0833333333333333d0*Fx*ux + 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fz*uz + &
			      0.0833333333333333d0*Fz + 0.25d0*(Fx + Fz)*(ux + uz))
			 SS(6) = t1*(-0.0833333333333333d0*Fx*ux + 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fz*uz - &
			      0.0833333333333333d0*Fz + 0.25d0*(Fx - Fz)*(ux - uz))
			 SS(7) = t1*(-0.0833333333333333d0*Fx*ux + 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy + 0.0833333333333333d0*Fy - &
			      0.0833333333333333d0*Fz*uz + 0.25d0*(Fx + Fy)*(ux + uy))
			 SS(8) = t1*(-0.0833333333333333d0*Fx*ux + 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fy - &
			      0.0833333333333333d0*Fz*uz + 0.25d0*(Fx - Fy)*(ux - uy))
			 SS(9) = t1*(-0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fy*uy + &
			      0.0833333333333333d0*Fy - 0.0833333333333333d0*Fz*uz + &
			      0.0833333333333333d0*Fz + 0.25d0*(Fy + Fz)*(uy + uz))
			 SS(10) = t1*(-0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fy*uy - &
			      0.0833333333333333d0*Fy - 0.0833333333333333d0*Fz*uz + &
			      0.0833333333333333d0*Fz + 0.25d0*(-Fy + Fz)*(-uy + uz))
			 SS(11) = t1*(0.333333333333333d0*Fx*ux - 0.166666666666667d0*Fx - &
			      0.166666666666667d0*Fy*uy - 0.166666666666667d0*Fz*uz)
			 SS(12) = t1*(-0.166666666666667d0*Fx*ux + 0.333333333333333d0*Fy*uy - &
			      0.166666666666667d0*Fy - 0.166666666666667d0*Fz*uz)
			 SS(13) = t1*(-0.166666666666667d0*Fx*ux - 0.166666666666667d0*Fy*uy + &
			      0.333333333333333d0*Fz*uz - 0.166666666666667d0*Fz)
			 SS(14) = t1*(-0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fz*uz - &
			      0.0833333333333333d0*Fz + 0.25d0*(-Fx - Fz)*(-ux - uz))
			 SS(15) = t1*(-0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fz*uz + &
			      0.0833333333333333d0*Fz + 0.25d0*(-Fx + Fz)*(-ux + uz))
			 SS(16) = t1*(-0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fy - &
			      0.0833333333333333d0*Fz*uz + 0.25d0*(-Fx - Fy)*(-ux - uy))
			 SS(17) = t1*(-0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy + 0.0833333333333333d0*Fy - &
			      0.0833333333333333d0*Fz*uz + 0.25d0*(-Fx + Fy)*(-ux + uy))
			 SS(18) = t1*(-0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fy*uy - &
			      0.0833333333333333d0*Fy - 0.0833333333333333d0*Fz*uz - &
			      0.0833333333333333d0*Fz + 0.25d0*(-Fy - Fz)*(-uy - uz))
			 SS(19) = t1*(-0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fy*uy + &
			      0.0833333333333333d0*Fy - 0.0833333333333333d0*Fz*uz - &
			      0.0833333333333333d0*Fz + 0.25d0*(Fy - Fz)*(uy - uz))
        End Subroutine
End Subroutine collide_srt
!
!
Subroutine  collide_trt(f, rho, poros, u, Fv, nodetype, tau_s, magic_para, forcing_model,periodicity, lz, ly, lx)
!subroutine to compute collision step for single relaxation time scheme
!for D3Q19 lattice with order 2 terms of u considered in feq
	Implicit None
	Real(kind = 8), Intent(Inout) :: f(lz, ly,lx,19)	!Distribution function
	Real(kind = 8), Intent(In) :: rho(lz,ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(In) :: poros(lz,ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(InOut) :: u(3,lz,ly,lx)	!velocity in z direction
    Real(kind = 8), Intent(In) :: Fv(3,lz,ly,lx) !Volumetric Forcing 
	Real(kind = 8), Intent(In) :: nodetype(lz,ly,lx)	!nodetype
	Real(kind = 8), Intent(In) :: magic_para !magic parameter
	Integer, Intent(In):: periodicity(3) !periodicity deifinition 
	Real(kind = 8), Intent(In) :: tau_s(lz,ly,lx)	!nodetype
	character(len = *), Intent(In):: forcing_model !forcing model for LB. Available models 'sc', 'guo', 'mlga'
	Integer, Intent(In) :: lz,ly,lx	!length in z, y and x direction respectively
	Integer :: i, j, k, n
	Real(kind = 8):: omega_s, omega_a, tau_a, pRho, feq(19), feq_a(19), feq_s(19),f_a(19), f_s(19), &
	& fswap, ss(19), utemp(3),Fpcorr(19),Fvlocal(3)
    ss = 0.d0
    Fpcorr = 0.d0
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, n, feq, feq_s, feq_a ,f_s ,f_a,  fswap, ss, omega_s, &
	!$OMP& omega_a, tau_a, utemp, pRho,Fpcorr,Fvlocal)
	!$OMP DO 
	Do k= 1, lx
	   Do j= 1, ly
	    Do i= 1, lz
	      If (nodetype(i,j,k) <= 0) Then
            tau_a = 0.5 + (magic_para/(tau_s(i,j,k)-0.5)) 
            omega_s = 1.0d0/tau_s(i,j,k)
            omega_a = 1.0d0/tau_a
        	pRho = poros(i,j,k)*rho(i,j,k)
        	!compute pressure corrrection force 
        	!modify velocity to calculate equilibrium distribution function for some forcing models  
        	Fvlocal = Fv(:,i,j,k)
        	SELECT CASE(forcing_model)
        	Case ('guo')
        		u(1,i,j,k) = u(1,i,j,k) + Fvlocal(1)/(2.0*pRho)
             	u(2,i,j,k) = u(2,i,j,k) + Fvlocal(2)/(2.0*pRho)
             	u(3,i,j,k) = u(3,i,j,k) + Fvlocal(3)/(2.0*pRho)
            END SELECT
  		    feq(1)=  pRho*(-1.0/2.0*u(1,i,j,k)**2 - 1.0/2.0*u(2,i,j,k)**2 - 1.0/2.0*u(3,i,j,k)**2 + &
      1.0/3.0)
	        feq(2)=  pRho*((1.0/6.0)*u(1,i,j,k)**2 + (1.0/6.0)*u(1,i,j,k) - 1.0/12.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(3)=  pRho*(-1.0/12.0*u(1,i,j,k)**2 + (1.0/6.0)*u(2,i,j,k)**2 + (1.0/6.0)*u(2,i,j,k) - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(4)=  pRho*(-1.0/12.0*u(1,i,j,k)**2 - 1.0/12.0*u(2,i,j,k)**2 + (1.0/6.0)*u(3,i,j,k)**2 + &
      (1.0/6.0)*u(3,i,j,k) + 1.0/18.0)
	        feq(5)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(u(1,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(6)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(u(1,i,j,k) - u(3,i,j,k) &
      )**2 + 1.0/36.0)
	        feq(7)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 + ( &
      1.0/12.0)*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(u(1,i,j,k) + &
      u(2,i,j,k))**2 + 1.0/36.0)
	        feq(8)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(u(1,i,j,k) - u(2,i,j,k) &
      )**2 + 1.0/36.0)
	        feq(9)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 + (1.0/12.0)*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(u(2,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(10)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 - 1.0/12.0*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(-u(2,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(11)=  pRho*((1.0/6.0)*u(1,i,j,k)**2 - 1.0/6.0*u(1,i,j,k) - 1.0/12.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(12)=  pRho*(-1.0/12.0*u(1,i,j,k)**2 + (1.0/6.0)*u(2,i,j,k)**2 - 1.0/6.0*u(2,i,j,k) - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(13)=  pRho*(-1.0/12.0*u(1,i,j,k)**2 - 1.0/12.0*u(2,i,j,k)**2 + (1.0/6.0)*u(3,i,j,k)**2 - &
      1.0/6.0*u(3,i,j,k) + 1.0/18.0)
	        feq(14)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(-u(1,i,j,k) - &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(15)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(-u(1,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(16)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(-u(1,i,j,k) - &
      u(2,i,j,k))**2 + 1.0/36.0)
	        feq(17)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 + ( &
      1.0/12.0)*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(-u(1,i,j,k) + &
      u(2,i,j,k))**2 + 1.0/36.0)
	        feq(18)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 - 1.0/12.0*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(-u(2,i,j,k) - &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(19)=  pRho*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 + (1.0/12.0)*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(u(2,i,j,k) - u(3,i,j,k) &
      )**2 + 1.0/36.0)
	        SELECT CASE(forcing_model)
	        	CASE('guo')
              	 call get_Guo_forcing(SS,Fv(1,i,j,k),Fv(2,i,j,k),Fv(3,i,j,k),u(1,i,j,k),u(2,i,j,k),u(3,i,j,k), &
              	 	& omega_s,omega_a)
              	CASE('mlga')
                	Call get_mLGA_forcing(SS,Fvlocal(1),Fvlocal(2),Fvlocal(3))
                	u(1,i,j,k) = u(1,i,j,k) + Fvlocal(1)/(2.0d0*pRho)
                	u(2,i,j,k) = u(2,i,j,k) + Fvlocal(2)/(2.0d0*pRho)
                	u(3,i,j,k) = u(3,i,j,k) + Fvlocal(3)/(2.0d0*pRho)
              	CASE('sc')
                	!only correct the velocity 
                	u(1,i,j,k) = utemp(1) + Fvlocal(1)/(2.0d0*pRho)
                	u(2,i,j,k) = utemp(2) + Fvlocal(2)/(2.0d0*pRho)
                	u(3,i,j,k) = utemp(3) + Fvlocal(3)/(2.0d0*pRho)
	        End SELECT
              !compute f_s, f_a, feq_s,feq_a
              feq_s = 0.d0
              feq_a = 0.d0
              f_s = 0.d0
              f_a = 0.d0
              feq_s(1)=feq(1)
              f_s(1)=f(i,j,k,1)
              feq_a(1)=0.d0
              feq_a(1)=0.d0
              Do n=2,10
                feq_s(n)= (feq(n) + feq(n+9))/2.d0
                feq_a(n)= (feq(n) - feq(n+9))/2.d0
                f_s(n) =  (f(i,j,k,n) + f(i,j,k,n+9))/2.d0
                f_a(n) =  (f(i,j,k,n) - f(i,j,k,n+9))/2.d0
              End Do
              Do n=11,19
                feq_s(n) = feq_s(n-9)
                feq_a(n) = -feq_a(n-9)
                f_s(n) = f_s(n-9)
                f_a(n) = -f_a(n-9)
              End Do
    	      f(i,j,k,:) = f(i,j,k,:) - omega_s*(f_s-feq_s) - omega_a*(f_a-feq_a) + SS
              call get_Fpcorr(Fpcorr,i,j,k,rho(i,j,k),poros, periodicity,lz,ly,lx)
             f(i,j,k,:)=f(i,j,k,:)+Fpcorr
	        Do n=2,10
	          fswap = f(i,j,k,n)
	          f(i,j,k,n)=f(i,j,k,n+9)
	          f(i,j,k,n+9)=fswap
	        End Do
	      End If
	    End Do
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
     Contains
     	Subroutine get_Fpcorr(Fpcorr,i,j,k,rho,poros, periodicity,lz,ly,lx)
     		Implicit None
     		Real(kind = 8), Intent(Inout)::Fpcorr(19)
        Real(kind = 8), Intent(In):: rho
     		Real(kind = 8), Intent(In):: poros(lz,ly,lx)
        Integer,Intent(In)::i,j,k,lz,ly,lx,periodicity(3)
     		Integer,parameter, dimension(19)::	ex = (/0, 1, 0, 0, 1, 1, 1, 1, 0, 0, -1, 0, 0, -1, -1, -1, -1, 0, 0/)
			  Integer,parameter, dimension(19)::	ey = (/0, 0, 1, 0, 0, 0, 1, -1, 1, -1, 0, -1, 0, 0, 0, -1, 1, -1, 1/)
			  Integer,parameter, dimension(19)::	ez = (/0, 0, 0, 1, 1, -1, 0, 0, 1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1/)
			  Real(8),parameter, dimension(19)::	w = (/0.333333333333333d0, 0.0555555555555556d0, 0.0555555555555556d0, &
			                                            0.0555555555555556d0, 0.0277777777777778d0, 0.0277777777777778d0, &
			                                             0.0277777777777778d0, 0.0277777777777778d0, 0.0277777777777778d0, &
			                                            0.0277777777777778d0, 0.0555555555555556d0, 0.0555555555555556d0, &
			                                            0.0555555555555556d0, 0.0277777777777778d0, 0.0277777777777778d0, &
			                                            0.0277777777777778d0, 0.0277777777777778d0, 0.0277777777777778d0, &
			                                            0.0277777777777778d0/)
			 Integer::n, nextX,nextY,nextZ
  			Fpcorr= 0.d0
     		Do n  = 2, 19
     			nextZ = i-ez(n)
     			nextY = j-ey(n)
     			nextX = k+ex(n)
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
                If (nextX > 0 .And. nextX < lx+1 &
 	          & .And. nextY > 0  .And. nextY < ly+1 &
	          & .And. nextZ > 0  .And. nextZ < lz+1) Then
     				Fpcorr(n)=  w(n)* (poros(nextZ,nextY,nextX)-poros(i,j,k)) *rho
     			End If
     		End Do
     	End Subroutine
        Subroutine get_mLGA_forcing(SS,Fx,Fy,Fz)
          Implicit None
          Real(kind = 8), Intent(InOut):: SS(19)
          Real(kind = 8), Intent(In):: Fx, Fy, Fz
          SS(1) = 0
          SS(2) = (1.0/6.0)*Fx
          SS(3) = (1.0/6.0)*Fy
          SS(4) = (1.0/6.0)*Fz
          SS(5) = (1.0/12.0)*Fx + (1.0/12.0)*Fz
          SS(6) = (1.0/12.0)*Fx - 1.0/12.0*Fz
          SS(7) = (1.0/12.0)*Fx + (1.0/12.0)*Fy
          SS(8) = (1.0/12.0)*Fx - 1.0/12.0*Fy
          SS(9) = (1.0/12.0)*Fy + (1.0/12.0)*Fz
          SS(10) = -1.0/12.0*Fy + (1.0/12.0)*Fz
          SS(11) = -1.0/6.0*Fx
          SS(12) = -1.0/6.0*Fy
          SS(13) = -1.0/6.0*Fz
          SS(14) = -1.0/12.0*Fx - 1.0/12.0*Fz
          SS(15) = -1.0/12.0*Fx + (1.0/12.0)*Fz
          SS(16) = -1.0/12.0*Fx - 1.0/12.0*Fy
          SS(17) = -1.0/12.0*Fx + (1.0/12.0)*Fy
          SS(18) = -1.0/12.0*Fy - 1.0/12.0*Fz
          SS(19) = (1.0/12.0)*Fy - 1.0/12.0*Fz
        End Subroutine
        Subroutine get_Guo_forcing(SS,Fx,Fy,Fz,ux,uy,uz,omega_s,omega_a)
        	Implicit None
        	Real(8), Intent(InOut):: SS(19)
        	Real(8), Intent(In)::Fx, Fy, Fz, ux, uy, uz, omega_s, omega_a
        	Real(8):: SS_a(19), SS_s(19)
        	Integer:: n
			 SS(1) = -1.0d0*Fx*ux - 1.0d0*Fy*uy - 1.0d0*Fz*uz
			 SS(2) = 0.333333333333333d0*Fx*ux + 0.166666666666667d0*Fx - 0.166666666666667d0 &
			      *Fy*uy - 0.166666666666667d0*Fz*uz
			 SS(3) = -0.166666666666667d0*Fx*ux + 0.333333333333333d0*Fy*uy + &
			      0.166666666666667d0*Fy - 0.166666666666667d0*Fz*uz
			 SS(4) = -0.166666666666667d0*Fx*ux - 0.166666666666667d0*Fy*uy + &
			      0.333333333333333d0*Fz*uz + 0.166666666666667d0*Fz
			 SS(5) = -0.0833333333333333d0*Fx*ux + 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fz*uz + &
			      0.0833333333333333d0*Fz + 0.25d0*(Fx + Fz)*(ux + uz)
			 SS(6) = -0.0833333333333333d0*Fx*ux + 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fz*uz - &
			      0.0833333333333333d0*Fz + 0.25d0*(Fx - Fz)*(ux - uz)
			 SS(7) = -0.0833333333333333d0*Fx*ux + 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy + 0.0833333333333333d0*Fy - &
			      0.0833333333333333d0*Fz*uz + 0.25d0*(Fx + Fy)*(ux + uy)
			 SS(8) = -0.0833333333333333d0*Fx*ux + 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fy - &
			      0.0833333333333333d0*Fz*uz + 0.25d0*(Fx - Fy)*(ux - uy)
			 SS(9) = -0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fy*uy + &
			      0.0833333333333333d0*Fy - 0.0833333333333333d0*Fz*uz + &
			      0.0833333333333333d0*Fz + 0.25d0*(Fy + Fz)*(uy + uz)
			 SS(10) = -0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fy*uy - &
			      0.0833333333333333d0*Fy - 0.0833333333333333d0*Fz*uz + &
			      0.0833333333333333d0*Fz + 0.25d0*(-Fy + Fz)*(-uy + uz)
			 SS(11) = 0.333333333333333d0*Fx*ux - 0.166666666666667d0*Fx - 0.166666666666667d0 &
			      *Fy*uy - 0.166666666666667d0*Fz*uz
			 SS(12) = -0.166666666666667d0*Fx*ux + 0.333333333333333d0*Fy*uy - &
			      0.166666666666667d0*Fy - 0.166666666666667d0*Fz*uz
			 SS(13) = -0.166666666666667d0*Fx*ux - 0.166666666666667d0*Fy*uy + &
			      0.333333333333333d0*Fz*uz - 0.166666666666667d0*Fz
			 SS(14) = -0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fz*uz - &
			      0.0833333333333333d0*Fz + 0.25d0*(-Fx - Fz)*(-ux - uz)
			 SS(15) = -0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fz*uz + &
			      0.0833333333333333d0*Fz + 0.25d0*(-Fx + Fz)*(-ux + uz)
			 SS(16) = -0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy - 0.0833333333333333d0*Fy - &
			      0.0833333333333333d0*Fz*uz + 0.25d0*(-Fx - Fy)*(-ux - uy)
			 SS(17) = -0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fx - &
			      0.0833333333333333d0*Fy*uy + 0.0833333333333333d0*Fy - &
			      0.0833333333333333d0*Fz*uz + 0.25d0*(-Fx + Fy)*(-ux + uy)
			 SS(18) = -0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fy*uy - &
			      0.0833333333333333d0*Fy - 0.0833333333333333d0*Fz*uz - &
			      0.0833333333333333d0*Fz + 0.25d0*(-Fy - Fz)*(-uy - uz)
			 SS(19) = -0.0833333333333333d0*Fx*ux - 0.0833333333333333d0*Fy*uy + &
			      0.0833333333333333d0*Fy - 0.0833333333333333d0*Fz*uz - &
			      0.0833333333333333d0*Fz + 0.25d0*(Fy - Fz)*(uy - uz)
        	SS_s = 0.0d0
        	SS_a = 0.0d0
        	SS_s(1) =  SS(1)
        	Do n = 2,10
          		SS_s(n) = 0.5*(SS(n) + SS(n+9))
          		SS_a(n) = 0.5*(SS(n) - SS(n+9))
        	End Do
        	Do n = 11,19
          		SS_s(n) = SS_s(n-9)
          		SS_a(n) = -SS_a(n-9)
        	End Do
        	SS = (1.0d0-0.5d0*omega_s)*SS_s + (1.0d0-0.5d0*omega_a)*SS_a  
        End Subroutine
End Subroutine collide_trt
!
Subroutine  stream_and_bounce(f, nodetype, periodicity, lz, ly, lx)
!subroutine for propogation step using latt(2007) swap algorithim
!for D3Q19 lattices
	Implicit None
	Real(kind = 8), Intent(Inout) :: f(lz,ly,lx,19)	!Distribution function
	Real(kind = 8), Intent(In) :: nodetype(lz,ly,lx)	!nodetype
     Integer, Intent(In) :: periodicity(3) !flag for periodicity 0 dimension for x, 1 dimension for y and 2 dimension for z 
	Integer, Intent(In) :: lz, ly, lx	!length in y and x direction respectively
	Integer :: ex(19)	!lattice velocity in x direction
	Integer :: ey(19)	!lattice velocity in y direction
	Integer :: ez(19)	!lattice velocity in z direction
	Integer :: i, j, k, n, half, nextX, nextY, nextZ
	Real(kind = 8) :: ftemp
	ex = (/0, 1, 0, 0, 1, 1, 1, 1, 0, 0, -1, 0, 0, -1, -1, -1, -1, 0, 0/)
	ey = (/0, 0, 1, 0, 0, 0, 1, -1, 1, -1, 0, -1, 0, 0, 0, -1, 1, -1, 1/)
	ez = (/0, 0, 0, 1, 1, -1, 0, 0, 1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1/)
	half= 9
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, n, nextX, nextY, nextZ, ftemp)
	!$OMP DO  
	Do k= 1, lx
	  Do j= 1, ly
	    Do i= 1, lz
	      If (nodetype(i,j,k) <= 0) Then
	        Do  n= 2,10
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
                If (nextX > 0 .And. nextX < lx+1 &
 	          & .And. nextY > 0  .And. nextY < ly+1 &
	          & .And. nextZ > 0  .And. nextZ < lz+1) Then
	            If (nodetype(nextZ, nextY, nextX) <= 0 .And. &
	            & nodetype(nextZ,j,k)<=0 .And. nodetype(i,nextY,k)<=0 .And. nodetype(i,j,nextX) <=0) Then
	            !If (nodetype(nextZ, nextY, nextX) <= 0) Then
	              ftemp = f (nextZ, nextY, nextX, n)
	              f(nextZ, nextY, nextX, n) = f (i, j, k, n+half)
	              f(i, j, k, n+half) = ftemp
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
Subroutine apply_bc(f,rho0,poros,nodetype,interp,left_type,left_val,right_type,right_val, &
&   top_type,top_val,bottom_type,bottom_val,front_type,front_val,back_type,back_val,&
&   lz,ly,lx)
    !Subroutine to apply boundary condition for navier stokes equation
    !for d3q19 lattice.  
    Implicit None
    Real(kind = 8), Intent(InOut):: f(lz,ly,lx,19) ! distribution function
    Real(kind = 8), Intent(In):: rho0 ! density of fluid
    Real(kind = 8), Intent(In):: poros(lz,ly,lx) ! density of fluid
    Real(kind = 8), Intent(In):: nodetype(lz,ly,lx) ! nodetype indicator
    Character(Len=*), Intent(In):: left_type ! left boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    Character(Len=*), Intent(In):: right_type ! right boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    Character(Len=*), Intent(In):: top_type ! top boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    character(Len=*), Intent(In):: bottom_type ! bottom boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    character(Len=*), Intent(In):: front_type ! front boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    character(Len=*), Intent(In):: back_type ! back boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    Real(kind = 8), Intent(In):: left_val(3) ! value for left boundary
    Real(kind = 8), Intent(In):: right_val(3) ! value for right boundary
    Real(kind = 8), Intent(In):: top_val(3) ! value for top boundary
    Real(kind = 8), Intent(In):: bottom_val(3) ! value for bottom boundary
    Real(kind = 8), Intent(In):: front_val(3) ! value for front boundary
    Real(kind = 8), Intent(In):: back_val(3) ! value for back boundary    
    Integer, Intent(In):: interp ! 1 to allow interpolation to set boundary condtion 0 to not interpolate at the boundary node
    Integer, Intent(In):: lz ! number of nodes in z direction
    Integer, Intent(In):: ly ! number of nodes in y direction
    Integer, Intent(In):: lx ! number of nodes in x direction
    !other variables...
    Integer:: i, j, k, inext,jnext
    Real(kind = 8):: bc_val1,bc_val2,bc_val3,l,ftemp,rhonbh,ubx,uby,ubz,es2
    es2= 0.333333333333
    !$OMP PARALLEL DEFAULT (SHARED) PRIVATE(i,j,k,inext,jnext,bc_val1,bc_val2,bc_val3,rhonbh,ubx,uby,ubz,l,ftemp)    
    !periodic boundary conditions
    If (left_type== 'periodic' .AND. right_type == 'periodic') Then
      !$OMP DO 
      Do i = 1,lz
        Do j= 1,ly
          inext = i-0
          jnext = j-0
          if (inext > lz) inext=1
          if (inext < 1) inext=lz
          if (jnext > ly) jnext=1
          if (jnext < 1) jnext=ly
          ftemp = f(i,j,1,2)
          f(i,j,1,2) = f(inext,jnext,lx,11) 
          f(inext,jnext,lx,11) = ftemp      
          inext = i+1
          jnext = j-0
          if (inext > lz) inext=1
          if (inext < 1) inext=lz
          if (jnext > ly) jnext=1
          if (jnext < 1) jnext=ly
          ftemp = f(i,j,1,5)
          f(i,j,1,5) = f(inext,jnext,lx,14) 
          f(inext,jnext,lx,14) = ftemp      
          inext = i-1
          jnext = j-0
          if (inext > lz) inext=1
          if (inext < 1) inext=lz
          if (jnext > ly) jnext=1
          if (jnext < 1) jnext=ly
          ftemp = f(i,j,1,6)
          f(i,j,1,6) = f(inext,jnext,lx,15) 
          f(inext,jnext,lx,15) = ftemp      
          inext = i-0
          jnext = j+1
          if (inext > lz) inext=1
          if (inext < 1) inext=lz
          if (jnext > ly) jnext=1
          if (jnext < 1) jnext=ly
          ftemp = f(i,j,1,7)
          f(i,j,1,7) = f(inext,jnext,lx,16) 
          f(inext,jnext,lx,16) = ftemp      
          inext = i-0
          jnext = j-1
          if (inext > lz) inext=1
          if (inext < 1) inext=lz
          if (jnext > ly) jnext=1
          if (jnext < 1) jnext=ly
          ftemp = f(i,j,1,8)
          f(i,j,1,8) = f(inext,jnext,lx,17) 
          f(inext,jnext,lx,17) = ftemp      
        End Do
      End Do
      !$OMP END DO 
    End if
    If (top_type== 'periodic' .AND. bottom_type == 'periodic') Then
      !$OMP DO 
      Do i = 1,lz
        Do j = 1,lx
          inext = i-0
          jnext = j+0
          if (inext > lz) inext=1
          if (inext < 1) inext=lz
          if (jnext > lx) jnext=1
          if (jnext < 1) jnext=lx
          ftemp = f(i,1,j,12)
          f(i,1,j,12) = f(inext,ly,jnext,3) 
          f(inext,ly,jnext,3) = ftemp      
          inext = i-0
          jnext = j+1
          if (inext > lz) inext=1
          if (inext < 1) inext=lz
          if (jnext > lx) jnext=1
          if (jnext < 1) jnext=lx
          ftemp = f(i,1,j,16)
          f(i,1,j,16) = f(inext,ly,jnext,7) 
          f(inext,ly,jnext,7) = ftemp      
          inext = i-1
          jnext = j+0
          if (inext > lz) inext=1
          if (inext < 1) inext=lz
          if (jnext > lx) jnext=1
          if (jnext < 1) jnext=lx
          ftemp = f(i,1,j,18)
          f(i,1,j,18) = f(inext,ly,jnext,9) 
          f(inext,ly,jnext,9) = ftemp      
          inext = i-0
          jnext = j-1
          if (inext > lz) inext=1
          if (inext < 1) inext=lz
          if (jnext > lx) jnext=1
          if (jnext < 1) jnext=lx
          ftemp = f(i,1,j,8)
          f(i,1,j,8) = f(inext,ly,jnext,17) 
          f(inext,ly,jnext,17) = ftemp      
          inext = i+1
          jnext = j+0
          if (inext > lz) inext=1
          if (inext < 1) inext=lz
          if (jnext > lx) jnext=1
          if (jnext < 1) jnext=lx
          ftemp = f(i,1,j,10)
          f(i,1,j,10) = f(inext,ly,jnext,19) 
          f(inext,ly,jnext,19) = ftemp      
        End Do
      End Do
      !$OMP END DO 
    End if
    If (front_type== 'periodic' .AND. back_type == 'periodic') Then
      !$OMP DO 
      Do i = 1,ly
        Do j = 1,lx
          inext = i-0
          jnext = j+0
          if (inext > ly) inext=1
          if (inext < 1) inext=ly
          if (jnext > lx) jnext=1
          if (jnext < 1) jnext=lx
          ftemp = f(1,i,j,13)
          f(1,i,j,13) = f(lz,inext,jnext,4) 
          f(lz,inext,jnext,4) = ftemp      
          inext = i-0
          jnext = j+1
          if (inext > ly) inext=1
          if (inext < 1) inext=ly
          if (jnext > lx) jnext=1
          if (jnext < 1) jnext=lx
          ftemp = f(1,i,j,14)
          f(1,i,j,14) = f(lz,inext,jnext,5) 
          f(lz,inext,jnext,5) = ftemp      
          inext = i-1
          jnext = j+0
          if (inext > ly) inext=1
          if (inext < 1) inext=ly
          if (jnext > lx) jnext=1
          if (jnext < 1) jnext=lx
          ftemp = f(1,i,j,18)
          f(1,i,j,18) = f(lz,inext,jnext,9) 
          f(lz,inext,jnext,9) = ftemp      
          inext = i+1
          jnext = j+0
          if (inext > ly) inext=1
          if (inext < 1) inext=ly
          if (jnext > lx) jnext=1
          if (jnext < 1) jnext=lx
          ftemp = f(1,i,j,19)
          f(1,i,j,19) = f(lz,inext,jnext,10) 
          f(lz,inext,jnext,10) = ftemp      
          inext = i-0
          jnext = j-1
          if (inext > ly) inext=1
          if (inext < 1) inext=ly
          if (jnext > lx) jnext=1
          if (jnext < 1) jnext=lx
          ftemp = f(1,i,j,6)
          f(1,i,j,6) = f(lz,inext,jnext,15) 
          f(lz,inext,jnext,15) = ftemp      
        End Do
      End Do
      !$OMP END DO 
    End if
    select case (left_type)
     case ('p')
        bc_val1 = left_val(1)/es2 + rho0
       !$OMP DO 
        Do j= 1,ly
          Do i=1,lz
            If (nodetype(i,j,1)<=0) Then
              !convert pressure to density
              If (interp==1) Then
                rhonbh = sum(f(i,j,1,:)) 
                bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
              End If 
              call pressure_extrapol_bc(f(i,j,1,:),f(i,j,2,:),bc_val1,poros(i,j,1),poros(i,j,2))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = left_val(1)
        bc_val2 = left_val(2)
        bc_val3 = left_val(3)
        !$OMP DO 
        Do j= 1, ly
          Do i=1, lz
            If (nodetype(i,j,1)<=0) Then
              If (interp==1) Then
                call left_vel_bb_bc(f(i,j,1,:),rho0,bc_val1,bc_val2,bc_val3,poros(i,j,1))
              Else 
                call vel_extrapol_bc(f(i,j,1,:),f(i,j,2,:),bc_val1,bc_val2,bc_val3,poros(i,j,1),poros(i,j,2))
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = left_val(1)
        bc_val2 = left_val(2)
        bc_val3 = left_val(3)
        if ((nint(lz*bc_val2)-lz*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(lz*bc_val2) + 1
          else
            l =nint(ly*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(lz*bc_val2)
          else
            l=nint(ly*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i=1, lz
            If (nodetype(i,j,1)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
                uby= 0.0d0
                ubz= 0.0d0
                call left_vel_bb_bc(f(i,j,1,:),rho0,bc_val1,bc_val2,bc_val3,poros(i,j,1))
              Else
                ubx= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
                uby= 0.0d0
                ubz= 0.0d0
                call vel_extrapol_bc(f(i,j,1,:),f(i,j,2,:),bc_val1,bc_val2,bc_val3,poros(i,j,1),poros(i,j,2))          
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i= 1,lz
            If (nodetype(i,j,1)<=0) Then
              f(i,j,1,:) = f(i,j,2,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    select case (right_type)
     case ('p')
        bc_val1 = right_val(1)/es2 + rho0
       !$OMP DO 
        Do j= 1,ly
          Do i=1,lz
            If (nodetype(i,j,lx)<=0) Then
              !convert pressure to density
              If (interp==1) Then
                rhonbh = sum(f(i,j,lx,:)) 
                bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
              End If 
              call pressure_extrapol_bc(f(i,j,lx,:),f(i,j,lx-1,:),bc_val1,poros(i,j,lx),poros(i,j,lx-1))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = right_val(1)
        bc_val2 = right_val(2)
        bc_val3 = right_val(3)
        !$OMP DO 
        Do j= 1, ly
          Do i=1, lz
            If (nodetype(i,j,lx)<=0) Then
              If (interp==1) Then
                call right_vel_bb_bc(f(i,j,lx,:),rho0,bc_val1,bc_val2,bc_val3,poros(i,j,lx))
              Else 
                call vel_extrapol_bc(f(i,j,lx,:),f(i,j,lx-1,:),bc_val1,bc_val2,bc_val3,poros(i,j,lx),poros(i,j,lx-1))
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = right_val(1)
        bc_val2 = right_val(2)
        bc_val3 = right_val(3)
        if ((nint(lz*bc_val2)-lz*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(lz*bc_val2) + 1
          else
            l =nint(ly*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(lz*bc_val2)
          else
            l=nint(ly*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i=1, lz
            If (nodetype(i,j,lx)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
                uby= 0.0d0
                ubz= 0.0d0
                call right_vel_bb_bc(f(i,j,lx,:),rho0,bc_val1,bc_val2,bc_val3,poros(i,j,lx))
              Else
                ubx= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
                uby= 0.0d0
                ubz= 0.0d0
                call vel_extrapol_bc(f(i,j,lx,:),f(i,j,lx-1,:),bc_val1,bc_val2,bc_val3,poros(i,j,lx),poros(i,j,lx-1))          
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i= 1,lz
            If (nodetype(i,j,lx)<=0) Then
              f(i,j,lx,:) = f(i,j,lx-1,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    select case (top_type)
     case ('p')
        bc_val1 = top_val(1)/es2 + rho0
       !$OMP DO 
        Do j= 1,lx
          Do i=1,lz
            If (nodetype(i,1,j)<=0) Then
              !convert pressure to density
              If (interp==1) Then
                rhonbh = sum(f(i,1,j,:)) 
                bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
              End If 
              call pressure_extrapol_bc(f(i,1,j,:),f(i,2,j,:),bc_val1,poros(i,1,j),poros(i,2,j))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = top_val(1)
        bc_val2 = top_val(2)
        bc_val3 = top_val(3)
        !$OMP DO 
        Do j= 1, lx
          Do i=1, lz
            If (nodetype(i,1,j)<=0) Then
              If (interp==1) Then
                call top_vel_bb_bc(f(i,1,j,:),rho0,bc_val1,bc_val2,bc_val3,poros(i,1,j))
              Else 
                call vel_extrapol_bc(f(i,1,j,:),f(i,2,j,:),bc_val1,bc_val2,bc_val3,poros(i,1,j),poros(i,2,j))
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = top_val(1)
        bc_val2 = top_val(2)
        bc_val3 = top_val(3)
        if ((nint(lz*bc_val2)-lz*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(lz*bc_val2) + 1
          else
            l =nint(lx*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(lz*bc_val2)
          else
            l=nint(lx*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i=1, lz
            If (nodetype(i,1,j)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= 0.0d0
                uby= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
                ubz= 0.0d0
                call top_vel_bb_bc(f(i,1,j,:),rho0,bc_val1,bc_val2,bc_val3,poros(i,1,j))
              Else
                ubx= 0.0d0
                uby= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
                ubz= 0.0d0
                call vel_extrapol_bc(f(i,1,j,:),f(i,2,j,:),bc_val1,bc_val2,bc_val3,poros(i,1,j),poros(i,2,j))          
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i= 1,lz
            If (nodetype(i,1,j)<=0) Then
              f(i,1,j,:) = f(i,2,j,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    select case (bottom_type)
     case ('p')
        bc_val1 = bottom_val(1)/es2 + rho0
       !$OMP DO 
        Do j= 1,lx
          Do i=1,lz
            If (nodetype(i,ly,j)<=0) Then
              !convert pressure to density
              If (interp==1) Then
                rhonbh = sum(f(i,ly,j,:)) 
                bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
              End If 
              call pressure_extrapol_bc(f(i,ly,j,:),f(i,ly-1,j,:),bc_val1,poros(i,ly,j),poros(i,ly-1,j))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = bottom_val(1)
        bc_val2 = bottom_val(2)
        bc_val3 = bottom_val(3)
        !$OMP DO 
        Do j= 1, lx
          Do i=1, lz
            If (nodetype(i,ly,j)<=0) Then
              If (interp==1) Then
                call bottom_vel_bb_bc(f(i,ly,j,:),rho0,bc_val1,bc_val2,bc_val3,poros(i,ly,j))
              Else 
                call vel_extrapol_bc(f(i,ly,j,:),f(i,ly-1,j,:),bc_val1,bc_val2,bc_val3,poros(i,ly,j),poros(i,ly-1,j))
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = bottom_val(1)
        bc_val2 = bottom_val(2)
        bc_val3 = bottom_val(3)
        if ((nint(lz*bc_val2)-lz*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(lz*bc_val2) + 1
          else
            l =nint(lx*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(lz*bc_val2)
          else
            l=nint(lx*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i=1, lz
            If (nodetype(i,ly,j)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= 0.0d0
                uby= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
                ubz= 0.0d0
                call bottom_vel_bb_bc(f(i,ly,j,:),rho0,bc_val1,bc_val2,bc_val3,poros(i,ly,j))
              Else
                ubx= 0.0d0
                uby= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
                ubz= 0.0d0
                call vel_extrapol_bc(f(i,ly,j,:),f(i,ly-1,j,:),bc_val1,bc_val2,bc_val3,poros(i,ly,j),poros(i,ly-1,j))          
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i= 1,lz
            If (nodetype(i,ly,j)<=0) Then
              f(i,ly,j,:) = f(i,ly-1,j,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    select case (front_type)
     case ('p')
        bc_val1 = front_val(1)/es2 + rho0
       !$OMP DO 
        Do j= 1,ly
          Do i=1,ly
            If (nodetype(1,i,j)<=0) Then
              !convert pressure to density
              If (interp==1) Then
                rhonbh = sum(f(1,i,j,:)) 
                bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
              End If 
              call pressure_extrapol_bc(f(1,i,j,:),f(2,i,j,:),bc_val1,poros(1,i,j),poros(2,i,j))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = front_val(1)
        bc_val2 = front_val(2)
        bc_val3 = front_val(3)
        !$OMP DO 
        Do j= 1, ly
          Do i=1, ly
            If (nodetype(1,i,j)<=0) Then
              If (interp==1) Then
                call front_vel_bb_bc(f(1,i,j,:),rho0,bc_val1,bc_val2,bc_val3,poros(1,i,j))
              Else 
                call vel_extrapol_bc(f(1,i,j,:),f(2,i,j,:),bc_val1,bc_val2,bc_val3,poros(1,i,j),poros(2,i,j))
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = front_val(1)
        bc_val2 = front_val(2)
        bc_val3 = front_val(3)
        if ((nint(ly*bc_val2)-ly*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(ly*bc_val2) + 1
          else
            l =nint(ly*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(ly*bc_val2)
          else
            l=nint(ly*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i=1, ly
            If (nodetype(1,i,j)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= 0.0d0
                uby= 0.0d0
                ubz= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
                call front_vel_bb_bc(f(1,i,j,:),rho0,bc_val1,bc_val2,bc_val3,poros(1,i,j))
              Else
                ubx= 0.0d0
                uby= 0.0d0
                ubz= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
                call vel_extrapol_bc(f(1,i,j,:),f(2,i,j,:),bc_val1,bc_val2,bc_val3,poros(1,i,j),poros(2,i,j))          
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i= 1,ly
            If (nodetype(1,i,j)<=0) Then
              f(1,i,j,:) = f(2,i,j,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    select case (back_type)
     case ('p')
        bc_val1 = back_val(1)/es2 + rho0
       !$OMP DO 
        Do j= 1,lx
          Do i=1,lx
            If (nodetype(lz,i,j)<=0) Then
              !convert pressure to density
              If (interp==1) Then
                rhonbh = sum(f(lz,i,j,:)) 
                bc_val1 = (bc_val1 + 0.5d0*rhonbh)/1.5d0
              End If 
              call pressure_extrapol_bc(f(lz,i,j,:),f(lz-1,i,j,:),bc_val1,poros(lz,i,j),poros(lz-1,i,j))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = back_val(1)
        bc_val2 = back_val(2)
        bc_val3 = back_val(3)
        !$OMP DO 
        Do j= 1, lx
          Do i=1, lx
            If (nodetype(lz,i,j)<=0) Then
              If (interp==1) Then
                call back_vel_bb_bc(f(lz,i,j,:),rho0,bc_val1,bc_val2,bc_val3,poros(lz,i,j))
              Else 
                call vel_extrapol_bc(f(lz,i,j,:),f(lz-1,i,j,:),bc_val1,bc_val2,bc_val3,poros(lz,i,j),poros(lz-1,i,j))
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = back_val(1)
        bc_val2 = back_val(2)
        bc_val3 = back_val(3)
        if ((nint(lx*bc_val2)-lx*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(lx*bc_val2) + 1
          else
            l =nint(lx*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(lx*bc_val2)
          else
            l=nint(lx*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i=1, lx
            If (nodetype(lz,i,j)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= 0.0d0
                uby= 0.0d0
                ubz= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
                call back_vel_bb_bc(f(lz,i,j,:),rho0,bc_val1,bc_val2,bc_val3,poros(lz,i,j))
              Else
                ubx= 0.0d0
                uby= 0.0d0
                ubz= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
                call vel_extrapol_bc(f(lz,i,j,:),f(lz-1,i,j,:),bc_val1,bc_val2,bc_val3,poros(lz,i,j),poros(lz-1,i,j))          
              End If
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i= 1,lx
            If (nodetype(lz,i,j)<=0) Then
              f(lz,i,j,:) = f(lz-1,i,j,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    !$OMP END PARALLEL
    Contains
    Subroutine pressure_extrapol_bc(f,fnbh,rho,poros,porosnbh)
      Implicit None
      Real(kind = 8), Intent(In):: fnbh(19),rho,poros,porosnbh
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8) :: feq(19),feqnbh(19),rhonbh,uxnbh,uynbh,uznbh
      rhonbh = sum(fnbh)/porosnbh
      call get_vel(uxnbh,uynbh,uznbh,fnbh,rhonbh,porosnbh)
      call get_feq(feq,rho,uxnbh,uynbh,uznbh,poros)
      call get_feq(feqnbh,rhonbh,uxnbh,uynbh,uznbh,porosnbh)
      f= feq + (fnbh-feqnbh)
    End Subroutine pressure_extrapol_bc
    Subroutine vel_extrapol_bc(f,fnbh,ux,uy,uz,poros,porosnbh)
      Implicit None
      Real(kind = 8), Intent(In):: fnbh(19),ux,uy,uz,poros,porosnbh
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8) :: feq(19),feqnbh(19),rhonbh,uxnbh,uynbh,uznbh
      rhonbh = sum(fnbh)/porosnbh
      call get_vel(uxnbh,uynbh,uznbh,fnbh,rhonbh,porosnbh)
      call get_feq(feq,rhonbh,ux,uy,uz,poros)
      call get_feq(feqnbh,rhonbh,uxnbh,uynbh,uznbh,porosnbh)
      f= feq + (fnbh-feqnbh)
    End Subroutine vel_extrapol_bc
    Subroutine left_vel_bb_bc(f,rho,ux,uy,uz,poros)
      Real(kind = 8), Intent(In):: rho,ux,uy,uz,poros
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8):: rhoPoros
      rhoPoros = rho * poros
      f(2)= f(2)-(-1.0d0/3.0d0*rhoPoros*ux)
      f(5)= f(5)-((1.0d0/6.0d0)*rhoPoros*(-ux - uz))
      f(6)= f(6)-((1.0d0/6.0d0)*rhoPoros*(-ux + uz))
      f(7)= f(7)-((1.0d0/6.0d0)*rhoPoros*(-ux - uy))
      f(8)= f(8)-((1.0d0/6.0d0)*rhoPoros*(-ux + uy))
    End Subroutine left_vel_bb_bc
    Subroutine right_vel_bb_bc(f,rho,ux,uy,uz,poros)
      Real(kind = 8), Intent(In):: rho,ux,uy,uz,poros
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8):: rhoPoros
      rhoPoros = rho * poros
      f(11)= f(11)-((1.0d0/3.0d0)*rhoPoros*ux)
      f(14)= f(14)-((1.0d0/6.0d0)*rhoPoros*(ux + uz))
      f(15)= f(15)-((1.0d0/6.0d0)*rhoPoros*(ux - uz))
      f(16)= f(16)-((1.0d0/6.0d0)*rhoPoros*(ux + uy))
      f(17)= f(17)-((1.0d0/6.0d0)*rhoPoros*(ux - uy))
    End Subroutine right_vel_bb_bc
    Subroutine top_vel_bb_bc(f,rho,ux,uy,uz,poros)
      Real(kind = 8), Intent(In):: rho,ux,uy,uz,poros
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8):: rhoPoros
      rhoPoros = rho * poros
      f(12)= f(12)-((1.0d0/3.0d0)*rhoPoros*uy)
      f(16)= f(16)-((1.0d0/6.0d0)*rhoPoros*(ux + uy))
      f(18)= f(18)-((1.0d0/6.0d0)*rhoPoros*(uy + uz))
      f(8)= f(8)-((1.0d0/6.0d0)*rhoPoros*(-ux + uy))
      f(10)= f(10)-((1.0d0/6.0d0)*rhoPoros*(uy - uz))
    End Subroutine top_vel_bb_bc
    Subroutine bottom_vel_bb_bc(f,rho,ux,uy,uz,poros)
      Real(kind = 8), Intent(In):: rho,ux,uy,uz,poros
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8):: rhoPoros
      rhoPoros = rho * poros
      f(17)= f(17)-((1.0d0/6.0d0)*rhoPoros*(ux - uy))
      f(19)= f(19)-((1.0d0/6.0d0)*rhoPoros*(-uy + uz))
      f(3)= f(3)-(-1.0d0/3.0d0*rhoPoros*uy)
      f(7)= f(7)-((1.0d0/6.0d0)*rhoPoros*(-ux - uy))
      f(9)= f(9)-((1.0d0/6.0d0)*rhoPoros*(-uy - uz))
    End Subroutine bottom_vel_bb_bc
    Subroutine front_vel_bb_bc(f,rho,ux,uy,uz,poros)
      Real(kind = 8), Intent(In):: rho,ux,uy,uz,poros
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8):: rhoPoros
      rhoPoros = rho * poros
      f(13)= f(13)-((1.0d0/3.0d0)*rhoPoros*uz)
      f(14)= f(14)-((1.0d0/6.0d0)*rhoPoros*(ux + uz))
      f(18)= f(18)-((1.0d0/6.0d0)*rhoPoros*(uy + uz))
      f(19)= f(19)-((1.0d0/6.0d0)*rhoPoros*(-uy + uz))
      f(6)= f(6)-((1.0d0/6.0d0)*rhoPoros*(-ux + uz))
    End Subroutine front_vel_bb_bc
    Subroutine back_vel_bb_bc(f,rho,ux,uy,uz,poros)
      Real(kind = 8), Intent(In):: rho,ux,uy,uz,poros
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8):: rhoPoros
      rhoPoros = rho * poros
      f(15)= f(15)-((1.0d0/6.0d0)*rhoPoros*(ux - uz))
      f(4)= f(4)-(-1.0d0/3.0d0*rhoPoros*uz)
      f(5)= f(5)-((1.0d0/6.0d0)*rhoPoros*(-ux - uz))
      f(9)= f(9)-((1.0d0/6.0d0)*rhoPoros*(-uy - uz))
      f(10)= f(10)-((1.0d0/6.0d0)*rhoPoros*(uy - uz))
    End Subroutine back_vel_bb_bc
    Subroutine get_vel(ux,uy,uz,f,rho,poros)
      Implicit None
      Real(kind = 8), Intent(In)::  f(19),rho,poros
      Real(kind = 8), Intent(InOut)::ux,uy,uz
      Real(kind = 8):: m1x,m1y,m1z
      m1x= -f(11) - f(14) - f(15) - f(16) - f(17) + f(2) + f(5) + f(6) + f(7) + f(8)
      m1y= -f(10) - f(12) - f(16) + f(17) - f(18) + f(19) + f(3) + f(7) - f(8) + f(9)
      m1z= f(10) - f(13) - f(14) + f(15) - f(18) - f(19) + f(4) + f(5) - f(6) + f(9)
      ux= m1x/(rho*poros)
      uy= m1y/(rho*poros)
      uz= m1z/(rho*poros)
    End Subroutine get_vel
    Subroutine get_feq(feq,rho,ux,uy,uz,poros)
      Implicit None
      Real(kind = 8), Intent(In)::  rho,ux,uy,uz,poros
      Real(kind = 8), Intent(InOut)::feq(19)
      Real(kind = 8):: rhoPoros
      rhoPoros  = rho * poros
      feq(1)=rhoPoros*(-0.5d0*ux**2 - 0.5d0*uy**2 - 0.5d0*uz**2 + 0.333333333333333d0)
      feq(2)=rhoPoros*(0.166666666666667d0*ux**2 + 0.166666666666667d0*ux - &
      0.0833333333333333d0*uy**2 - 0.0833333333333333d0*uz**2 + &
      0.0555555555555556d0)
      feq(3)=rhoPoros*(-0.0833333333333333d0*ux**2 + 0.166666666666667d0*uy**2 + &
      0.166666666666667d0*uy - 0.0833333333333333d0*uz**2 + &
      0.0555555555555556d0)
      feq(4)=rhoPoros*(-0.0833333333333333d0*ux**2 - 0.0833333333333333d0*uy**2 + &
      0.166666666666667d0*uz**2 + 0.166666666666667d0*uz + &
      0.0555555555555556d0)
      feq(5)=rhoPoros*(-0.0416666666666667d0*ux**2 + 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 - 0.0416666666666667d0*uz**2 + &
      0.0833333333333333d0*uz + 0.125d0*(ux + uz)**2 + &
      0.0277777777777778d0)
      feq(6)=rhoPoros*(-0.0416666666666667d0*ux**2 + 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 - 0.0416666666666667d0*uz**2 - &
      0.0833333333333333d0*uz + 0.125d0*(ux - uz)**2 + &
      0.0277777777777778d0)
      feq(7)=rhoPoros*(-0.0416666666666667d0*ux**2 + 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 + 0.0833333333333333d0*uy - &
      0.0416666666666667d0*uz**2 + 0.125d0*(ux + uy)**2 + &
      0.0277777777777778d0)
      feq(8)=rhoPoros*(-0.0416666666666667d0*ux**2 + 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 - 0.0833333333333333d0*uy - &
      0.0416666666666667d0*uz**2 + 0.125d0*(ux - uy)**2 + &
      0.0277777777777778d0)
      feq(9)=rhoPoros*(-0.0416666666666667d0*ux**2 - 0.0416666666666667d0*uy**2 + &
      0.0833333333333333d0*uy - 0.0416666666666667d0*uz**2 + &
      0.0833333333333333d0*uz + 0.125d0*(uy + uz)**2 + &
      0.0277777777777778d0)
      feq(10)=rhoPoros*(-0.0416666666666667d0*ux**2 - 0.0416666666666667d0*uy**2 - &
      0.0833333333333333d0*uy - 0.0416666666666667d0*uz**2 + &
      0.0833333333333333d0*uz + 0.125d0*(-uy + uz)**2 + &
      0.0277777777777778d0)
      feq(11)=rhoPoros*(0.166666666666667d0*ux**2 - 0.166666666666667d0*ux - &
      0.0833333333333333d0*uy**2 - 0.0833333333333333d0*uz**2 + &
      0.0555555555555556d0)
      feq(12)=rhoPoros*(-0.0833333333333333d0*ux**2 + 0.166666666666667d0*uy**2 - &
      0.166666666666667d0*uy - 0.0833333333333333d0*uz**2 + &
      0.0555555555555556d0)
      feq(13)=rhoPoros*(-0.0833333333333333d0*ux**2 - 0.0833333333333333d0*uy**2 + &
      0.166666666666667d0*uz**2 - 0.166666666666667d0*uz + &
      0.0555555555555556d0)
      feq(14)=rhoPoros*(-0.0416666666666667d0*ux**2 - 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 - 0.0416666666666667d0*uz**2 - &
      0.0833333333333333d0*uz + 0.125d0*(-ux - uz)**2 + &
      0.0277777777777778d0)
      feq(15)=rhoPoros*(-0.0416666666666667d0*ux**2 - 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 - 0.0416666666666667d0*uz**2 + &
      0.0833333333333333d0*uz + 0.125d0*(-ux + uz)**2 + &
      0.0277777777777778d0)
      feq(16)=rhoPoros*(-0.0416666666666667d0*ux**2 - 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 - 0.0833333333333333d0*uy - &
      0.0416666666666667d0*uz**2 + 0.125d0*(-ux - uy)**2 + &
      0.0277777777777778d0)
      feq(17)=rhoPoros*(-0.0416666666666667d0*ux**2 - 0.0833333333333333d0*ux - &
      0.0416666666666667d0*uy**2 + 0.0833333333333333d0*uy - &
      0.0416666666666667d0*uz**2 + 0.125d0*(-ux + uy)**2 + &
      0.0277777777777778d0)
      feq(18)=rhoPoros*(-0.0416666666666667d0*ux**2 - 0.0416666666666667d0*uy**2 - &
      0.0833333333333333d0*uy - 0.0416666666666667d0*uz**2 - &
      0.0833333333333333d0*uz + 0.125d0*(-uy - uz)**2 + &
      0.0277777777777778d0)
      feq(19)=rhoPoros*(-0.0416666666666667d0*ux**2 - 0.0416666666666667d0*uy**2 + &
      0.0833333333333333d0*uy - 0.0416666666666667d0*uz**2 - &
      0.0833333333333333d0*uz + 0.125d0*(uy - uz)**2 + &
      0.0277777777777778d0)
    End Subroutine get_feq
End subroutine 


