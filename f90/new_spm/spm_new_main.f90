! SPM implementation for the simulation of a zero dimensional reactor.
! This implements the version of the SPM found in Sasha's SPM Paper.
! Alastair J. Smith (ajs224@cam.ac.uk)
! 20/06/04

PROGRAM spm_new_main

  USE global_data

  USE nrtype

  IMPLICIT NONE
  
  REAL(DP)::z(N) ! Field particle masses
  REAL(DP)::x ! test particle mass
  REAL(DP)::z_min
  INTEGER::i,j,event
  REAL(DP)::Phi ! Total mass density of field particles
  REAL(DP)::M ! Total mass density of field particles

  REAL(DP)::rate(6)
  REAL(DP)::p(6)
  INTEGER::events(6)
  REAL(DP)::rate_total
  REAL(DP)::t,tau,t_stop
  REAL(DP)::tau_N
 
  REAL(DP)::coll_rate, f_coll_part_denom, f_coll_part_numer

  LOGICAL::chosen_i

  REAL(DP),EXTERNAL::theta
  REAL(DP),EXTERNAL::K


  ! Get seeds for the random number generator
  OPEN(UNIT=5, FILE="rand_seeds", STATUS="old")
  READ(UNIT=5,FMT=*) n0, n1, n2
  CLOSE(UNIT=5)
  
  ! Open file for data output
  OPEN(UNIT=9, FILE=fp_file, STATUS="unknown")
  OPEN(UNIT=10, FILE=mass_dens_file, STATUS="unknown")
 
  ! Initialise the distribution of field particles (monodispersed)
  DO i=1,N
     z(i)=one
  ENDDO
  
  ! Initialise total mass densities
  Phi=SUM(z)/N
  
  M=one

  t=zero
  t_stop=200d0
  
  tau_N=1!N
  
!  x=1e-6
  x=rnd240(n0,n1,n2)
  

  events=0
  
  DO WHILE (t<t_stop)

     ! Compute rates (could shove this in a subroutine)
     ! Compute collision rate
     coll_rate=0.0
     z_min=1.0e8
     DO i=1, N
        z_min=MIN(z_min, z(i))
        coll_rate=coll_rate+Phi/N*K(x,z(i))/z(i)
     ENDDO

     rate(1)=coll_rate
     rate(2)=1/theta(x)
     rate(3)=Q_in/M
     rate(4)=M*N/(Phi*tau_N)
     rate(5)=N/theta(x)
     rate(6)=N/tau_N

     rate_total=SUM(rate)

     ! Compute exponentially distributed waiting time
     tau=-LOG(rnd240(n0,n1,n2))/(rate_total)

!     WRITE(*,*) M,N,Phi,tau_N
!     WRITE(*,'(6F16.8)') (rate(i),i=1,6)
!     WRITE(*,*) rate_total, rnd240(n0,n1,n2), tau
!     PAUSE


     DO i=1,6
        p(i)=rate(i)/rate_total
     ENDDO

!     WRITE(*,'(6F16.8)') (p(i),i=1,6)


     ! Sample from this discrete distribution
!!$     DO i=1,6
!!$        u1=rnd240(n0,n1,n2)
!!$        total=0
!!$        index=1
!!$        DO j=1,6
!!$           total=total+p(j)
!!$           IF(u1<total) THEN
!!$              index=j
!!$              BREAK
!!$           ENDIF
!!$        ENDDO


     chosen_i=.FALSE.
     
     DO WHILE(chosen_i.EQV..FALSE.)
        event=FLOOR(6*rnd240(n0,n1,n2))+1
        IF(rnd240(n0,n1,n2)<p(event)) THEN
           chosen_i=.TRUE.
        ENDIF
     ENDDO

!     WRITE(*,*) event


     SELECT CASE (event)
       
     CASE(1)

           !Coagulation
           
           chosen_i=.FALSE.
           
           ! Denominator of collision partner distribution
           f_coll_part_denom=c1*Phi*K(x,z_min)/z_min
           
           ! Choose a particle 
           DO WHILE(chosen_i.EQV..FALSE.)
              !i=1.0+N*rnd240(n0,n1,n2)
              i=FLOOR(N*rnd240(n0,n1,n2))+1
              ! Numerator of collision partner distribution
              f_coll_part_numer=c1*Phi*K(x,z(i))/z(i)
              IF (rnd240(n0,n1,n2)<f_coll_part_numer/f_coll_part_denom) chosen_i=.TRUE.
           ENDDO
           
           x=x+z(i)
           
           events(1)=events(1)+1
           
        CASE(2)
           
           ! TP exchanged with FP
           i=FLOOR(N*rnd240(n0,n1,n2))+1
           x=z(i)
           events(2)=events(2)+1
           
        CASE(3)

           ! TP inflow
           !x=FLOOR(x_max*rnd240(n0,n1,n2))+1
           x=rnd240(n0,n1,n2)
           events(3)=events(3)+1
           
        CASE(4)
           
           ! FP replaced by TP
           i=FLOOR(N*rnd240(n0,n1,n2))+1
           z(i)=x
           events(4)=events(4)+1
           
        CASE(5)
           
           ! Update TP mass density
           M=(one-one/N)*M+Q_in*theta(x)/N
           events(5)=events(5)+1
           
        CASE(6)
           
           ! Update FP mass density
           Phi=(one-one/N)*Phi+M/N
           events(6)=events(6)+1
           
        END SELECT
        
     
!!$     IF (rnd240(n0,n1,n2)<(rate(3)/rate_total)) THEN
!!$        ! TP inflow
!!$        x=FLOOR(x_max*rnd240(n0,n1,n2))+1
!!$        events(3)=events(3)+1
!!$        
!!$     ELSEIF(rnd240(n0,n1,n2)<(rate(2)/rate_total)) THEN
!!$        ! TP exchanged with FP
!!$        i=FLOOR(N*rnd240(n0,n1,n2))+1
!!$        x=z(i)
!!$        events(2)=events(2)+1
!!$     ELSEIF(rnd240(n0,n1,n2)<(rate(1)/rate_total)) THEN
!!$        !Coagulation
!!$ 
!!$        chosen_i=.FALSE.
!!$
!!$        ! Denominator of collision partner distribution
!!$        f_coll_part_denom=c1*Phi*K(x,z_min)/z_min
!!$
!!$        ! Choose a particle 
!!$        DO WHILE(chosen_i.EQV..FALSE.)
!!$           !i=1.0+N*rnd240(n0,n1,n2)
!!$           i=FLOOR(N*rnd240(n0,n1,n2))+1
!!$           ! Numerator of collision partner distribution
!!$           f_coll_part_numer=c1*Phi*K(x,z(i))/z(i)
!!$           IF (rnd240(n0,n1,n2)<f_coll_part_numer/f_coll_part_denom) chosen_i=.TRUE.
!!$        ENDDO
!!$
!!$        x=x+z(i)
!!$
!!$        events(1)=events(1)+1
!!$        
!!$     ELSEIF(rnd240(n0,n1,n2)<(rate(4)/rate_total)) THEN
!!$        ! FP replaced by TP
!!$        i=FLOOR(N*rnd240(n0,n1,n2))+1
!!$        z(i)=x
!!$        events(4)=events(4)+1
!!$         
!!$     ELSEIF(rnd240(n0,n1,n2)<(rate(5)/rate_total)) THEN
!!$        ! Update TP mass density
!!$        M=(one-one/N)*M+Q_in/N*theta(x)
!!$        events(5)=events(5)+1
!!$     
!!$     ELSE
!!$        ! Update FP mass density
!!$        M=(one-one/N)*Phi+M/N
!!$        events(5)=events(5)+1
!!$        
!!$     ENDIF

     t=t+tau

     WRITE(*,'(5F16.8,2I10)') tau, t, M, Phi, x, SUM(events),event
     WRITE(10,*) tau, t, SUM(events), M, Phi, x

  ENDDO

  ! Summarise events
  WRITE(*,*)
  WRITE(*,'(6I15)') (events(i),i=1,6)
  

  DO i=1,N
     WRITE(9,*) z(i)
  ENDDO




  CLOSE(UNIT=10)
  CLOSE(UNIT=9)
  
  OPEN(UNIT=5, FILE="rand_seeds", STATUS="unknown")
  WRITE(UNIT=5,FMT=*) n0,n1,n2
  CLOSE(UNIT=5)
  
  ! Print how long the algorithm has taken
  !	xt2=dtime(ta)
  !	PRINT *, 'Time elapsed: ', xt2
  
  
END PROGRAM spm_new_main
