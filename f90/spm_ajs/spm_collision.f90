SUBROUTINE spm_collision(x,w,write_to_file)

	USE global_data
	IMPLICIT NONE

	INTEGER, INTENT(IN)::write_to_file
	REAL(KIND=8), EXTERNAL::rnd240

	REAL::x(np),x1(np),f1, w
	REAL(KIND=8)::y
	INTEGER:: i, j, n
	REAL::t1, t2, tout, tin, xmin, xmax
	REAL::tcoll, fcoll, fcoll_1, tau, dt, t

	w=1.0 !Q_in
	n=0
      

	t=0.0
	y=rnd240(n0,n1,n2)
	f1=0.0


10	CONTINUE
	tcoll=0.0
	xmin=1.0e8


	DO i=1,np
		xmin=MIN(xmin,x(i))
		tcoll=tcoll+c1*f*(y**a)*(x(i)**(a-1.0))/np
	ENDDO

	tout=c2*y
	fcoll=c1*f*(y**a)*(xmin**(a-1.0))
	dt=-LOG(rnd240(n0,n1,n2))/(tcoll+tout)
	t=t+dt


	! Start of particle replacement 
	f1=f1+w*dt
	IF (rnd240(n0,n1,n2)<eps2*w*dt/f) THEN
		i=1.0+np*rnd240(n0,n1,n2)
		x1(i)=y
	ENDIF
	! End of particle replacement

	
	IF (rnd240(n0,n1,n2)<tcoll/(tcoll+tout)) THEN

1		CONTINUE
		j=1.0+np*rnd240(n0,n1,n2)
		fcoll_1=c1*f*(y**a)*(x(j)**(a-1.0))
		IF (rnd240(n0,n1,n2)>fcoll_1/fcoll) GOTO 1
		y=y+x(j)
		n=n+1
      
		GOTO 10
	ENDIF

	IF (write_to_file>0) WRITE(14,*) y,f,n,t	

	PRINT *, 1.0*n,f

	f=eps1*f1+(1.0-eps1)*f
	DO i=1,np
		x(i)=x1(i)
	ENDDO

END SUBROUTINE spm_collision
