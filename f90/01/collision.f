      subroutine collision(x,w,number)

	INCLUDE 'param.h'

	real x(np),x1(np),f1

	integer*4 n0, n1, n2

	common /random/ n0,n1,n2

      if (ijk .lt. 2) then
	f=1.0
	do i=1,np
	x(i)=1.0
	x1(i)=1.0
	enddo
	ijk=5
	endif

	a=1.0/3.0
	c1=1.0
	c2=0.1
	w=1.0
	n=0
      

	t=0.0
      y=rnd240(n0,n1,n2)
      f1=0.0


10     continue
      tcoll=0.0
	xmin=1.0e8
	do i=1,np
	xmin=min(xmin,x(i))
	tcoll=tcoll+c1*f*(y**a)*(x(i)**(a-1.0))/np
	enddo
	tout=c2*y
      fcoll=c1*f*(y**a)*(xmin**(a-1.0))
	dt=-log(rnd240(n0,n1,n2))/(tcoll+tout)
	t=t+dt


c      REPLACEMENT
      f1=f1+w*dt
	if (rnd240(n0,n1,n2) .lt. eps2*w*dt/f) then
	i=1.0+np*rnd240(n0,n1,n2)
      x1(i)=y
	endif

c      END OF REPLACEMENT

	
      IF (rnd240(n0,n1,n2) .lt. tcoll/(tcoll+tout)) THEN

1     continue
      j=1.0+np*rnd240(n0,n1,n2)
      fcoll_1=c1*f*(y**a)*(x(j)**(a-1.0))
      if (rnd240(n0,n1,n2) .gt. fcoll_1/fcoll) goto 1
	y=y+x(j)
	n=n+1
      
      goto 10
	ENDIF

      if (number .gt. 0) write(14,*) y,f,n,t	

	print *, 1.0*n,f

      f=eps1*f1+(1.0-eps1)*f
	do i=1,np
	x(i)=x1(i)
	enddo

      END
