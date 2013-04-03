PROGRAM binarywriter
IMPLICIT NONE
INTEGER X(10),Y(10), a
INTEGER :: i,j
REAL*8 :: rnd !random real
REAL*8, DIMENSION(10,15) ::  topog

DO a = 1, 10
X(a)=a
Y(a)=a+100
END DO

DO i=1,10
DO j=1,15

CALL random_number(rnd)
topog(i,j)  =0.2d0*dble(i)+0.2d0*dble(j) + 1000.d0

END DO
END DO

OPEN(10,status='unknown',file='bin.dat',form='UNFORMATTED')
write(10) X,Y
CLOSE(10)

OPEN(11, status='unknown', file='bin2.dat', form='Unformatted')
write(11) topog
CLOSE(11)

OPEN(12, file='num.dat', delim='apostrophe')
write(12,'(e13.6)') topog
CLOSE(12)
END PROGRAM binarywriter

