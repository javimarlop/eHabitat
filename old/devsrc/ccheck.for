
      LOGICAL FUNCTION OUTREG (npoly,X,Y)
C        CHECK WHETHER POINT (X,Y) IS WITHIN REGION POLY
C        GIVEN BY ITS INDICES 
C        .TRUE. IF POINT (X,Y) OUT  OF REGION 
C        .TRUE. IF POINT (X,Y) ON BOUNDARY 
C        THIS IS NOT EXACTLY THE SAME AS COELEC

      INTEGER MAXPKT,MAXPOL
      PARAMETER (maxpkt=20000,maxpol=500)
      REAL XP,YP,EPSOR,EPSOR2,EPSOR3  
      INTEGER IXMIN,IXMAX,IYMIN,IYMAX   
      INTEGER I,I1,I2,ICOUNT,NINTS
      REAL XMIN,XMAX,YMIN,YMAX
      REAL DIST,S1,S2,X1,X2,Y1,Y2,SS,XOUT,YOUT,X,Y
      COMMON /PKTLI/ XP(MAXPKT),YP(MAXPKT)
      COMMON /PARAM/ XMIN,XMAX,YMIN,YMAX,IXMIN,IXMAX,IYMIN,IYMAX
                          
C        IMPROVE DEFFINING OF EPSILON

      EPSOR=0.00001
      EPSOR2=0.0000001
      EPSOR3=0.01  

      call GETMAX(NPOLY)

C        SELECT POINT WHICH IS OUT OF REGION BY SURE  
      XOUT=ABS(XMAX)
      IF (ABS(XMIN).GT.XOUT) XOUT=ABS(XMIN)  
      YOUT=ABS(YMAX)
      IF (ABS(YMIN).GT.YOUT) YOUT=ABS(YMIN)
      XOUT=2.1*XOUT
      YOUT=2.2*YOUT
      ICOUNT=0
444   CONTINUE
C        NUMBER OF INTERSECTION POINTS
      NINTS=0   
C        TURN LINE
      YOUT=YOUT*0.999  
      ICOUNT=ICOUNT+1   
      IF (ICOUNT.EQ.10) THEN
C         PRINT*,'ICOUNT=10 - NO INTERSECTION'
C           (ASSUMED TO BE ON BOUNDARY)
         GOTO 666     
      END IF   
C        INTERSECT SEGMENT (X,Y)-(XOUT,YOUT) WITH POLY
      DO 100 I=1,NPOLY
C          NODES OF THE SEGMENT TO BE ADDRESSED
         I1=i
         IF (I.EQ.NPOLY) THEN 
            I2=1
         ELSE
            I2=i+1  
         END IF
C           LOOK POINTS OF POLY UP IN XP LIST (MUST BE EXISTING POINTS)
         X1=XP(I1)
         Y1=YP(I1)
         X2=XP(I2)
         Y2=YP(I2)                    

         SS=(XOUT-X)*(Y2-Y1)-(X2-X1)*(YOUT-Y)         
         IF (ABS(SS).LE.EPSOR3) THEN
C              PARALLEL      
C              CHECK IF IDENTICAL OR NOT  
            DIST=(Y1*(X2-X1)-X1*(Y2-Y1))-
     +                 (YOUT*(XOUT-X)-XOUT*(YOUT-Y))
            IF (ABS(DIST) .LT. EPSOR3) THEN
C              LINES ARE IDENTICAL  
C                 TRY AGAIN
               GOTO 444
            ELSE  
C              THIS IS A VALID NO CUT 
            END IF         
         ELSE     
C              INTERSECTION IN   I1,I2
            S1=((X1-X)*(YOUT-Y)-(Y1-Y)*(XOUT-X))/SS  
C              INTERSECTION IN (X,Y)-(XOUT,YOUT)
            S2=((X1-X)*(Y2-Y1)-(Y1-Y)*(X2-X1))/SS                   
C              POINT (X,Y) IS RIGHT ON A VERTEX
            IF ((ABS(S1).LE.EPSOR2.OR.ABS(S1-1.).LE.EPSOR2)
     +       .AND.(ABS(S2).LE.EPSOR2.OR.ABS(S2-1.).LE.EPSOR2))THEN                 
C                 MAY NOT BE ON ANY VERTEX -
C                 ASSUMED TO BE ON BOUNDARY
               GOTO 666
            END IF                             
            IF (ABS(S1).LE.EPSOR2.OR.ABS(S1-1.).LE.EPSOR2) THEN
C                 LINE  (X,Y)-(XOUT,YOUT) GOES RIGHT THROUGH VERTEX                      
C                 TRY AGAIN
               GOTO 444   
            ELSE
               IF ((ABS(S2).LE.EPSOR2).AND. 
     +               (S1.GT.0.AND.S1.LT.1.)) THEN  
C                    POINT (X,Y) IS RIGHT ON POLYGON, GIVE RESULTS
                  GOTO 666  
               END IF                 
C                 NOTE: S2=1 MEANS XYOUT IS ON POLYGON (MAY NOT BE) 
               IF ((S1.GT.0.AND.S1.LT.1.).AND.
     +              (S2.GT.0.AND.S2.LT.1.)) THEN
C                     THIS IS A REGULAR CUT 
                  NINTS=NINTS+1
               ELSE 
C                     THIS IS A VALID NO CUT
               END IF
            END IF
         END IF
99     continue
100   CONTINUE   
C        COUNT INTERSECTION POINTS :EVEN  
C        .TRUE. IF POINT (X,Y) OUT  OF REGION
      IF (MOD(NINTS,2).EQ.0) THEN 
         OUTREG=.TRUE.
      ELSE
         OUTREG=.FALSE.
      END IF   
      RETURN
666   CONTINUE 
C        POINT IS ON THE BOUNDARY
      OUTREG=.TRUE.
      RETURN
      END FUNCTION OUTREG




      SUBROUTINE GETMAX(NPOLY)
C
CC    SUBROUTINE GETMAX BERECHNET DIE
CC    MIN. UND MAX. X- UND Y-KOORDINATEN EINES POLYGONS 
C     .. WITH INDICES OF THE VERTICES GIVEN
C
      INTEGER MAXPKT
      PARAMETER (maxpkt=20000)
      INTEGER NPOLY
      INTEGER IXMIN,IXMAX,IYMIN,IYMAX
      REAL XP,YP,zp,XMIN,XMAX,YMIN,YMAX
      COMMON /PKTLI/ XP(MAXPKT),YP(MAXPKT),ZP(MAXPKT)
      COMMON /PARAM/ XMIN,XMAX,YMIN,YMAX,IXMIN,IXMAX,IYMIN,IYMAX          
C
      ip=1
      XMIN=XP(ip)
      XMAX=XP(ip)
      IXMIN=ip
      IXMAX=ip
      YMIN=YP(ip)
      YMAX=YP(ip)
      IYMIN=ip
      IYMAX=ip
      DO 1000 I=2,NPOLY
         IP=i
         IF(XP(IP).LT.XMIN)THEN
            XMIN=XP(IP)
            IXMIN=IP
         ENDIF
         IF(XP(IP).GT.XMAX)THEN
            XMAX=XP(IP)
            IXMAX=IP
         ENDIF
         IF(YP(IP).LT.YMIN)THEN
            YMIN=YP(IP)
            IYMIN=IP
         ENDIF
         IF(YP(IP).GT.YMAX)THEN
            YMAX=YP(IP)
            IYMAX=IP
         ENDIF
 1000 CONTINUE
      RETURN
      END SUBROUTINE GETMAX      
