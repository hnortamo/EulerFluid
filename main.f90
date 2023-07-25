MODULE UTIL
        IMPLICIT NONE
contains
SUBROUTINE PRINT_FIELD(F,FN)
        REAL, DIMENSION(:,:),intent(inout) :: F
        INTEGER :: I,J
        CHARACTER(*),intent(in) :: FN
        OPEN(1, file=FN,status='replace')
        DO I=1,SIZE( F, 2 )
                DO J=1,SIZE(F,1)
                       WRITE(1,*) I-1,J-1, F(J,I)
                END DO 
                WRITE(1,*)
        END DO 
        CLOSE(1)
END SUBROUTINE
      SUBROUTINE WRITEFIELD(F,FN)
        REAL, DIMENSION(:,:),intent(inout) :: F
        INTEGER :: I,J
        CHARACTER(*),intent(in) :: FN
        OPEN(1, file=FN,status='replace')
            DO J=1,SIZE(F,1)
                WRITE(1,*) (F(J,I),I=1,SIZE(F,2))
            END DO
      END
      SUBROUTINE WRITEFIELD_I(F,FN)
        INTEGER, DIMENSION(:,:),intent(inout) :: F
        INTEGER :: I,J
        CHARACTER(*),intent(in) :: FN
        OPEN(1, file=FN,status='replace')
            DO J=1,SIZE(F,1)
                WRITE(1,*) (F(J,I),I=1,SIZE(F,2))
            END DO
      END

SUBROUTINE ADD_OBSTACLE_CIRCLE(STATE,x0,y0,r,h)
        INTEGER, DIMENSION(:,:),intent(inout) :: STATE
        REAL :: x0,y0,r,h,dx,dy
        INTEGER :: I,J
        DO I=1,SIZE( STATE,2 )
                DO J=1,SIZE(STATE,1)
                        dx = h*(I+0.5)-x0
                        dy = h*(J+0.5) -y0
                        IF ( dx*dx +dy*dy  .LT. r*r ) THEN
                             STATE(J,I) = 0.0
                        ENDIF

                END DO
        END DO
        
END SUBROUTINE

SUBROUTINE LEFT_VELOCITY_INLET(VX,STATE,INITIAL_VELOCITY)
        IMPLICIT NONE
        REAL, DIMENSION(:,:),intent(inout) :: VX
        INTEGER, DIMENSION(:,:),intent(inout) :: STATE
        REAL,intent(in) :: INITIAL_VELOCITY
        INTEGER :: I,J
        INTEGER :: s

         
        DO I=1,SIZE( VX,2 )
                DO J=1,SIZE(VX,1)
                        IF ( I .EQ. 2     ) THEN
                               VX(J,I) = INITIAL_VELOCITY  
                        ENDIF
                        s = 1.0
                        IF ( I .EQ. 1 .OR. J .EQ. 1 .OR. J .EQ. SIZE(VX,1) ) THEN
                                s = 0
                        ENDIF
                        STATE(J,I) = s 

                END DO
        END DO
END SUBROUTINE 

FUNCTION AVGVX(VX,I,J)
          INTEGER,intent(in) :: I,J
          REAL,DIMENSION(:,:),intent(in) :: VX 
          REAL :: AVGVX
          AVGVX = ( VX(J-1,I) +VX(J,I) + VX(J-1,I+1) + VX(J,I+1) ) *0.25
END FUNCTION
FUNCTION AVGVY(VY,I,J)
          REAL :: AVGVY
          INTEGER,intent(in) :: I,J
          REAL,DIMENSION(:,:),intent(in) :: VY 
          AVGVY = (VY(J,I-1) +VY(J,I) + VY(J+1,I-1) + VY(J+1,I) ) *0.25
END FUNCTION

SUBROUTINE EXTRAPOLATE(VX,VY)
        REAL,DIMENSION(:,:),intent(inout) :: VX,VY 
        INTEGER I,J
        DO I=1,SIZE( VX,2 )
               VX(1,I) = VX(2,I)
               VX(SIZE(VX,1),I) = VX(SIZE(VX,1)-1,I)
        END DO
        DO J=1,SIZE(VX,1)
               VY(J,1) = VY(J,2)
               VY(J,SIZE(VX,2)) = VY(J,SIZE(VX,2)-1)
        END DO 
END SUBROUTINE 

FUNCTION SAMPLE_FIELD(F,x,y,h,typ)
      character(*),intent(in) :: typ
      REAL :: SAMPLE_FIELD
      REAL,DIMENSION(:,:),intent(in) :: F 
      REAL,intent(in) :: h
      REAL,intent(in) :: x,y
      REAL :: xp,yp
      REAL ::dx,dy,h1,h2  
      REAL :: tx,sx,sy,ty
      INTEGER :: x0,x1,y0,y1
      h1 = 1.0/h
      h2 = 0.5*h
      !write(*,*) x,y
      xp = MAX(MIN(x,SIZE(F,2)*h),h)
      yp = MAX(MIN(y,SIZE(F,1)*h),h)
      dx = 0
      dy = 0
      SELECT CASE (typ)
         CASE("VX")
                 dy = h2
         CASE("VY")
                 dx = h2
         CASE("Smoke")        
                 dx = h2
                 dy = h2
      END SELECT
      x0 = MIN(FLOOR((xp-dx)*h1)+1,size(F,2))
      y0 = MIN(FLOOR((yp-dy)*h1)+1,size(F,1))
      tx = ((xp-dx) - (x0-1)*h) * h1
      ty =  ((yp-dy) - (y0-1)*h) * h1
      x1 = MIN(x0+1,size(F,2))
      y1 = MIN(y0+1,size(F,1))
      sx = 1- tx
      sy = 1 - ty
      !write(*,*) x0,y0,x,y 
      SAMPLE_FIELD=sx*sy*F(y0,x0) + tx*sy*F(Y0,X1) + tx*ty*F(Y1,X1) + sx*ty*F(Y1,X0)

END FUNCTION
  subroutine swap_fields(curr, prev)


    REAL,dimension(:,:),allocatable, intent(inout) :: curr, prev
    real, allocatable, dimension(:,:) :: tmp

    call move_alloc(curr, tmp)
    call move_alloc(prev, curr)
    call move_alloc(tmp, prev)
  end subroutine swap_fields

SUBROUTINE ADVECTSMOKE(VX,VY,SMOKE,SMOKE_T,S,dt,h)
        REAL,allocatable,dimension(:,:),intent(in) :: VX,VY
        REAL,dimension(:,:),allocatable,intent(inout) :: SMOKE,SMOKE_T 
        INTEGER,dimension(:,:),intent(in) ::S
        INTEGER I,J
        REAL :: dt,h
        REAL :: x,y,h2,u,v
        h2 = h*0.5
        DO I=2,SIZE( VX,2 )-1
                DO J=2,SIZE(VX,1)-1
                        IF ( S(J,I) .NE. 0 )  THEN
                                u = (VX(J,I) + VX(J,I+1) )*0.5
                                v = (VY(J,I) + VY(J+1,I) )*0.5 
                                x = (i-1)*h + h2 - dt*u
                                y = (j-1)*h + h2 - dt*v
                                SMOKE_T(J,I) = SAMPLE_FIELD(SMOKE,x,y,h,"Smoke")
                        ENDIF  
                END DO
        END DO

        call swap_fields(SMOKE,SMOKE_T)
END SUBROUTINE

SUBROUTINE ADVECTVEL(VX,VY,VX_T,VY_T,S,dt,h)
        REAL,allocatable,dimension(:,:) :: VX,VX_T,VY,VY_T
        INTEGER,allocatable,dimension(:,:) ::S
        INTEGER I,J
        REAL :: dt,h
        REAL :: x,y,h2,u,v
        h2 = h*0.5
        DO  I=2,SIZE( VX,2 )-1
                DO J=2,SIZE(VX,1)-1
                        IF (S(J,I) .NE. 0 .AND. S(J,I-1) .NE. 0 .AND. J .LT. SIZE(VX,1) ) THEN
                                x = (I-1)*h
                                y = (J-1)*h+h2
                                u = VX(J,I)
                                v = AVGVY(VY,I,J)
                                x = x -dt*u
                                y = y -dt*v
                                u = SAMPLE_FIELD(VX,x,y,h,"VX")
                                VX_T(J,I) =  u
                        ENDIF 
                        IF (S(J,I) .NE. 0 .AND. S(J-1,I) .NE. 0 .AND. I .LT. SIZE(VX,2) ) THEN
                                x = (I-1)*h +h2
                                y = (J-1)*h
                                u = AVGVX(VX,I,J)
                                v = VY(J,I)
                                x = x -dt*u
                                y = y -dt*v
                                v = SAMPLE_FIELD(VY,x,y,h,"VY")
              
                                VY_T(J,I) =  v
                        ENDIF 
                END DO
        END DO
        call swap_fields(VX,VX_T)
        call swap_fields(VY,VY_T)

END SUBROUTINE

SUBROUTINE solveIncompressibility(VX,VY,State,niters,relax,density,dt,h)
        REAL,dimension(:,:) :: VX,VY
        INTEGER,dimension(:,:) :: STATE
        REAL :: relax,dt,h,density
        REAL :: cp,sx0,sx1,sy0,sy1,p,div
        INTEGER :: I,J,K,niters,s
        cp = density *h /dt 

        DO K=1,niters
        DO I=2,SIZE( VX,2 )-1
                DO J=2,SIZE(VX,1)-1
                        IF (STATE(J,I) .EQ. 0) CYCLE 
                        s = STATE(J,I)
                        sx0 = STATE(J,I-1)
                        sx1 = STATE(J,I+1)
                        sy0 = STATE(J-1,I)
                        sy1 = STATE(J+1,I)
                        s = sx0 + sx1 + sy0 + sy1
                        IF ( s .EQ. 0 ) CYCLE
                        div = VX(J,I+1) - VX(J,I) + VY(J+1,I) - VY(J,I)
                        p = (-div /s)*relax
                        VX(J,I)=VX(J,I)-sx0*p
                        vx(j,i+1)=vx(j,i+1)+sx1*p
                        vy(j,i)  =vy(j,i)-sy0*p
                        vy(j+1,i)=vy(j+1,i)+sy1*p
                END DO
        END DO
        END DO
END SUBROUTINE

END MODULE

PROGRAM FLUIDSIMULATION
        USE util
        IMPLICIT NONE
        REAL :: density,dx,dt
        REAL :: velocity,overRelaxation 
        INTEGER :: NX,NY,res,steps,niters,IT
        REAL :: simulationHeight,simulationWidth 
        REAL,DIMENSION(:,:), allocatable :: VX,VY,pressure,smoke
        INTEGER,DIMENSION(:,:), allocatable :: state
        real :: start, finish
        REAL,DIMENSION(:,:), allocatable :: VX_T,VY_T,smoke_t 
        CHARACTER(100) :: FN



        density=1000
        velocity=2


        steps=1000
        res=200
        overRelaxation = 1.9
        niters=100
        simulationHeight=1.0
        simulationWidth=3.0
        dt = 1.0 / 60.0

        dx = simulationHeight / res
        NX = FLOOR(simulationWidth / dx) 
        NY = FLOOR(simulationHeight / dx )
        allocate( VX(NY,NX),VY(NY,NX),state(NY,NX),VX_T(NY,NX),VY_T(NY,NX),smoke(NY,NX),smoke_t(NY,NX))
        

        write(*,*) "Starting simulation"
        write(*,*) "Dimensions: ", NX , "x" , NY
        write(*,*) "RES and dt: ",dx,dt,dt*velocity
        CALL LEFT_VELOCITY_INLET(VX,state,velocity)
        CALL LEFT_VELOCITY_INLET(VX_T,state,velocity)
        CALL ADD_OBSTACLE_CIRCLE(state,1.5,0.5,0.15,dx) 
!        CALL WRITEFIELD(VX,"VelX_.out")
!        CALL WRITEFIELD_I(STATE,"state.out")
        SMOKE=0
        SMOKE_T=0
        DO IT=40,60
                SMOKE(IT,1) = 4
                SMOKE_T(IT,1) = 4
        END DO
        DO IT=1,steps
                WRITE(*,*) "STEP: ", IT
          WRITE(FN,'(A5,I4.4,A4)') 'VelX_',IT,'.out'
          CALL WRITEFIELD(VX,FN) 
          WRITE(FN,'(A6,I4.4,A4)') 'Smoke_',IT,'.out'
          CALL WRITEFIELD(SMOKE,FN) 
          WRITE(FN,'(A5,I4.4,A4)') 'VelY_',IT,'.out'
          CALL WRITEFIELD(VY,FN) 
          CALL solveIncompressibility(VX,VY,State,niters,overRelaxation,density,dt,dx)
          CALL EXTRAPOLATE(VX,VY)
          VX_T = VX
          VY_T = VY
          CALL ADVECTSMOKE(VX,VY,smoke,smoke_t,state,dt,dx)
          CALL ADVECTVEL(VX,VY,VX_T,VY_T,State,dt,dx)
        END DO


        deallocate(VX,VY,state,VX_T,VY_T,smoke,smoke_t)
END PROGRAM FLUIDSIMULATION
