PROGRAM Sensitivity

!use csv_file

IMPLICIT NONE

REAL,PARAMETER  :: pi = 3.14159
INTEGER :: m, n, mn, infOrder, vegOrder, sEmerge, tPersist, tsteps, nSteps, precip
INTEGER :: i, j, k, l, esteps 
INTEGER, ALLOCATABLE ::clock(:)
INTEGER :: infOrdermin,infOrdermax,vegOrdermin,vegOrdermax 
REAL*8 :: roughness
REAL*8 :: kv, kb, Dv, Db, rnd
REAL*8 :: K0, Kmax, kf, rf, alpha
REAL*8 :: Emax, kc,rc,beta
REAL*8 :: gamma, Esb, Esv, bav, Psv, Psb
REAL*8 :: dx, dy  !spatial dimesions of lattice cells
REAL*8 :: pa, pbar, ie !pa = mean annual precip, pbar = water in each particle, ie = mean rain intensity
REAL *8 :: ts, te  !ts = duration raining, te = time step for evap calculations
INTEGER :: np, ne !number of particles, ne = number of evap iterations
INTEGER :: vegmax, isSEmerge, IOStatus
REAL*8 :: pc
INTEGER, DIMENSION(6) :: simflags
REAL*8, DIMENSION(4) :: climParams
Real*8, DIMENSION(6) :: infiltParams
REAL*8, DIMENSION(8) :: evapParams
REAL*8, DIMENSION(5) ::  vegParams
REAL*8, DIMENSION(4) :: erParams
REAL*8,DIMENSION(14) :: readRealParams
Integer, DIMENSION(9) :: readIntParams

CHARACTER(LEN=8) :: fid, fid2

!check how many parameter sets are in params.txt
!Check params file for length
PRINT *,'Reading file "params.txt" ...'
Open(1,FILE="params.txt", status='old') !parameter settings are in params.txt
!first three rows are headers
READ(1,*)
READ(1,*)
READ(1,*)

!the rest are parameter sets
k=0	!count for the rows
IOStatus=0
DO WHILE (IOStatus.eq.0)
READ(1,*,IOSTAT=IOStatus)
k=k+1
END DO
CLOSE(1)

k=k-1
k=k/2
PRINT *,k," parameter sets"

!Read input parameter file
!put program in loop to run through parameter file
Open(99,FILE="params.txt", status='old')
READ(99,*)
READ(99,*)
READ(99,*)

DO i=1,k  

  READ(99,'(9i9)') readIntParams(:)
  READ(99,'(14f23.7)') readRealParams(:)

   PRINT *,readIntParams
   PRINT *,readRealParams
  m=readIntParams(1)  !numbers of rows and columns
  n=readIntParams(2)
  mn = m*n
  np = readIntParams(3) !number of particles of rain falling
  !PRINT *,"m,n",m,n
  dx = readRealParams(1)  !dimensions of cells
  dy = readRealParams(1)  !
  pa = readRealParams(2) !400.0    !annual precip
  !PRINT *,"dx,dy,pa",dx,dy,pa
  ts = readRealParams(3) !12.0/365.0 !duration in years of time raining
  pbar = pa/dble(np) !mm 
  ie = pa / ts  !effective rainfall intensity mm / year
  !PRINT *,"ts,pbar,ie",ts,pbar,ie
  K0 = readRealParams(4) * (365.d0* 24.d0)  !from mm/hr to mm/year
  Kmax = readRealParams(5) * (365.d0* 24.d0)  !from mm/hr to mm/year
  kf = readRealParams(6) !0.9d0 !per meter : rate of decline in plant effect on infiltration
  rf = readRealParams(7) !10.d0  ! meter: maximum length for plant effect on infiltration
  alpha = (Kmax-K0)/K0
  !PRINT *,"K0,Kmax,kf,rf,alpha",K0,Kmax,kf,rf,alpha
 !K(r) = K0 +K0 *alpha*exp(-kf*r)
 !K[r_] := K0 + If[r < rf, K0*alpha*Exp[-kf*r], 0]
 !PInf[r_] := Min[1, K[r]/ie]

  Emax = readRealParams(8) !0.5 * (365.0* 24.0)  !from mm/hr to mm/year
  kc = readRealParams(9) !1.8d0 !per meter : rate of decline in plant water uptake with distance
  rc = readRealParams(10) !10.0 !meter: maximum length for plant water uptake
  beta = Emax
 !PRINT*,"Emax,kc,rc,beta",Emax,kc,rc,beta
 !ET(r) =  beta*exp(-kc*r)*dx*dy
  gamma = readRealParams(11) !0.1 !relative reduction of soil evap under canopy
  bav = readRealParams(12) !0.03  !scaling factor for Esb calc
  Esb = bav*Emax  !maximum bare soil evap rate
  Esv = gamma*Esb   !soil evap from under canopy
  !PRINT*,"gamma,bav,Esb,Esv",gamma,bav,Esb,Esv
  ne = Int(Ceiling(Max((1.0-ts)*Esb/pbar,(1.0-ts)*dx*dy*beta/pbar)))
  !ne is the number of times required to divide 1-ts by in order to ensure only one particle   removed
  !by one evap process at any one time
  te = (1-ts)/ne  !years: length of a time step in evap calcs

  Psb = te * Esb / pbar
  Psv = te * Esv / pbar
  !Pet(r) = te * dx * dy*beta*exp(- kc*r)/pbar

  nSteps = readIntParams(4) !100  !total number of "years"
  !esteps = 1
  tSteps = readIntParams(5) !200 !ne !number of iterations for evap calcs between veg change
  !PRINT*,"te,Psb,Psv,nSteps,tSteps",te,Psb,Psv,nSteps,tSteps

  precip = np

  sEmerge = readIntParams(6) !120
  vegmax = readIntParams(7) !9
  pc = readRealParams(13) !0.01d0
  tPersist = readIntParams(8) !200
  isSEmerge = readIntParams(9)
  roughness = readRealParams(14) !0.00005d0
  !PRINT*,"sEmerge,vegmax,pc,tPersist,isSEmerge,roughness",sEmerge,vegmax,pc,tPersist,isSEmerge,roughness
  !READ*,
 
  CALL RANDOM_SEED(size=l)
  ALLOCATE(clock(l))
  DO j=1,l
    CALL SYSTEM_CLOCK(COUNT=clock(j))
    clock(j)=clock(j)+j
  END DO
  CALL RANDOM_SEED(PUT = clock)
  DEALLOCATE(clock)

  kv = 10.0d0 !readRealParams(15)
  Dv = 1.0d0 !readRealParams(16)
  !i=3
  j=1
  WRITE(fid,'(i4)') i
  WRITE(fid2,'(i4)') j
  fid=trim(trim(adjustl(fid))//'_'//trim(adjustl(fid2)))
  !PRINT*,"fid",fid
  simflags = (/ 1,0,1,1,1, nSteps/)
  climParams = (/dble(np), pbar, ie, roughness  /)
  infiltParams = (/K0,rf,kf, Kmax, dx,dy /)
  evapParams = (/rc, kc,dble(tSteps),te,Psb,Psv,beta ,dble(ne)/)
  vegParams = (/ dble(vegmax), dble(sEmerge),dble(tPersist),pc,dble(isSEmerge) /)
  erParams = (/0.1d0, 0.1d0, 0.1d0, 0.1d0 /)
   
  CALL SimCODE(m,n,mn,simflags,climParams,infiltParams,evapParams,vegParams,erParams,fid)

END DO  !end of loop through parameter file
CLOSE(99)

CONTAINS

!SimCode
SUBROUTINE SimCODE(m,n,mn,simflags,climParams,infiltParams,evapParams,vegParams,erParams, resultsFID)
 
IMPLICIT NONE

INTEGER, INTENT(IN) :: m !# of rows
INTEGER, INTENT(IN) :: n !# of columns
INTEGER,INTENT(IN) :: mn !m*n
INTEGER, DIMENSION(6), INTENT(IN) :: simflags
REAL*8, DIMENSION(4), INTENT(IN) :: climParams  !climate parameters
REAL*8, DIMENSION(6), INTENT(IN) :: infiltParams !veg and soil kernel parameters 
REAL*8, DIMENSION(8), INTENT(IN) :: evapParams !evaporation parameters
REAL*8, DIMENSION(5), INTENT(IN) :: vegParams !veg and soil kernel parameters 
REAL*8, DIMENSION(4), INTENT(IN) :: erParams !veg and soil kernel parameters 
CHARACTER(LEN=6), INTENT(IN) :: resultsFID  !results file id code
!*************************************************************************************
CHARACTER(len=10) :: fid
CHARACTER(len=3) :: mc
CHARACTER delimiter
CHARACTER*150 command
CHARACTER(len=255) :: cwd, savedir
INTEGER :: i,j,k, m1,n1,outflow, sold
INTEGER :: nSteps, eSteps, flag, solMax, ne

REAL*8 :: K0, rfx, rfy, kf, Kmax, dx ,dy !infilt variables
REAL*8 :: rcx, rcy, kc, te, Psb, Psv, beta  !evap variables
REAL*8,DIMENSION(7) :: eparams !evap input array

REAL*8 :: rnd
REAL*8 :: roughness
REAL*8 :: kv, kb, Dv, Db

INTEGER :: sEmerge, tPersist, vegmax  !veg change variables
REAL*8 :: pc  !veg change variable

INTEGER, DIMENSION(mn,2) :: randOrder
REAL*8, DIMENSION(m,n) ::  infiltKern, storeKern, topog, flowResistance0,flowResistance1
REAL*8, DIMENSION(m,n) :: manningsN, infx  !kinematic water variables
REAL*8, DIMENSION(m,n) :: Ksat, wfs, cumInfilt,  cumInfiltOld, inflow, infex !kinematic water variables
REAL*8, DIMENSION(m,n,2) :: disOld,  disNew, iex !kinematic water variables

INTEGER, DIMENSION(m,n) :: precip, newflowdirns, store,discharge,veg, eTActual
INTEGER, DIMENSION(m,n) :: flowdirns, lakes, bareE, dummyveg, solutionOrder
INTEGER, DIMENSION(m,n) :: solOrder
INTEGER, DIMENSION(m,n,9) :: mask
REAL*8, DIMENSION(m,n) :: alpha, deltax

LOGICAL :: lexist


!simflags(1) = 1 !then use topography to route flows
!simflags(2) = 1 and simflags(1) = 1 then simulate erosion and update flow pathways
!simflags(3) = 1 then simulate evaporation
!simflags(4) = 1 and simflags(1) = 1 then simulate evolving vegetation
!simflags(5) = 1 allow veg to be randomly distributed initially, othewrwise set all veg to 0
flag = simflags(4)
nSteps = simflags(6)
!**************************************************************************************

!*************************************************************************************
!*** Model Parameters ****************************************************************
!************************************************************************************

!PRINT*,"in Sim"
!Climate /landscape parameters
precip = int(climParams(1))
pbar = climParams(2)
ie = climParams(3)
roughness = climParams(4)

!infiltration kernel parameters
K0 = infiltParams(1)
rfx = infiltParams(2)
rfy = infiltParams(2)
kf =  infiltParams(3)
Kmax = infiltParams(4)
dx = infiltParams(5)
dy = infiltParams(6)

!evaporation kernel parameters
rcx = evapParams(1)
rcy = evapParams(1)
kc = evapParams(2)
eSteps = int(evapParams(3))
te = evapParams(4)
Psb = evapParams(5)
Psv = evapParams(6)
beta = evapParams(7)
ne = int(evapParams(8))

eSteps = max(min(ne,eSteps),1)

eparams = (/dx,dy,te,pbar,Psb,Psv,beta /)

!vegetation parameters
!vegmax = int(vegParams(1)) ! parameter set to define maximum biomass was 9
!sEmerge = int(vegParams(2))
!etPersist = int(vegParams(3))
!pc = vegParams(4)
!usePCFlag = int(vegParams(5))

!erosion parameters 
kv = erParams(1)
kb = erParams(2)
Dv = erParams(3)
Db = erParams(4)


!**************************************************************************************
!**** Initial conditions ****************************************************************
!*************************
DO i=1,m
  DO j=1,n
   CALL random_number(rnd)
   flowResistance0(i,j) = (dble(rnd)+1.0d0)*kb
   CALL random_number(rnd)
   !If(i<j) THEN
   !If[i > j, 60 - ((i 1 + j 6) ), 60 - (( i 6 + j 1))]
   topog(i,j)  =0.2d0*dble(i)+0.2d0*dble(j) + 1000.d0
   !topog(i,j) = 1000.d0 +0.01d0*dble(i)+0.01d0*dble(j) + dble(rnd -0.5d0)*roughness
   !ELSE
   !  topog(i,j)  = 0.0*Sin(real(i)/5.0)*Sin(real(j)/5.0) +0.2d0*dble(i)+0.1d0*dble(j) + 1000.d0
   !ENDIF
   IF (simflags(5).gt.0) THEN
     IF (rnd>0.9) THEN
       veg(i,j) = 1
     ELSE
       veg(i,j) = 0
     END IF
   ELSE
    veg(i,j) = 0
   END IF
 END DO
END DO
store=0
lakes = 0
flowdirns = -2
 
 !manningsN = 0.02d0
 !iex = 0.0d0
 !disOld = 1e-20
 !disNew = 1e-20
 !PRINT*,"Running GD8"
 CALL GD8(topog,flowdirns, m,n)
 !PRINT*,"Finished GD8"
 !READ*,
 !CALL KWPrimer(m,n,topog, manningsN, mask,solOrder, flowdirns,alpha,deltax, solMax)
 !Ksat=10.d0  
 !wfs = 90.d0 * 0.434d0
 !cumInfilt = 0.d0
 
 !OPEN(3,file='Discharge.out')
 !DO j=1,100
 !  CALL random_number(rnd)
 !  inflow = 5.0 *rnd + disNew(:,:,2)  !precip + runon
 !  cumInfiltOld = cumInfilt
 !  CALL GAInfilt (m,n,0.1d0,inflow, Ksat, wfs, cumInfilt, infex)
 !  iex(:,:,1) = iex(:,:,2)
 !  iex(:,:,2) = 15.0 - (cumInfilt-cumInfiltOld)/1.0d0 !/dt
 !  CALL KinematicWave(m,n,2,0.1d0, iex,flowdirns,solOrder,solMax,mask,alpha,deltax,disOld,disNew)
 !DO i=1,m
 !  WRITE(3,'(32f16.6)') disNew(i,:,2)
 !END DO
 !  PRINT *,j,'----------------------'
 !END DO
 !iex = 0.0
 !DO j=1,100
 !  inflow = 0.0 + disNew(:,:,2)  !precip + runon
 !  cumInfiltOld = cumInfilt
 !  CALL GAInfilt (m,n,1.0d0,inflow, Ksat, wfs, cumInfilt, infex)
 !  iex(:,:,1) = iex(:,:,2)
 !  iex(:,:,2) = 0.d0 !- (cumInfilt-cumInfiltOld)/1.0d0 !/dt
 !  CALL KinematicWave(m,n,2,0.5d0, iex,flowdirns,solOrder,solMax,mask,alpha,deltax,disOld,disNew)
  
 
 !DO i=1,m
 !  WRITE(3,'(32f16.6)') disNew(i,:,2)
 !END DO
 !  PRINT *,j + 100 ,'----------------------'
 !END DO
 
write(fid,'(i3)') n  !internal write for colum number
!PRINT *,"simflags",simflags
IF (simflags(1).eq.1) THEN
  ! CALL GD8(topog,flowdirns, m,n)
  !CALL KMWOrder(flowdirns,m,n,solutionOrder)
  !DO i=1,m
  !  WRITE(*,'('//fid//'i3,a2,'//fid//'i3)') flowdirns(i,:),"--",solutionOrder(i,:)
  !END DO
   CALL InfiltProb(veg,m,n,K0,ie,rfx,rfx,kf,Kmax,dx,dy,infiltKern)
   !PRINT*,"out InfiltProb"
   !READ*,
   !PRINT *, infiltKern
   !READ*,
END IF
!WHERE (veg>0)
!  flowResistance1= flowResistance0 + kb
!ELSEWHERE
!  flowResistance1= flowResistance0
!END WHERE
!infiltKern =   infiltKern * slopeFactor
storeKern = 99999.d0


!inquire(file='TStepResults'//trim(ADJUSTL(resultsFID))//'.out',exist=lexist)
!IF (lexist) THEN
!  OPEN(1,file='TStepResults'//trim(ADJUSTL(resultsFID))//'.out',POSITION='APPEND')
!  OPEN(2,file='SummaryResults'//trim(adjustl(resultsFID))//'.out',POSITION='APPEND')
!ELSE
  OPEN(1,file='TStepResults'//trim(ADJUSTL(resultsFID))//'.out')
  OPEN(2,file='SummaryResults'//trim(adjustl(resultsFID))//'.out')
  write(1,*) m,n,'veg,flowdirns,store,discharge,eTActual,bareE,topog,flowResistance1'
  write(2,*) 'timeStep,vegDensity,totalET,totalBE,totalStore, totalDischarge, totalOutflow'
!END IF

!***********************************************************************************
!Timestep iterations
DO j=1,nSteps
 !IF ((j>int(nSteps/2)).and.(maxval(veg).eq.0)) THEN
 !  CLOSE(1)
 !  CLOSE(2)
 !  RETURN
 !END IF
  eTActual = 0
  bareE = 0
  !outflow = 0
  DO k=1,1  !erosion, routing, evaporation, loop
    IF (simflags(1).eq.1) THEN
      
      !PRINT*,"in Routing"
      CALL RoutingWithKernel(m,n,mn,precip, infiltKern, storeKern, flowdirns,topog, store,discharge,outflow)
      !PRINT*,"out Routing"
    END IF
   
    IF ((simflags(1).eq.1).and.(simflags(2).eq.1)) THEN
      CALL Erosion(discharge,topog,flowdirns,veg,flowResistance1,m,n)
      CALL Splash(topog,veg, Dv, Db,m,n)
      CALL GD8(topog,flowdirns, m,n)
    END IF
    IF (simflags(3).eq.1) THEN
      !PRINT*,"in Evap"
      CALL Evaporation(veg,eTActual,bareE,store,eSteps,rcx,rcy,kc,dx,dy,eparams)
       !PRINT *, "End evap"
       dummyveg = veg
       CALL VegChange(dummyveg,m,n,vegParams, store, eTActual,1,0)
       !PRINT *, "End vegChange1: Emerging Veg"
       dummyveg = dummyveg - veg  !identify just the new veg
      If (ne > eSteps) THEN
        !PRINT *, "in evap2",ne-eSteps
        CALL Evaporation(veg,eTActual,bareE,store,ne - eSteps,rcx,rcy,kc,dx,dy,eparams)
        !PRINT *, "out evap2"
      END IF
    END IF
  END DO
  
IF (MOD(j,1).eq.0) THEN
  !flag = 1
  WRITE(1,*) "time step = ",j
  DO i=1,m
    WRITE(1,'('//fid//'i3,'//fid//'i4,'//fid//'i12,'//fid//'i12,'//fid//'i12,'//fid//'i12,'//fid//'e14.6,'//fid//'e14.6)') & 
      veg(i,:),flowdirns(i,:),store(i,:),discharge(i,:),eTActual(i,:),bareE(i,:),topog(i,:),infiltKern(i,:)
  END DO
  WRITE(2,'(i3,e14.6,5i10)') j,dble(sum(veg))/dble((m*n)), sum(ETActual), sum(bareE), sum(store), sum(discharge), outflow
ELSE
  !flag = 0
END IF  


 
  flag = 1
  IF ((simflags(3).eq.1).and.(flag.eq.1)) THEN
     CALL VegChange(veg,m,n,vegParams, store, eTActual,1,1)
    veg = veg + dummyveg  !add on emerging vegetation
     CALL InfiltProb(veg,m,n,K0,ie,rfx,rfy,kf,Kmax,dx,dy,infiltKern)
     !infiltKern =   infiltKern * slopeFactor

     WHERE (veg>0)
        flowResistance1 = flowResistance0 + kv
     ELSEWHERE
        flowResistance1 = flowResistance0
      END WHERE
  END IF

  WRITE(*,'("% Complete : ",F6.2)') Real(j)/REAL(nSteps) *100.0
END DO

 CLOSE(1)
 CLOSE(2)

END SUBROUTINE SimCODE

!Subroutine GAInfilt
SUBROUTINE GAInfilt (m,n,dt,inflow, Ksat, wfs, cumInfilt, infex)
!code for green-ampt infiltration calculation
IMPLICIT NONE 

INTEGER, INTENT(IN) :: m, n !dimensions of arrays
REAL*8, INTENT(IN) :: dt !time step
REAL*8, DIMENSION(m,n), INTENT(IN) :: inflow !surface runon + precipitation
REAL*8, DIMENSION(m,n), INTENT(IN) :: Ksat !saturated hydraulic conductivity
REAL*8,DIMENSION(m,n), INTENT(IN) :: wfs !wetting front suction * moisture deficit
REAL*8, DIMENSION(m,n), INTENT(OUT) :: infex !infiltration excess
REAL*8, DIMENSION(m,n), INTENT(INOUT) :: cumInfilt !cumulative infiltration

REAL*8 :: Fold2, Fnew !dummy variables for iterative calculation of infiltration
REAL*8 :: F1, F2 !initial and final cumulative infiltration in time step
REAL*8 :: ip, i  !infiltration excess and inflow rates
REAL*8 :: ff2,  f !final and initial infiltration rates
REAL*8 :: Ks, MS !sat hydraulic conductivity, wetting front suction
REAL*8 :: error !criteria for stopping iterations
REAL*8 :: tp, Fp ! time to ponding an cum infiltratio when ponding occurs
INTEGER :: x, y, cc !counters

DO x=1,m !loop through spatial domain
DO y=1,n

  i    = inflow(x,y)
  F1 = cumInfilt(x,y)
  Ks   = Ksat(x,y)
  MS = wfs(x,y)
  !IF ((x.eq.1).AND.(y.eq.1)) THEN
  !PRINT *,F1
  !END IF
  
  !calculation of initial infiltration capacity
  IF (F1.gt.1e-3) THEN
    f = Ks*(1.d0 + MS/F1)
  ELSE
    f = 100.0*Ks
  END IF
  
  If(i > f) THEN
 
  !instant ponding
    Fold2 = 0.0d0
    error = 1.d0
    cc=1
    DO While((error > 0.001d0).AND.(cc<1000))
      Fnew = F1 + Ks*dt + MS*Log((Fold2 + MS)/(F1 + MS))
      error = abs((Fnew - Fold2)/Fnew)
      Fold2 = Fnew
      cc = cc+1
    END DO
    If (cc.ge.1000) THEN
     PRINT *,"GAInflt recurrsion error 1"
    END IF
    F2 = Fnew
    ip = i - (F2 - F1)/dt
    !PRINT *,"1 ip, i, f, dt, F2, F1, tp",ip, i, f, dt, F2, F1, tp
  ELSE
    F2 = F1 + i*dt !cumulative infilt if all infiltrates
    ff2 = Ks*(1.d0 + MS/F2) 
    tp = 0.d0
    If(i < ff2) THEN !check if ponding occurs in time step
      ip = 0.d0
      tp = 0.d0
    ELSE
      Fp = Ks*MS /(i - Ks)
      tp = (Fp - F1)/i
      IF (tp>dt) THEN
        tp = 0.d0
      END IF
      Fold2 = 0.0 !first guess
      error = 1.d0
      cc=1
      DO While((error > 0.001d0).AND.(cc<1000))
        Fnew = F1 + Ks*(dt-tp) + MS*Log((Fold2 + MS)/(F1 + MS))
        error = (Fnew - Fold2)/Fnew
        Fold2 = Fnew
        cc = cc +1
      END DO
       If (cc.ge.1000) THEN
        PRINT *,"GAInflt ecurrsion error 2"
       END IF
      F2 = Fnew
      ip = (i*dt - (F2 - F1))/(dt - tp)
      !PRINT *,"ip, i, dt, F2, F1, tp",ip, i, dt, F2, F1, tp
    END IF
  END IF
  
  cumInfilt(x,y) = F2
  infex(x,y) = ip
 
 ! (* End of Green Ampt infiltration for time step *)

  
END DO
END DO

!DO x=1,m
!PRINT *, infex(x,:)
!END DO

END SUBROUTINE GAInfilt

!TwoDRandPos
SUBROUTINE TwoDRandPos(randOrder,m, n,mn) 
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: m, n, mn
  INTEGER, DIMENSION(mn,2), INTENT(INOUT) :: randOrder
  INTEGER :: i, j, k
    
  k=1
  DO i=1,m
    DO j=1,n
       randOrder(k,:) =(/ i, j /)
       k = k +1
    END DO
  END DO
  
  CALL OneDRandList(randOrder,mn)
end SUBROUTINE TwoDRandPos

!KWPrimer
SUBROUTINE KWPrimer(m,n,topog, manningsN, mask,solOrder,flowdirns,alpha,deltax, solMax)
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: m, n
  REAL*8, DIMENSION(m,n), INTENT(IN) :: topog, manningsN
  INTEGER, INTENT(OUT) :: solMax
  REAL*8, DIMENSION(m,n), INTENT(OUT) :: alpha, deltax
  INTEGER, DIMENSION(m,n),INTENT(OUT) :: solOrder, flowdirns
  INTEGER, DIMENSION(m,n,9),INTENT(OUT) :: mask
  
  
  INTEGER, DIMENSION(9,2) :: neighbs
  INTEGER :: x, y, i ,j, dirn, dx, dy
  REAL*8 :: cellLength = 1.d0, slope


  CALL GD8(topog,flowdirns,m,n)
  CALL KMWOrder(flowdirns,m,n,solOrder)
  solMax = maxval(solOrder)

  !iterate over space to generate calculation mask
  ! mask is a 9 x 1 array where 1s indicate that a neighbour (as defined by the position in mask(x,y) ) 
  !discharges to the cell at x,y
  DO x=1,m  
  DO y=1,n
    IF (solOrder(x,y).eq.1) THEN
      mask(x,y,:) = (/ 0,0,0,0,0,0,0,0,0/)
    ELSE
      CALL Neighbours(1,(/x,y/), (/m,n/),neighbs)
      DO j=1,9
        IF (neighbs(j,1).ne.-99) THEN
          dirn = flowdirns(neighbs(j,1),neighbs(j,2))
          CALL Lookupfdir(dirn, dx,dy)
          If(((neighbs(j,1)-dx).eq.x).and.((neighbs(j,2)-dy).eq.y)) THEN
            mask(x,y,j) = 1
          ELSE
            mask(x,y,j) = 0
          END IF
        ELSE
          mask(x,y,j) = 0
        END IF
      END DO
    ENDIF
    !WRITE(*,'(a5,i2,i2,a1,9i2)') 'mask(',x,y,')',mask(x,y,:)
    !calculate deltax(x,y)
    If((flowdirns(x,y).eq.1).or.(flowdirns(x,y).eq.3).or.(flowdirns(x,y).eq.5).or.(flowdirns(x,y).eq.7)) THEN
      deltax(x,y) = (2.d0**0.5d0) * cellLength
    ELSE
      deltax(x,y) = cellLength
    END IF
    deltax = deltax 
    
    !calculate alpha(x,y)
    CALL Lookupfdir(flowdirns(x,y),dx,dy)
   ! PRINT *, x,y,dx,dy
    IF ((dx==-2).AND.(dy==-2)) THEN
      alpha(x,y) = (ManningsN(x,y) * (cellLength**(0.6666667d0)) / (0.0001d0)**(0.5d0))**0.6d0
      !use a default slope of 0.0001
    ELSE
      IF ((dx==0).AND.(dy==0)) THEN
        alpha(x,y) = (ManningsN(x,y) * (cellLength**(0.6666667d0)) / (0.0001d0)**(0.5d0))**0.6d0
         !use a default slope of 0.0001
      ELSE
        slope = ((topog(x,y)-topog(x - dx,y-dy))/deltax(x,y))
        !WRITE(*,'(2i3,5e14.7)') x, y, slope, topog(x,y),topog(x - dx,y-dy),deltax(x,y),REAL(ndx)
        IF (abs(slope)<0.0001d0) THEN
          slope = 0.0001d0
        END IF
        !assuming water3 height is small compared to cell size
        alpha(x,y) = (ManningsN(x,y) * cellLength**0.6666667d0 / abs(slope)**0.5d0)**0.6d0
      END IF
    END IF
  
  END DO
  END DO
  
END SUBROUTINE KWPrimer
!KinematicWave
SUBROUTINE KinematicWave(m,n,ndx, dt, iex,flowdirns,solOrder,solMax,mask,alpha,deltax,disOld,disNew)
!This subroutine calculates the kinematic wave equation for surface runoff 
!in a network for a single time step
!The flow network is given by the GD8 algorithm calculated on the topog
!we assume no interaction between adjacent flow pathways i.e. water does 
!not overflow in a direction not specified by the GD8 network
!iex is the excess water from precip - infiltration + GW discharge
!iex(m,n,1) = qit
!iex(m,n,2) = qitPlus
!manningsN the roughness coefficient
!init a value passed to initialise the variables

IMPLICIT NONE

INTEGER, INTENT(IN) :: m, n, ndx, solMax
REAL*8, DIMENSION(m,n,2),INTENT(IN) :: iex
REAL*8, DIMENSION(m,n,ndx), INTENT(INOUT) :: disOld, disNew
REAL*8, INTENT(IN) :: dt !time step
INTEGER, DIMENSION(m,n),INTENT(IN) :: flowdirns, solOrder
INTEGER,DIMENSION(m,n,9),INTENT(IN) :: mask
REAL*8, DIMENSION(m,n),INTENT(IN) :: alpha, deltax

REAL*8 :: eta = 1.e-6, cellLength = 1.d0, beta = 0.6d0
REAL*8,DIMENSION(m,n,ndx) :: dummydisOld
REAL*8 :: QQi1t2,QQi2t1, QQi2t2
REAL*8 :: fQ22, dfQ22
REAL*8 :: qi2t2, qi2t1
INTEGER :: a, i, j, k, l, x, y
LOGICAL :: IsConverged

dummydisOld = disNew

!Calculation for one time step
  DO a=1,solMax  !iterate over the solution order
    DO x=1,m  !iterate over space
    DO y=1,n
      IF (solOrder(x,y).eq.a) THEN
      
        DO l=1,ndx-1  !iterate over  deltax
          qi2t2 = iex(x,y,2) ! check
          qi2t1 = iex(x,y,1) !same ?
          IF (maxval(mask(x,y,:)).gt.0) THEN
            IF(l.eq.1) THEN
              QQi1t2 = 0.d0
              k = 0
              DO i=-1,1
              DO j = -1,1
                k = k + 1
                If(mask(x,y,k).gt.0) THEN
                  QQi1t2 = QQi1t2 + disNew(x+i,y+j,ndx)
		  !PRINT *,'(x,y)',x,y,',(x+i,y+j)=',x+i,y+j, '(k)=',k
                END IF
              END DO
              END DO
            ELSE
              QQi1t2 = disNew(x,y,l)
            END IF
          ELSE
            QQi1t2 = disNew(x,y,l)
          END IF
  
          QQi2t1 = disOld(x,y,l+1) 
          !initial guess from linear scheme
          QQi2t2 = (dt/(deltax(x,y)/real(ndx))*QQi1t2 + alpha(x,y)*beta*QQi2t1*&
                    ((QQi2t1 + QQi1t2)/2.d0)**(beta - 1) + &
                    dt*(qi2t2 + qi2t1)/2.d0)/(dt/ (deltax(x,y)/real(ndx)) + &
                    alpha(x,y)*beta*((QQi2t1 + QQi1t2)/2.d0)**(beta - 1.d0))
          !nonlinear KMW discretization
          !PRINT *,'QQi2t2',QQi2t2
          isConverged=.FALSE.
          DO While (isConverged .eqv. .FALSE.)
            
            fQ22 = dt/ (deltax(x,y)/real(ndx))*QQi2t2 + alpha(x,y)*(QQi2t2**beta) &
                    - (dt/ (deltax(x,y)/real(ndx))*QQi1t2 + alpha(x,y)*(QQi2t1**beta) + & 
                    dt*(qi2t2 + qi2t1)/2.d0)
            dfQ22 = dt /  (deltax(x,y)/real(ndx)) + alpha(x,y)*beta*(QQi2t2**(beta - 1.d0))
            !IF ((QQi2t2 - fQ22/dfQ22).lt.0) THEN
            ! PRINT *, (QQi2t2)
            !END IF
            QQi2t2 =  QQi2t2 - fQ22/dfQ22
            
            isConverged=Abs(fQ22).le.eta
             IF (QQi2t2.lt.0) THEN
             !PRINT *, QQi2t2, QQi1t2, QQi2t1, qi2t2, qi2t1
             !READ*,
              QQi2t2 = QQi1t2
              isConverged=.TRUE.
            END IF
            !problem here.  The outflow cell causes Nans .
            !IF ((x.eq.y)) THEN
            !  PRINT *,' fQ22, QQi2t2', fQ22,QQi2t2
            !end IF
           
          END DO
         disNew(x,y,l+1) = QQi2t2
        END DO
      END IF
    END DO
    END DO   
  END DO
  !disNew = dummydisNew
  disOld = dummydisOld

END SUBROUTINE KinematicWave
!KMWOrder
SUBROUTINE KMWOrder(flowdirns,m,n,solutionOrder)
! subroutine assigns values to the matrix solutionOrder to tell the Kinematic Wave subroutine
!in which order to solve the KM equation on the drainage network
!values of 1 assigned to top of catchment, 2 to 1st downstream node etc.
! value at a point is the maximum travel distance to that point from all points above it 

IMPLICIT NONE

INTEGER, INTENT(IN) :: m, n
INTEGER, DIMENSION(m,n), INTENT(IN) :: flowdirns
INTEGER, DIMENSION(m,n), INTENT(OUT) :: solutionOrder

LOGICAL :: isRoute
INTEGER :: x, y, dx, dy, i, j, particle
INTEGER, DIMENSION(m,n) :: usca
solutionOrder = 1

DO i=1,m
DO j=1,n  !loop though positions in matrix
  x=i
  y=j
  isRoute = .TRUE.
  particle = 1
  DO WHILE(isRoute)
    CALL Lookupfdir(flowdirns(x,y),dx,dy)
    IF ((dx.eq.-2).AND.(dy.eq.-2)) THEN
      isRoute = .FALSE.
    ELSE 
      IF ((dx.eq.0).AND.(dy.eq.0)) THEN
        isRoute = .FALSE.
      ELSE
        x = x - dx
        y = y - dy
        particle = particle + 1
        solutionOrder(x,y) =  max(solutionOrder(x,y),particle)
        IF (flowdirns(x,y).eq.-1) THEN
           isRoute = .FALSE.
        END IF
       !WRITE(*,'(a5,i3,a5,i3,a5,i3,a5,i3,a5,i3,a5,i3,a3)') 'x',x+dx,'y',y+dy,'newx',x ,'newy', y ,'dx',dx,'dy',dy,isRoute
      ENDIF
    ENDIF
  END DO
END DO
END DO

END SUBROUTINE KMWOrder
!OneDRandList
SUBROUTINE OneDRandList(a,mn)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: mn
  integer, DIMENSION(mn,2), intent(inout) :: a
  integer :: i, randpos, temp(2)
  real :: r
  
  do i = mn, 2, -1
    call random_number(r)
    randpos = int(r * i) + 1
    temp = a(randpos,1:2)
    a(randpos,1:2) = a(i,1:2)
    a(i,1:2) = temp
  end do
 
END SUBROUTINE OneDRandList
!Lookupfdir
SUBROUTINE Lookupfdir(dirn, dx,dy)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: dirn
  INTEGER, INTENT(OUT) :: dx, dy
  
  SELECT CASE(dirn)
    CASE(1)
      dx = -1
      dy = -1
    CASE(2) 
      dx = -1
      dy =  0
    CASE(3) 
      dx = -1
      dy =  1
    CASE(4) 
      dx = 0
      dy = 1
    CASE(5) 
      dx = 1
      dy = 1
    CASE(6) 
      dx = 1
      dy = 0
    CASE(7) 
      dx = 1
      dy = -1
    CASE(8) 
      dx = 0
      dy =  -1
    CASE(9) 
      dx = 0
      dy = 0
    CASE DEFAULT
      dx = -2
      dy =  -2
  END SELECT
END SUBROUTINE Lookupfdir
!RoutingWithKernel
SUBROUTINE RoutingWithKernel(m,n,mn,precip, infiltKern, storeKern, newflowdirns, topog, store,discharge,outflow)

IMPLICIT NONE

INTEGER, INTENT(IN) :: m, n, mn
INTEGER, INTENT(INOUT) :: outflow
REAL*8, DIMENSION(m,n), INTENT(IN) :: infiltKern,  storeKern , topog
INTEGER, DIMENSION(m,n), INTENT(IN) :: precip 
INTEGER, DIMENSION(m,n),INTENT(INOUT) :: store,newflowdirns
INTEGER, DIMENSION(m,n),INTENT(OUT) :: discharge

INTEGER, DIMENSION(m,n) :: tempPrecip, tempSurfaceStore
INTEGER, DIMENSION(mn,2) :: randOrder
LOGICAL :: isRoute
INTEGER :: x, y, dx, dy, i, j, k, nIter
REAL :: rnum
INTEGER, DIMENSION(9,2) :: neighbs

!initialise parameters

tempSurfaceStore = 0
discharge = 0
!outflow = 0
!PRINT *, "in Routing", infiltKern
!READ*,
CALL TwoDRandPos(randOrder,m,n,mn)  !randomly shuffel positions
nIter = maxval(precip)

donumParticles : DO i=1,nIter  !loop over nIter rainfall particles
  !PRINT *, "Precip number = ",i
  dodomain : DO j=1,mn  !loop though positions in matrix
    x = randOrder(j,1)
    y = randOrder(j,2)
    tempPrecip(x, y) = tempPrecip(x, y) - 1
    isRoute = .TRUE.

 dowhileisRoute :  DO WHILE(isRoute)
    !PRINT *, "x,y", x,y
    !READ*,
      CALL random_number(rnum)
     Firstif : IF ((rnum<infiltKern(x,y)).AND.(store(x,y)<storeKern(x,y))) THEN
        store(x,y) = store(x,y) + 1
        isRoute = .FALSE.
      ELSE
       ifnotPeriodicBC : IF (newflowdirns(x,y).ne.-3) THEN !if not a periodic boundary condition then
       CALL Lookupfdir(newflowdirns(x,y),dx,dy)
        discharge(x,y) = discharge(x,y) + 1
        ifoutflow : IF ((dx==-2).AND.(dy==-2)) THEN
          outflow = outflow + 1
          isRoute = .FALSE.
        ELSE 
          ifdepression : IF ((dx==0).AND.(dy==0)) THEN
            !need to better account for depressions
            tempSurfaceStore(x,y) = tempSurfaceStore(x,y) + 1
            CALL Neighbours(1,(/x,y/), (/m,n/),neighbs) !find next lowest cell
            dy = 1
            rnum = topog(neighbs(1,1),neighbs(1,2))+ 0.001*dble(tempSurfaceStore(neighbs(1,1),neighbs(1,2)))
            DO k=2,9
              IF ((neighbs(k,1).ne.-99).AND.(k.ne.5)) THEN
                IF ((topog(neighbs(k,1),neighbs(k,2))+ 0.001*dble(tempSurfaceStore(neighbs(k,1),neighbs(k,2)))).lt.rnum) THEN
                    dy = k
                    rnum = topog(neighbs(k,1),neighbs(k,2))+ 0.001*dble(tempSurfaceStore(neighbs(k,1),neighbs(k,2)))
                END IF
              END IF
            END DO
            
            ifoverflow : IF ((topog(x,y)+ dble(tempSurfaceStore(x,y))*0.004)>rnum) THEN
              !assumed 1 particle = 4 mm
              CALL fdirLookup((/x-neighbs(k,1) , y-neighbs(k,2)/), dy)
              newflowdirns(x,y)  =  dy !value reassigned 
              discharge(x,y) = discharge(x,y) + 1
              tempSurfaceStore(x,y) =  tempSurfaceStore(x,y) - 1
              x = neighbs(k,1)
              y = neighbs(k,2)
            ELSE 
              isRoute = .FALSE.
            END IF ifoverflow
            
          ELSE
            x = x - dx
            y = y - dy
            IF (newflowdirns(x,y).eq.-1) THEN
              isRoute = .FALSE.
            END IF
          ENDIF ifdepression
        ENDIF ifoutflow
      
      ELSE
        discharge(x,y) = discharge(x,y) + 1
        !periodic boundary conditionas can only really be defined simply for an inclined plane
        !with two adjacent edges defined as the boundary 
        !from which particles are routed to the opposite boundary
        !for simplicity I'll assume the landscape slopes downwards in the direction of lower x and y
       IF ((x.eq.1).AND.(y.gt.1)) THEN
         x = m
         y = y
       END IF
       IF ((y.eq.1).AND.(x.gt.1)) THEN
         x = x
         y = n
       END IF
       IF ((y.eq.1).AND.(x.eq.1)) THEN
         x = m
         y = n
       END IF
        
      END IF ifnotPeriodicBC
      END IF Firstif
    END DO dowhileisRoute
  END DO dodomain
CALL OneDRandList(randOrder,mn)
END DO donumParticles

DO i=1,m
DO j=1,n
  IF ((tempSurfaceStore(i,j)>0).AND.(store(i,j)<storeKern(i,j))) THEN
    store(i,j) = min(store(i,j)+tempSurfaceStore(i,j),INT(storeKern(i,j)))
  END IF
END DO
END DO

END SUBROUTINE RoutingWithKernel
!InfiltProb
SUBROUTINE InfiltProb(veg,m,n,K0,ie,rfx,rfy,kf,Kmax,dx,dy,iProb)
  IMPLICIT NONE         
  
  INTEGER, INTENT(IN) :: m, n  !dimensions of the spatial arrays
  REAL*8, INTENT(IN) :: rfx, rfy, kf !maximum radius for plant effects on soil properties 
  INTEGER, DIMENSION(m,n), INTENT(IN) :: veg  !vegetation matrix
  REAL*8, INTENT(IN) :: K0, ie, Kmax  !parameters defining infiltration properties
  REAL*8, INTENT(IN) :: dx,dy !length scales of lattice
  REAL*8, DIMENSION(m,n), INTENT(OUT) :: iProb  !infiltration probability matrix
  
  REAL*8, DIMENSION(-1*int(Floor(rfx/dx)):int(Floor(rfx/dx)),-1*int(Floor(rfy/dy)):int(Floor(rfy/dy))) :: kern  !infiltration kernel
  REAL*8 :: radius
  Integer :: i, j, di, dj, m1, n1 !counters
 
  m1 = int(Floor(rfx/dx))
  n1 = int(Floor(rfy/dy))
  !PRINT *, "in infiltProb, K0, ie, rfx,rfy,alpha",K0, ie, rfx,rfy,alpha
  !READ*,
  DO i=-1*m1,m1
    Do j=-1*n1,n1
      radius = ((i*dx)**2.0d0+(j*dy)**2.0d0)**(0.5d0)
      IF (radius.le.rfx) THEN
         kern(i,j) = dexp(-1.d0*(kf**2.d0)*(radius**2.d0))
      ELSE
         kern(i,j) = 0.d0
      END IF
    END DO
  END DO
  
  kern = kern / sum(kern)  *Kmax/ie
 
  
   iProb = K0/ie  !base level infiltration probability
   DO i=1,m
   DO j =1,n
     DO di = -1*m1,m1
     DO dj = -1*n1,n1
        IF (veg(Modulo(i + di - 1, m) + 1,Modulo(j + dj - 1, n) + 1)>0) THEN
            iProb(i,j) =  iProb(i,j) + kern(-1*di,-1*dj)
        END IF
     END DO
     END DO
   ENd DO
   END DO
  !PRINT*, "iprob",iProb
  !READ*,
END SUBROUTINE InfiltProb
!ListConvolve
SUBROUTINE ListConvolve(base,kernel,convol,m,n,m1,n1)
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: m, n, m1, n1
  INTEGER, DIMENSION(m,n), INTENT(IN) :: base
  REAL*8, DIMENSION(m1,n1), INTENT(in) :: kernel
  REAL*8, DIMENSION(m,n), INTENT(OUT) :: convol
  
  REAL*8,DIMENSION(:,:),ALLOCATABLE :: dummybase,dummybase2
  INTEGER :: m2, n2, i, j, dm,dn
  
  dm = floor(real(m1)/2.0)
  dn = floor(real(n1)/2.0)
  m2 = m+m1-1
  n2 = n+n1-1
  !PRINT *,"m2,n2",m2,n2
  ALLOCATE(dummybase(m2,n2))
  ALLOCATE(dummybase2(m2,n2))
  
 
  !creat dummy veg matrix with boundaries a reflection of inside
  dummybase = 0.d0
  dummybase2 = 0.d0
  dummybase(dm+1:m+dm,dn+1:n+dn)=DBLE(base)/1.d00
  
  DO i=1,dm
      dummybase(dm+1-i,dn+1:n+dn)=base(i,:)
      dummybase(dm-dm+i,dn+1:n+dn)=base(m-i+1,:)
  END DO
  DO i=1,dn
      dummybase(dm+1:m+dm,dn+1-i)=base(:,i)
      dummybase(dm+1:m+dm,n2-dn+i)=base(:,n-i+1)
  END DO
 
  DO i=dm+1,m2-dm
    DO j=dn+1,n2-dn
        dummybase2(i-dm:i+dm,j-dn:j+dn)= dummybase2(i-dm:i+dm,j-dn:j+dn)+dummybase(i,j)*kernel
    END DO
  END DO
  
  
  convol = dummybase2(dm+1:m+dm,dn+1:n+dn)
  
 DEALLOCATE(dummybase)
 DEALLOCATE(dummybase2)
END SUBROUTINE ListConvolve
!Erosion   
SUBROUTINE Erosion(discharge,topog,newflowdirns,veg,flowResistance,m,n)
  IMPLICIT NONE
  Integer, INTENT(IN) :: m,n
  Integer, DIMENSION(m,n), INTENT(IN) :: discharge, newflowdirns,veg
  REAL*8, DIMENSION(m,n), INTENT(IN) :: flowResistance
  REAL*8, DIMENSION(m,n), INTENT(INOUT) :: topog
  
  Real*8, DIMENSION(m,n) :: flow,  dx, deltaz, Qs_in, Qs_out
  INTEGER :: holes
  REAL*8 :: time !erosion resistance parameters
  REAL*8 :: timestep, bulkdensity, porosity, celllength, wbar, specificweight, slope
  INTEGER ::  i, j, fDirRef, deltax,deltay, posx,posy
  LOGICAL :: test1
 
  !Parameters
  timestep = 30.d0*24.0d0*60.d0*60.d0  !  (* effective flow duration in seconds p*)
  time = timestep
  bulkdensity = 2650.d0 ! (* (kg/m3) *) 
  porosity = 0.4d0  
  cellLength = 1.d0
  wbar = cellLength !; (* lets assume a fixed width of flow *)
  specificweight = (1.d0 - porosity)*bulkdensity*wbar
 
  
  !Initialise variables
  flow = DBLE(discharge)*(cellLength)**2/(1000.d0)/(365.d0*24.d0*60.d0*60.d0) !  conversion of mm/cell/year to m^3/sec
  
  Qs_in = 0.d0  ! (* sediment flux in *)
  Qs_out = 0.d0 !sediment flux out
  dx = 0.d0 !flow length, dx
   
  Do i=1,m
  DO j=1,n
    fDirRef = newflowdirns(i,j)
    CALL LookUpFdir(fDirRef,deltax,deltay)
    
    If ((fDirRef.gt.0).and.(fDirRef.lt.9)) THEN
     posx=i-deltax
     posy=j-deltay
     dx(i, j) = (deltax**2.d0 + deltay**2.d0)**0.5
     slope = (topog(i, j) - topog(posx, posy))/dx(i, j)/cellLength

     Qs_out(i, j) = 7000.d0*(slope**1.5d0) * (cellLength**2.d0)*flow(i, j) / flowResistance(i, j)**0.5d0
     Qs_in(posx, posy) = Qs_out(i, j) + Qs_in(posx, posy)
     !PRINT*, "posx, posy" , posx, posy,"i, j", i, j , "slope", slope, "dx", dx(i,j), "Qs_out(i, j)",Qs_out(i, j)
    ELSEIF (fDirRef.le.0) THEN
      dx(i, j) = cellLength
      slope = 1.d0
      Qs_out(i, j) = 7000.d0*(slope**1.5d0) * (cellLength**2.d0)*flow(i, j) / flowResistance(i, j)**0.5d0
    ELSE
      dx(i, j) = cellLength
      Qs_out(i, j) = 0.d0
    END IF
    
  END DO
  END DO
  
  !boundary condition if ouflow cells undergoes no lowering
  DO i=1,m
  DO j=1,n
    fDirRef = newflowdirns(i,j)
    IF (fDirRef.le.0) THEN  
      Qs_out(i, j) = Qs_in(i, j)
    END IF
  END DO
  END DO
  
  

  !(* calc of change in elevation dz *)
  deltaz = 0.d0
  Do i=1,m
  DO j=1,n
    fDirRef = newflowdirns(i,j)
    If(fDirRef.eq.0) THEN 
     deltaz(i, j) = -1.d0*(( - Qs_in(i, j))/dx(i, j)/specificweight)
    ELSEIF (fDirRef.eq.-1) THEN 
      deltaz(i, j) = -0.05d0 / timestep  !base level lowering rate
    ELSE
     deltaz(i, j) = -1.d0*((Qs_out(i, j) - Qs_in(i, j))/dx(i, j)/specificweight)
    END IF
  END DO
  !PRINT *, deltaz(i,:)
  !PRINT*,
  !PRINT*,topog(i,:)+deltaz(i,:)
  !  deltaz(i,1) = 0.0d0  !base level lowering rate
  END DO
  !READ*,
  
  
  deltaz = deltaz * timestep
  topog = topog + deltaz
  
END SUBROUTINE Erosion
!Splash
SUBROUTINE Splash(topog,veg, Dv, Db,m,n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: m, n
  INTEGER, DIMENSION(m,n), INTENT(IN) :: veg
  REAL*8, DIMENSION(m,n), INTENT(INOUT) :: topog
  REAL*8, INTENT(IN) :: Dv, Db  !m^2/kyr
  
  INTEGER :: i, j, k
  REAL*8, DIMENSION(m,n) :: D, deltaz !D = diffusion coeff
  INTEGER, DIMENSION(9,2) :: neighbs
  REAL*8 :: dx, deltat, cellLength
  
  cellLength = 1.0d0
  deltaz = 0.0d0
  deltat = (1.0e-3)*30.d0/365.d0 !kyr
  
  WHERE (veg>0) 
    D = Dv
  ELSEWHERE
    D = Db
  END WHERE
  
  DO i=1,m
  DO j=1,n
    CALL Neighbours(1,(/i,j/), (/m,n/),neighbs)
    DO k = 1,9
      IF (k.ne.5) THEN
        IF (neighbs(k,1).ne.-99) THEN
          dx = ((i-neighbs(k,1))**2 + (j-neighbs(k,2))**2)**0.5
          deltaz(i,j) = deltaz(i,j) + deltat*D(neighbs(k,1),neighbs(k,2))* &
          (topog(neighbs(k,1),neighbs(k,2))-topog(i,j))/((dx*cellLength)**2)
        END IF
      END IF
    END DO
  END DO
  END DO
  
  topog = topog + deltaz
  
  
END SUBROUTINE Splash
!FindHoles
SUBROUTINE FindHoles(newtopog,holes)
IMPLICIT NONE

REAL*8, DIMENSION(:,:), INTENT(IN) :: newtopog
INTEGER, INTENT(OUT) :: holes

INTEGER :: m,n, i,j,k,l

m=SIZE(newtopog,1)
n=SIZE(newtopog,2)

DO i=2,m-1
DO j=2,n-1
  holes=0
  DO k=-1,1
  DO l=-1,1
    IF (newtopog(i,j)<newtopog(i + k, j + l)) THEN
      holes=holes+1
    END IF
  END DO
  END DO
  IF (holes.ge.8) THEN
    holes = 1
    RETURN
  END IF
ENDDO
END DO

END SUBROUTINE FindHoles
!Evaporation
SUBROUTINE Evaporation(veg,eTActual,bareE,store,tsteps,rcx,rcy,kc,dx,dy,params)
  !This version cycles through sites and evaporates water from site and neighbouring sites if vegetated
  IMPLICIT NONE
   
  INTEGER, DIMENSION(:,:),  INTENT(IN) :: veg
  INTEGER, DIMENSION(:,:),  INTENT(INOUT) :: store
  INTEGER, DIMENSION(:,:),  INTENT(INOUT) :: eTActual
  INTEGER, DIMENSION(:,:),  INTENT(INOUT) :: bareE
  REAL*8, INTENT(IN) :: rcx,rcy, kc
  INTEGER, INTENT(IN) :: tsteps
  REAL*8, DIMENSION(7),  INTENT(IN) :: params
  REAL*8, INTENT(IN) :: dx, dy
  
  REAL*8 :: beta,te,pbar,Psb,Psv
  INTEGER :: n, m, mn, i, j, k, a, b, di,dj, m1, n1, storeCounter
  REAL*8 :: rnd, baresoilprob, radius
  REAL*8, DIMENSION(-1*int(Floor(rcx/dx)):int(Floor(rcx/dx)),-1*int(Floor(rcx/dy)):int(Floor(rcx/dy))) :: uptakeprob
  INTEGER, DIMENSION(SIZE(store),2) :: posxy
  storeCounter = sum(store)
  !dx = params(1)
  !dy = params(2)
  te = params(3)
  pbar = params(4)
  Psb = params(5)
  Psv = params(6)
  beta = params(7)
  
  n = SIZE(store,1) 
  m = SIZE(store,2)
  mn = SIZE(store)
  
  m1 = int(Floor(rcx/dx))
  n1 = int(Floor(rcy/dy))
 
  Do i=(-1*m1),m1,1
  DO j =(-1*n1),n1,1
     radius = sqrt((i*dx)**2.0+(j*dy)**2.0)
     IF (radius<=rc) THEN
        uptakeprob(i,j) = dexp(-1.d0*(kc**2.0)*(radius**2.0))
     ELSE
        uptakeprob(i,j) = 0.d0
     END IF
  END DO
  END DO
  
   uptakeprob= uptakeprob/sum(uptakeprob)*dx*dx*beta*te/pbar
  k=1
  Do i=1,m
  Do j=1,n
    posxy(k,:)=(/i,j/)
    k=k+1
  END DO
  END DO
  
  !PRINT *, "in evap, tsteps ",tsteps
  i=1
  Do WHILE ((i.le.tsteps).AND.(storeCounter.gt.0))
     i=i+1
    CALL TwoDRandPos(posxy,m, n,mn) 
    Do j =1, mn
     a = posxy(j,1)
     b = posxy(j,2)
     
     If (veg(a,b)>0) THEN
        
        DO di=-1*m1,m1
        DO dj=-1*n1,n1
          CALL random_number(rnd)
          IF ((store(Modulo(a + di - 1, m) + 1,Modulo(b + dj - 1, n) + 1).ne.0).and. &
            (rnd<uptakeprob(-1*di,-1*dj))) THEN
              eTActual(a,b)=eTActual(a,b)+1
	      storeCounter=storeCounter-1
              store(Modulo(a + di - 1, m) + 1,Modulo(b + dj - 1, n) + 1)= &
              store(Modulo(a + di - 1, m) + 1,Modulo(b + dj - 1, n) + 1)-1
          END IF
       END DO
       END DO
     END IF
     
     If (veg(a,b)>0) THEN
       baresoilprob = Psv
     ELSE
       baresoilprob = Psb
     END IF
     
     CALL random_number(rnd)
     IF ((store(a,b).ne.0).and.(rnd<baresoilprob)) THEN
        bareE(a,b) = bareE(a,b)+ 1
        store(a,b)=store(a,b)-1
	storeCounter = storeCounter-1
      ENDIF
     END DO
  END DO
  
END SUBROUTINE Evaporation
!Neighbours
SUBROUTINE Neighbours(order,posij, dom,neighbs)
  ! (* returns the positions of neighbouring cells depending upon the extent of the neighbourhood  *)
  ! (* order=1 is the 1st moore neighbourhood *)
  !(retunrs (-99,-99) for positions outside domain
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: order
  INTEGER, DIMENSION(2), INTENT(IN) :: posij, dom
  INTEGER, DIMENSION((order*2+1)**2,2), INTENT(out) :: neighbs
  
  INTEGER :: xmax, ymax, x0, y0, i, j, k
  
  neighbs = -99
  xmax = dom(1)
  ymax = dom(2)
  x0 = posij(1)
  y0 = posij(2)
 
     k=1
     Do i=(-1*order),order 
     DO j=(-1*order),order
      If (((x0+i).ge.1).and.((x0+i).le.xmax).and.((y0+j).ge.1).and.((y0+j).le.ymax)) THEN
        neighbs(k,1)=x0+i
        neighbs(k,2)=y0+j
      END IF
      k=k+1
     END DO
     END DO
    
   END SUBROUTINE Neighbours
!LSDs   
SUBROUTINE LSDs(order, posxy, topog,m,n,lsdList)
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: order,m ,n
      INTEGER, DIMENSION(2), INTENT(IN) :: posxy
      REAL*8, DIMENSION(m,n), INTENT(IN) :: topog
      INTEGER, DIMENSION(3,2), INTENT(OUT) :: lsdList
      
      !(* This function returns the matrix positions;
      !LSD1: the position of the neighbouring cell with the steepest slope downhill;
      !LSD2: the position of the neighbouring cell with second steepest slope adjacent LSD1;
      !otherpos: the position of the other neighbouring cell the mirror reflection about LSD1 of LSD2;
      
      REAL*8 :: z0, slope1, slope2
      INTEGER :: ymax, xmax, n1,i,posy
      INTEGER, DIMENSION(1) :: posx
      INTEGER, DIMENSION(2) :: dom, LSD1, LSD2,clockpos, anticlockpos, otherpos
      INTEGER, DIMENSION((order*2+1)**2,2) :: pos2
      INTEGER, DIMENSION((order*2+1)**2) :: locs
      REAL*8, DIMENSION((order*2+1)**2-1) :: slps
      INTEGER, DIMENSION((order*2+1)**2-1,2) :: pos3
      
      ymax = SIZE(topog,1)
      xmax = SIZE(topog,2)
      dom = (/ymax, xmax /)
      n1 = 8*order
      
    If(posxy(1).ne.-99)  THEN 
     z0 = topog(posxy(1), posxy(2))
     CALL Neighbours(order, posxy, dom, pos2)
     pos3=pos2((/1,2,3,6,9,8,7,4/),:)
     !DO i=1,n1
     !PRINT *, "pos3", pos3(i,:)
     !END DO
     If(maxval(pos3) .eq.-99) THEN
       lsdList=RESHAPE((/ -99,-99,-99,-99,-99,-99 /),(/3,2/))  
       RETURN
     END IF
     Do i=1,n1
      If(pos3(i,1).ne.-99) THEN 
        slps(i) = 1.0*(z0 - topog(pos3(i,1), pos3(i,2)))/((posxy(1) - pos3(i,1))**2 + (posxy(2) - pos3(i,2))**2)**0.5
      ELSE
        slps(i) = -99999.d0
      end if
      !  PRINT*, "slps",slps(i)
     end do 
     
     !PRINT *, "maxval(slps)",maxval(slps)
     If(maxval(slps) < 0) THEN
      lsdList=RESHAPE((/ -99,-99,-99,-99,-99,-99 /),(/3,2/))     ! (* i.e. all neighbouring cells higher *)
      RETURN
     ELSE
      posx = maxloc(slps)
      !PRINT *, "posx",posx
      LSD1 = pos3(posx(1),:)
      !PRINT *,"LSD1",LSD1
      !CALL Pos1d(pos2,(order*2+1)**2,2,LSD1,locs)
      !posy = locs(1)
      !Do i=1,size(locs)
      !PRINT *, "locs",locs(i)
      !END DO
      ! PRINT *,"posy",posy
      CALL RotateArray(pos3,(order*2+1)**2-1,2,-1)
      clockpos = pos3(posx(1),:)
      ! PRINT *, "clockpos", clockpos
      CALL RotateArray(pos3,(order*2+1)**2-1,2,2)
      !CALL RotateArray(pos2,(order*2+1)**2,1,2)
      anticlockpos = pos3(posx(1),:)
      ! PRINT *, "anticlockpos", anticlockpos
       
      If((clockpos(1) .ne. -99).and.(anticlockpos(1) .ne. -99)) THEN
       Slope1 = 1.0*(z0 - topog(clockpos(1), clockpos(2)))/ &
       ((posxy(1) - clockpos(1))**2 + (posxy(2) - clockpos(2))**2)**0.5
      !PRINT *, "slope1", slope1
       Slope2 = 1.0*(z0 - topog(anticlockpos(1), anticlockpos(2)))/ &
       ((posxy(1) - anticlockpos(1))**2 + (posxy(2) - anticlockpos(2))**2)**0.5
      !PRINT *, "slope2", slope2
       If(Slope1 > Slope2) THEN
         LSD2 = clockpos
         otherpos = anticlockpos
       ELSE
         LSD2 = anticlockpos
         otherpos = clockpos
        END IF !(* end of If Slope1>Slope2 *)
         !PRINT *, "LSD2", LSD2, "otherpos", otherpos
      ELSE
       If((clockpos(1) .eq.-99).and.(anticlockpos(1).eq.-99)) THEN
         LSD2 = (/-99,-99/)
       ELSE
         If(clockpos(1).eq.-99) THEN
           LSD2 = anticlockpos
           otherpos = clockpos
         ELSE
           LSD2 = clockpos
           otherpos = anticlockpos
         END IF
       END IF
       END IF ! (* end of If clockpos!={"Null","Null"}anticlockpos!={"Null","Null"}, *)
      !PRINT *,"here"
    END IF !(* end of If Max[slps]<0 *)
     END IF !  end of if posxy!=(/-99,-99/)
  !PRINT *,"LSD1",LSD1
!PRINT *,"LSD2",LSD2
!PRINT *,"otherpos",otherpos
  lsdList(1,:) = LSD1
 ! PRINT *,"here"
  lsdList(2,:)= LSD2
  !PRINT *,"here"
  lsdList(3,:)=otherpos
  !PRINT *,"here"
   END SUBROUTINE LSDs
!RotateArray   
SUBROUTINE RotateArray(list,m,n,leftorRight)
     IMPLICIT NONE
     
     INTEGER, INTENT(IN) :: m,n
     INTEGER, DIMENSION(m,n), INTENT(INOUT) :: list
     INTEGER, INTENT(IN) :: leftorRight
     
     !Roates the rows of a 2 D array to the left<0 or
     !vice versa if rotate right > 0
     !if left = -1 then rotate a list left
     !if left != -1 then rotate right
     
     INTEGER,DIMENSION(n) :: temp
     INTEGER :: i
     
     DO i=1,abs(leftorRight)
     If(leftorRight.lt.0) THEN
       temp = list(1,:)
       list(1:m-1,:)=list(2:m,:)
       list(m,:)=temp
     ELSE
      IF (leftorRight.gt.0) THEN
       temp = list(m,:)
       list(2:m,:)=list(1:m-1,:)
       list(1,:)=temp
      END IF
     ENDIF
     END DO
   END SUBROUTINE RotateArray
!Pos1d   
SUBROUTINE Pos1d(list,m,n,match,rownum)
     IMPLICIT NONE
     !return the row position that is the same as match , returns 0,0 if no match 
     INTEGER, INTENT(IN) :: m,n
     INTEGER, DIMENSION(m,n),INTENt(in) :: list
     integer, dimension(n), intent(in) :: match
     INTEGER, DIMENSION(m),INTENT(OUT) :: rownum

    Integer :: i, j,k
    INTEGER :: test
    rownum=0
    k=1
    !PRINT*, "match", match
    DO i=1,m
    test = 1
    DO j=1,n
    !PRINT*, "list",list(i,:)
    IF (list(i,j).eq.match(j)) THEN
     test=test*1
    ELSE
     test = test*0
    END IF
    END DO
      IF(test.eq.1) THEN
        rownum(k)=i
        k=k+1
      END IF
    END DO
    
   END SUBROUTINE Pos1D 
!GD8   
Subroutine GD8(topog,flowdirns, m,n)
   
   IMPLICIT NONE
   
   INTEGER,INTENT(IN) :: m,n
   !INTEGER, DIMENSION(m,n), INTENT(IN) :: lakes
   REAL*8, DIMENSION(m,n), INTENT(IN) :: topog
   INTEGER, DIMENSION(m,n), INTENT(INOUT) :: flowdirns
   
   
   INTEGER :: mn, ordsCounter, order, counter, i,j,k, idirn, dx ,dy
   INTEGER, DIMENSION(SIZE(topog),2) :: ords
   REAL *8, DIMENSION(SIZE(topog)) :: oneDtopog
   INTEGER, DIMENSION(SIZE(topog)) :: ind2, rownum
   INTEGER, DIMENSION(2) :: startcell,  upstreamcell, upflowdirn
   INTEGER, DIMENSION(2) :: LSD1, secondary1,otherAnglePosition1
   INTEGER, DIMENSION(2) :: LSD2, secondary2,otherAnglePosition2
   INTEGER, DIMENSION(2) :: upsecondary, targetcell, dirn11, dirn12, dirn21, dirn22
   INTEGER, DIMENSION(3,2) :: lsdList
   REAL*8 :: z02, z0, z1, z2, slp1, slp2
   LOGICAL :: test1, test2, test3
   INTEGER :: isPeriodic  !flag for periodic boundary conditions
   isPeriodic = 1
   
   mn = m*n
   flowdirns = -2  !default value use for checking is a dirn has been assigned yet

   z02 = Floor(Minval(topog)) - 10
   order = 1
   CALL makeOrds(topog,ords, m,n)
   startcell = ords(1,:)
   ordsCounter = 2
   CALL Pos1d(ords,mn,2,startcell,rownum)
   IF (rownum(1).ne.0) THEN
      CALL endShift(ords,rownum(1), mn,2)
      ords(mn,:)=(/0,0/)
  END IF
     
   CALL LSDs(1, startcell, topog,m,n,lsdList)
   LSD1 = lsdList(1,:)
   secondary1 = lsdList(2,:)
   otherAnglePosition1 = lsdList(3,:)
   !PRINT *, " LSD1, secondary1, otherAnglePosition1", LSD1, secondary1, otherAnglePosition1
   !PRINT *, "startcell-LSD1", startcell-LSD1
   !READ*,
   
   CALL fdirLookup(startcell-LSD1, idirn)
   !PRINT *, "idirn",idirn
   flowdirns(startcell(1), startcell(2)) = idirn
   CALL Pos1d(ords,mn,2,startcell,rownum)
   IF (rownum(1).ne.0) THEN
      CALL endShift(ords,rownum(1), mn,2)
      ords(mn,:)=(/0,0/)
   END IF
    ! DO i=1,mn
    ! PRINT *,"ords0",ords(i,:)
    ! END DO
    !READ*,
     
   !PRINT *, "flowdirns(startcell(1), startcell(2))",flowdirns(startcell(1), startcell(2))
   upstreamcell = startcell
   upflowdirn = startcell - LSD1
   upsecondary = secondary1
   targetcell = LSD1
   counter = 1
   
   ordsDo : DO While(ords(1,1).gt.0)
    !Do i=1,m
    !PRINT *,"1: ",flowdirns(i,:)
    !END DO
    !READ*,
   ! PRINT *, "ords(1)",ords(1,:)
   ! PRINT *, "startcell,targetcell",startcell,targetcell
    
    mainif : If(targetcell(1).eq.-99) THEN
     startcell = ords(1,:)
     CALL LSDs(1, startcell, topog,m,n,lsdList)
     LSD1 = lsdList(1,:)
     secondary1 = lsdList(2,:)
     otherAnglePosition1 = lsdList(3,:)
   
     upstreamcell = startcell
     upflowdirn = startcell - LSD1
     upsecondary = secondary1
     targetcell = LSD1
     
     CALL fdirLookup(startcell-LSD1, idirn)
     flowdirns(startcell(1), startcell(2)) = idirn
     
      CALL Pos1d(ords,mn,2,startcell,rownum)
         IF (rownum(1).ne.0) THEN
           CALL endShift(ords,rownum(1),mn,2)
           ords(mn,:)=(/0,0/)
         END IF
       !DO i=1,mn
        ! PRINT *, "ords1",ords(i,:)
       !  END DO
         
     !If(topog(targetcell(1), targetcell(2)).gt. topog(startcell(1), startcell(2))) THEN
     !   Print *,"1:"
     !END IF
     
  ELSE !mainif :
  
     
     !READ*,
     ifAlreadyAsigned : If(flowdirns(targetcell(1), targetcell(2)).gt. -2) THEN 
       
       !assign flow dirn to upstreamcell
       CALL Lookupfdir(idirn, dx,dy)
       !PRINT *, "idirn",idirn
       !upstreamcell = targetcell +(/dx,dy/)
       flowdirns(upstreamcell(1), upstreamcell(2)) = idirn
       !PRINT *, "upstreamcell , targetcell",upstreamcell , targetcell
       CALL Pos1d(ords,mn,2,upstreamcell,rownum)
         IF (rownum(1).ne.0) THEN
           CALL endShift(ords,rownum(1),mn,2)
           ords(mn,:)=(/0,0/)
         END IF
         !DO i=1,mn
         !PRINT *, "ords2",ords(i,:)
         !END DO
         !READ*,
       startcell = ords(1,:)
       
       CALL LSDs(1, startcell, topog,m,n,lsdList)
       LSD1 = lsdList(1,:)
       secondary1 = lsdList(2,:)
       otherAnglePosition1 = lsdList(3,:)
       
       upstreamcell = startcell
       upflowdirn = startcell - LSD1
       upsecondary = secondary1
       targetcell = LSD1
       CALL fdirLookup(startcell-LSD1, idirn)
       flowdirns(startcell(1), startcell(2)) = idirn
       
       CALL RotateArray(ords,mn,2,-1)
       ords(mn,:)=(/0,0/)
       
       !DO i=1,mn
       !PRINT *,"ords3",ords(i,:)
       !END DO
       !READ*,
     
       ELSE !ifAlreadyAsigned 
       
        CALL LSDs(1, targetcell, topog,m,n,lsdList)
        LSD2 = lsdList(1,:)
        secondary2 = lsdList(2,:)
        otherAnglePosition2 = lsdList(3,:)
       
       If((secondary2(1) .eq. -99).OR.(LSD2(1) .eq. -99)) THEN
        flowdirns(targetcell(1), targetcell(2)) = -1
        
       CALL Pos1d(ords,mn,2,targetcell,rownum)
       IF (rownum(1).ne.0) THEN
         CALL endShift(ords,rownum(1),mn,2)
         ords(mn,:)=(/0,0/)
       END IF
       
        !DO i=1,mn
        ! PRINT *, "ords3.1",ords(i,:)
        ! END DO
         
        startcell = ords(1,:)
        CALL LSDs(order, startcell, topog,m,n,lsdList)
        LSD1 = lsdList(1,:)
        secondary1 = lsdList(2,:)
        otherAnglePosition1 = lsdList(3,:)
       
        upstreamcell = startcell
        targetcell = LSD1
        upflowdirn = startcell - LSD1
        upsecondary = secondary1
        CALL fdirLookup(startcell-LSD1, idirn)
        flowdirns(startcell(1), startcell(2)) = idirn
        
        CALL Pos1d(ords,mn,2,startcell,rownum)
         IF (rownum(1).ne.0) THEN
           CALL endShift(ords,rownum(1),mn,2)
           ords(mn,:)=(/0,0/)
         END IF
         !DO i=1,mn
         !PRINT *, "ords3.2",ords(i,:)
         !END DO
         !READ*,
         
        !If(topog(targetcell(1), targetcell(2)) > topog(startcell(1), startcell(2))) THEN
        !  Print*, "3:"
        !END IF
        
      ELSE !if ((secondary2(1) .eq. -99).OR.(LSD2(1) .eq. -99))
        
        z0 = topog(targetcell(1), targetcell(2))
        
        If(secondary2(1) .ne.-99) THEN 
         z1 = topog(secondary2(1), secondary2(2))
        ELSE
         z1 = z02
        END IF
        
        If(otherAnglePosition2(1).ne.-99) THEN
          z2 = topog(otherAnglePosition2(1),otherAnglePosition2(2))
        ELSE
          z2 = z02
        END IF
        
        !PRINT *,"z0,z1,z2",z0,z1,z2
        !READ*,
        
        If((z1 > z0).AND.(z2 > z0)) THEN
          CALL fdirLookup(targetcell-LSD2, idirn)
          flowdirns(targetcell(1), targetcell(2)) = idirn
         
          !If(topog(LSD2(1), LSD2(2)) > topog(targetcell(1), targetcell(2))) THEN
          !  Print*, "4:"
          !END IF
         
         CALL Pos1d(ords,mn,2,targetcell,rownum)
         IF (rownum(1).ne.0) THEN
           CALL endShift(ords,rownum(1),mn,2)
           ords(mn,:)=(/0,0/)
         END IF
         
        !  DO i=1,mn
        ! PRINT *, "ords4",ords(i,:)
        ! END DO
        ! READ*,
         
          
          startcell = LSD2
          CALL LSDs(1, startcell, topog,m,n,lsdList)
          LSD1 = lsdList(1,:)
          secondary1 = lsdList(2,:)
          otherAnglePosition1 = lsdList(3,:)
         
         If(LSD1(1) .ne.-99) THEN
          CALL fdirLookup(startcell-LSD1, idirn)
          flowdirns(startcell(1), startcell(2)) = idirn
          
         ! If(topog(LSD1(1), LSD1(2)) > topog(startcell(1), startcell(2))) THEN
         !   Print*,"5:"
         ! END IF
          
          CALL Pos1d(ords,mn,2,startcell,rownum)
          IF (rownum(1).ne.0) THEN
            CALL endShift(ords,rownum(1),mn,2)
            ords(mn,:)=(/0,0/)
          END IF
           ! DO i=1,mn
           !PRINT *, "ords5",ords(i,:)
           !END DO
           !READ*,
          
          
          upflowdirn = startcell - LSD1
          upsecondary = secondary1
         
         ELSE !if (LSD1(1) .ne.-99)
          
          flowdirns(startcell(1), startcell(2)) = -1
          
          CALL Pos1d(ords,mn,2,startcell,rownum)
          IF (rownum(1).ne.0) THEN
            CALL endShift(ords,rownum(1),mn,2)
            ords(mn,:)=(/0,0/)
          END IF
          !DO i=1,mn
         !PRINT *, "ords6",ords(i,:)
         !END DO
         !READ*,
           
          upflowdirn = (/-99, -99/)
          upsecondary = (/-99, -99/)
        END IF !(LSD1(1) .ne.-99) 
         
         targetcell = LSD1
         CALL LSDs(1, targetcell, topog,m,n,lsdList)
         LSD2 = lsdList(1,:)
         secondary2 = lsdList(2,:)
         otherAnglePosition2 = lsdList(3,:)
         upstreamcell = startcell
        ELSE !if ((z1 > z0).AND.(z2 > z0))
         dirn11 = upstreamcell - targetcell
         dirn21 = targetcell - LSD2
         test1 = ((dirn11(1) .eq. dirn21(1)).AND.(dirn11(2) .eq. dirn21(2)))
         dirn12 = upstreamcell - upsecondary
         dirn22 = targetcell - secondary2
         test2 = ((dirn12(1) .eq. dirn22(1)).AND.(dirn12(2) .eq. dirn22(2)))
         slp1 = (topog(startcell(1), startcell(2)) - & 
                    topog(LSD2(1), LSD2(2)))/Sqrt(DBLE(((startcell(1) - & 
                 LSD2(1))**2 + (startcell(2) - LSD2(2))**2)))
                 
         slp2 = (topog(startcell(1), startcell(2)) - & 
                    topog(secondary2(1), secondary2(2)))/&
                    Sqrt(DBLE(((startcell(1) - secondary2(1))**2 + (startcell(2) - & 
                 secondary2(2))**2)))
         
         test3 = (slp1 < slp2)
         
         If((test1.AND.test2).AND.test3) THEN
          !PRINT *, "test1,2,and 3"
          If(topog(secondary2(1), secondary2(2)).le. topog(targetcell(1), targetcell(2))) THEN
                CALL fdirLookup(targetcell-secondary2, idirn)
                flowdirns(targetcell(1), targetcell(2)) = idirn
          ELSE
            CALL fdirLookup(targetcell-LSD2, idirn)
            flowdirns(targetcell(1), targetcell(2)) = idirn
           
           !If(topog(LSD2(1), LSD2(2)) > topog(targetcell(1), targetcell(2))) THEN
           !   Print*, "6b:"
           !END IF
           
          END IF !(topog(secondary2(1), secondary2(2)).le. topog(targetcell(1), targetcell(2)))
          
          CALL Pos1d(ords,mn,2,targetcell,rownum)
          IF (rownum(1).ne.0) THEN
            CALL endShift(ords,rownum(1),mn,2)
            ords(mn,:)=(/0,0/)
          END IF
          ! DO i=1,mn
         !PRINT *, "ords7",ords(i,:)
         !END DO
         !READ*,
          
          upstreamcell = targetcell
          upflowdirn = targetcell - secondary2
          upsecondary = secondary2
          targetcell = secondary2
        ELSE !if ((test1.AND.test2).AND.test3)
           CALL fdirLookup(targetcell-LSD2, idirn)
           flowdirns(targetcell(1), targetcell(2)) = idirn
          
          If(topog(LSD2(1), LSD2(2)) > topog(targetcell(1), targetcell(2))) THEN 
          ! Print*, "7:"
          END IF
          
          CALL Pos1d(ords,mn,2,targetcell,rownum)
          IF (rownum(1).ne.0) THEN
            CALL endShift(ords,rownum(1), mn,2)
            ords(mn,:)=(/0,0/)
          END IF
          
          ! DO i=1,mn
         !PRINT *, "ords8",ords(i,:)
         !END DO
         !READ*,
         
          
          upstreamcell = targetcell
          upflowdirn = targetcell - LSD2
          upsecondary = secondary2
          targetcell = LSD2
          END IF !((test1.AND.test2).AND.test3)
        END IF  !((z1 > z0).AND.(z2 > z0))
      END IF !((secondary2(1) .eq. -99).OR.(LSD2(1) .eq. -99))
      END IF ifAlreadyAsigned 
     
    END IF  mainif
    
       !DO i=1,mn
       !  PRINT *, "ords end",ords(i,:)
       !  END DO
       !  READ*,
         
 END DO ordsDo
 
  !for perodic boundary conditions set lower
  ! boudning cells to have a flow dirn of 
  !-3 i.e. assuming topography slopes downward in the direction of 
  !lower i and j
  If (isPeriodic.eq.1) THEN
  DO i=1,m
    flowdirns(i, 1) = -3
  END DO
  DO j=1,n
    flowdirns(1,j) = -3
  END DO
  END IF 
   
END SUBROUTINE GD8   
!NewGD8
Subroutine NewGD8(topog,lakes,flowdirns, m,n)
   IMPLICIT NONE
   
   INTEGER,INTENT(IN) :: m,n
   INTEGER, DIMENSION(m,n), INTENT(IN) :: lakes
   REAL*8, DIMENSION(m,n), INTENT(IN) :: topog
   INTEGER, DIMENSION(m,n), INTENT(out) :: flowdirns
   
   REAL*8 :: zdefault, z0, z1, z2, slp1, slp2
   INTEGER, DIMENSION(m,n) :: isAssigned
   INTEGER :: mn, num2Assign, order, i,j,k, idirn
   LOGICAL :: isAss, isOutflow, reinit, test1, test2, test3
   INTEGER, DIMENSION(SIZE(topog),2) :: ords
   INTEGER, DIMENSION(SIZE(topog)) :: rownum
   
   INTEGER, DIMENSION(2) :: startcell, targetcell
   INTEGER, DIMENSION(2) :: LSD1, secondary1, otherAnglePosition1 
   INTEGER, DIMENSION(2) :: LSD2, secondary2, otherAnglePosition2 
   INTEGER, DIMENSION(2) :: upstreamcell, upflowdirn, upsecondary
   INTEGER, DIMENSION(2) ::  dirn11, dirn21, dirn12, dirn22
   INTEGER, DIMENSION(9,2) :: neighbs
   INTEGER, DIMENSION(3,2) :: lsdList
    
   flowdirns = -2
   WHERE (lakes>0)
     isAssigned =1
     flowdirns = -1
   ELSEWHERE
     isAssigned =0
   END WHERE
   
   mn = m*n
   
   zdefault = 0.d0
   order = 1
   CALL makeOrds(topog,ords, m,n)  !orders 2D positions of topography
   
   startcell=minloc(topog)
   PRINT *,"minloc",startcell
   flowdirns(startcell(1),startcell(2))=-1  !make lowest point an outflow point
   CALL Pos1d(ords,mn,2,startcell,rownum)
   IF (rownum(1).ne.0) THEN
     CALL endShift(ords,rownum(1),mn,2)
     ords(mn,:)=(/0,0/)
   END IF
   isAssigned(startcell(1),startcell(2))=1
   
   !DO i=1,mn
   !PRINT *,"ords",ords(i,:), " topog",topog(ords(i,1),ords(i,2))
   !END DO
   isAss=.TRUE.   !dummy assignment to begin loop
   isOutflow=.TRUE.   !dummy assignment to begin loop
   reinit = .TRUE.
   startcell=(/1,1/)  !dummy assignment to begin loop
   startcell = ords(1,:)
   10 num2AssignGT0 : DO WHILE(startcell(1).ne.0)
    
    !DO i=1,m
    !PRINT *,"1: flowdirns-",flowdirns(i,:)
    !END DO
    
    !PRINT *, "2: isAss, isOutflow",isAss,isOutflow
    !READ*,
    
    toReinit : IF(isAss.OR.isOutflow) THEN
       !This bit reinitialises all variables if beginning a new flow pathway
        !PRINT *,"3: reinit", reinit
        IF (reinit) THEN
         startcell = ords(1,:)
       END IF
       PRINT *,"4: startcell", startcell
       IF (startcell(1).ne.0) THEN
         CALL Pos1d(ords,mn,2,startcell,rownum)
         IF (rownum(1).ne.0) THEN
           CALL endShift(ords,rownum(1),mn,2)
           ords(mn,:)=(/0,0/)
         END IF
       
       !READ*,
       !DO i=1,mn
       ! PRINT *,"ords", ords(i,:)
       !END DO
       !READ*,
       
        CALL LSDs(1, startcell, topog,m,n,lsdList)
        LSD1 = lsdList(1,:)
        secondary1 = lsdList(2,:)
        otherAnglePosition1 = lsdList(3,:)
        PRINT *, "5: LSD1, secondary1, otherAnglePosition1/", LSD1, "/",secondary1, "/",otherAnglePosition1
        CALL fdirLookup(startcell-LSD1, idirn)
        flowdirns(startcell(1), startcell(2)) = idirn
        isAssigned(startcell(1),startcell(2))=1
        upstreamcell = startcell
        upflowdirn = startcell - LSD1
        upsecondary = secondary1
        targetcell = LSD1
         PRINT *,"6: targetcell", targetcell
         
        CALL Pos1d(ords,mn,2,targetcell,rownum)
        IF (rownum(1).ne.0) THEN
          CALL endShift(ords,rownum(1),mn,2)
          ords(mn,:)=(/0,0/)
        END IF
        
        !READ*,
        !DO i=1,mn
        !PRINT *,"ords", ords(i,:)
        !END DO
        !READ*,
        
        PRINT *,"7: isAssigned(targetcell(1),targetcell(2)).eq.1",isAssigned(targetcell(1),targetcell(2)).eq.1
        ifassigned1 : IF (isAssigned(targetcell(1),targetcell(2)).eq.1) THEN
            isAss=.TRUE.
            isOutflow=.FALSE.
            reinit = .TRUE.
            !change so this so that if all neighbours uphill and or boundary cells then isOutflow=.TRUE.
            CALL Neighbours(1,targetcell, (/m,n/),neighbs)
            DO i=1,9
              If (i.ne.5) THEN
                IF (neighbs(i,1).ne.-99) THEN
                  isOutflow = (isOutflow.and.(topog(neighbs(i,1),neighbs(i,2))>topog(targetcell(1),targetcell(2))))
                END IF
              END IF
            END DO
            
             PRINT *,"8: isOutflow",isOutflow
             
            IF (isOutflow) THEN
              isAssigned(targetcell(1),targetcell(2))=1
              CALL Pos1d(ords,mn,2,targetcell,rownum)
              IF (rownum(1).ne.0) THEN
                CALL endShift(ords,rownum(1), mn,2)
                ords(mn,:)=(/0,0/)
              END IF
              reinit = .TRUE.
            END IF
             PRINT *,"9: reinit",reinit
          ELSE !if isAssigned
            isAssigned(targetcell(1),targetcell(2))=1
            isAss=.FALSE.  !set condition to move program into fdirAss  do while loop
            isOutflow=.FALSE.
            reinit = .FALSE.
          END IF  ifassigned1
          PRINT *,  "10: isAss, isOutflow, reinit", isAss, isOutflow, reinit
          !READ*,
      ELSE !if startcell(1).ne.0
        RETURN
      END IF
    ELSE  ! toReinit
      
      PRINT *,  "11: isAss, isOutflow, reinit",  isAss, isOutflow, reinit
      PRINT *, "12: (.NOT. isAss) .AND. (.NOT. isOutflow) .AND. (.NOT. reinit)", &
                    (.NOT. isAss) .AND. (.NOT. isOutflow) .AND. (.NOT. reinit)
      !READ*,
      Do WHILE ( (.NOT. isAss) .AND. (.NOT. isOutflow) .AND. (.NOT. reinit)) 
        
        DO i=1,m
          PRINT *,"12:0: flowdirns-",flowdirns(i,:)
        END DO
        !READ*,
        PRINT *,"13.0: startcell, targetcell",startcell,targetcell
        CALL LSDs(1, targetcell, topog,m,n,lsdList)
        PRINT *,"13.1 lsdList", lsdList
        !READ*,
        LSD2 = lsdList(1,:)
        secondary2 = lsdList(2,:)
        otherAnglePosition2 = lsdList(3,:)  
        PRINT *, "13: LSD2, secondary2, otherAnglePosition2,targetcell/", LSD2, "/",secondary2, "/",otherAnglePosition2,"/" &
           ,targetcell
        boundarycheck : IF (LSD2(1).eq.-99) THEN
          flowdirns(targetcell(1), targetcell(2)) = -1
          isAssigned(targetcell(1),targetcell(2))=1
          CALL Pos1d(ords,mn,2,targetcell,rownum)
          IF (rownum(1).ne.0) THEN
           CALL endShift(ords,rownum(1),mn,2)
           ords(mn,:)=(/0,0/)
          END IF
          reinit = .TRUE.
          isOutFlow = .TRUE.
          Print *, "Gone to 10, 1"
          GOTO 10
        EnD IF boundarycheck
        
        
        !get local elevations
        z0 = topog(LSD2(1), LSD2(2))
        If(secondary2(1) .ne.-99) THEN 
         z1 = topog(secondary2(1), secondary2(2))
        ELSE
         z1 = zdefault
        END IF
        If(otherAnglePosition2(1).ne.-99) THEN
          z2 = topog(otherAnglePosition2(1),otherAnglePosition2(2))
        ELSE
          z2 = zdefault
        END IF
         PRINT *, "14: z0, z1, z2",z0, z1, z2 
         PRINT *, "15: (z1 > z0).AND.(z2 > z0)",(z1 > z0).AND.(z2 > z0) 
        !valley check
         valleyCheck : If((z1 > z0).AND.(z2 > z0)) THEN
         PRINT *, "Valley"
         !it is a valley so assign flow dirn to LSD2 and reset startcell
          CALL fdirLookup(targetcell-LSD2, idirn)
          flowdirns(targetcell(1), targetcell(2)) = idirn
          isAssigned(targetcell(1),targetcell(2))=1
          CALL Pos1d(ords,mn,2,targetcell,rownum)
          IF (rownum(1).ne.0) THEN
           CALL endShift(ords,rownum(1),mn,2)
           ords(mn,:)=(/0,0/)
          END IF
          startcell = LSD2
           PRINT *, "16: startcell", startcell
          isAss = .TRUE.  !dummy call to reinitialise
          PRINT *, "17: isAss", isAss 
          !programming comment: need to make sure if this condition 
          !met then next evaluation passed to end of this do while loop
          ELSE !if valleyCheck
          PRINT *,"Not a valley"
          !begin three tests for flow dirn adjustment
          
            dirn11 = upstreamcell - targetcell
            dirn21 = targetcell - LSD2
            dirn12 = upstreamcell - upsecondary
            dirn22 = targetcell - secondary2
            
            slp1 = (topog(startcell(1), startcell(2)) - & 
                    topog(LSD2(1), LSD2(2)))/Sqrt(DBLE(((startcell(1) - & 
                   LSD2(1))**2 + (startcell(2) - LSD2(2))**2)))
            slp2 = (topog(startcell(1), startcell(2)) - & 
                    topog(secondary2(1), secondary2(2)))/&
                    Sqrt(DBLE(((startcell(1) - secondary2(1))**2 + (startcell(2) - & 
                 secondary2(2))**2)))
                 
            test1 = ((dirn11(1) .eq. dirn21(1)).AND.(dirn11(2) .eq. dirn21(2)))
            test2 = ((dirn12(1) .eq. dirn22(1)).AND.(dirn12(2) .eq. dirn22(2)))
            test3 = (slp1 < slp2)
            PRINT *, "18: test1,test2,test3", test1,test2,test3
            PRINT *,"*****************************************"
            test123 : IF ((test1.and.test2).and.test3) THEN
            
              If(topog(secondary2(1), secondary2(2)).le. topog(targetcell(1), targetcell(2))) THEN
                CALL fdirLookup(targetcell-secondary2, idirn)
                flowdirns(targetcell(1), targetcell(2)) = idirn
                isAssigned(targetcell(1),targetcell(2))=1
                PRINT *,"-----------------------------------------------------"
              ELSE
                CALL fdirLookup(targetcell-LSD2, idirn)
                flowdirns(targetcell(1), targetcell(2)) = idirn
                isAssigned(targetcell(1),targetcell(2))=1
                PRINT *,"+++++++++++++++++++++++++++++++++++++++++"
              END IF
              !READ*,
              !clean up ords
              CALL Pos1d(ords,mn,2,targetcell,rownum)
              IF (rownum(1).ne.0) THEN
                CALL endShift(ords,rownum(1),mn,2)
                ords(mn,:)=(/0,0/)
              END IF
              upstreamcell = targetcell
              upflowdirn = targetcell - secondary2
              upsecondary = secondary2
              targetcell = secondary2
              PRINT *, "19:  startcell, targetcell",  startcell, targetcell
            ELSE !test123
              CALL fdirLookup(targetcell-LSD2, idirn)
              flowdirns(targetcell(1), targetcell(2)) = idirn
              isAssigned(targetcell(1),targetcell(2))=1
              CALL Pos1d(ords,mn,2,targetcell,rownum)
              IF (rownum(1).ne.0) THEN
                CALL endShift(ords,rownum(1), mn,2)
                ords(mn,:)=(/0,0/)
              END IF
              upstreamcell = targetcell
              upflowdirn = targetcell - LSD2
              upsecondary = secondary2
              targetcell = LSD2
              PRINT *, "20:  startcell, targetcell",  startcell, targetcell
          END IF test123 !((test1.AND.test2).AND.test3)
          END IF valleyCheck
          
          PRINT *, "21:  (isAssigned(targetcell(1),targetcell(2)).eq.1) ",  (isAssigned(targetcell(1),targetcell(2)).eq.1) 
          ifassigned2 : IF (isAssigned(targetcell(1),targetcell(2)).eq.1) THEN
            isAss=.TRUE.
            isOutflow=.TRUE.
            reinit = .FALSE.
            !change so this so that if all neighbours uphill and or boundary cells then isOutflow=.TRUE.
            CALL Neighbours(1,targetcell, (/m,n/),neighbs)
            DO i=1,9
              !PRINT *,"21.00 neighbs(i)",neighbs(i,:)
              If (i.ne.5) THEN
                IF (neighbs(i,1).ne.-99) THEN
                  !PRINT *,"21.01: isOutflow",isOutflow 
                   !PRINT *,"21.01:topog1,topog2", topog(neighbs(i,1),neighbs(i,2)),topog(targetcell(1),targetcell(2))
                  isOutflow = (isOutflow.and.(topog(neighbs(i,1),neighbs(i,2))>topog(targetcell(1),targetcell(2))))
                END IF
              END IF
            END DO
            PRINT *, "21.0: isOutflow",isOutflow
            IF (isOutflow) THEN
              isAssigned(targetcell(1),targetcell(2))=1
              CALL Pos1d(ords,mn,2,targetcell,rownum)
              IF (rownum(1).ne.0) THEN
                CALL endShift(ords,rownum(1), mn,2)
                ords(mn,:)=(/0,0/)
              END IF
              reinit = .TRUE.
            END IF
          ELSE !if isAssigned
            isAssigned(targetcell(1),targetcell(2))=1
            isAss=.FALSE.  !set condition to move program into fdirAss  do while loop
            isOutflow=.FALSE.
            reinit = .FALSE.
          END IF ifassigned2 !isAssigned
          PRINT *, "22:  isAss, isOutflow, reinit",isAss, isOutflow, reinit
      END DO
    END IF toReinit
    PRINT *, "23:  isAss, isOutflow, reinit",isAss, isOutflow, reinit
   END DO num2AssignGT0
    
END SUBROUTINE NewGD8
!fdirLookup
SUBROUTINE fdirLookup(dirnxy, idirn)
  IMPLICIT NONE
  
  INTEGER, DIMENSION(2), INTENT(IN) :: dirnxy
  INTEGER, INTENT(OUT) :: idirn
  
  IF ((dirnxy(1).eq.-1).AND.(dirnxy(2).eq.-1)) THEN 
    idirn =  1
  ELSEIF ((dirnxy(1).eq.-1).AND.(dirnxy(2).eq.0)) THEN 
    idirn = 2
  ELSEIF ((dirnxy(1).eq.-1).AND.(dirnxy(2).eq.1)) THEN
    idirn = 3
  ELSEIF ((dirnxy(1).eq.0).AND.(dirnxy(2).eq.1)) THEN
    idirn = 4
  ELSEIF ((dirnxy(1).eq.1).AND.(dirnxy(2).eq.1)) THEN
    idirn = 5
  ELSEIF ((dirnxy(1).eq.1).AND.(dirnxy(2).eq.0)) THEN
    idirn = 6
  ELSEIF ((dirnxy(1).eq.1).AND.(dirnxy(2).eq.-1)) THEN
    idirn = 7
  ELSEIF ((dirnxy(1).eq.0).AND.(dirnxy(2).eq.-1)) THEN
    idirn = 8
  ELSEIF ((dirnxy(1).eq.0).AND.(dirnxy(2).eq.0)) THEN
    idirn = 9
  ELSE
    idirn = -1  
  END IF
  !PRINT *,"idrn in fdirLookup", idirn
  END SUBROUTINE fdirLookup
!qsortd  
SUBROUTINE qsortd(x,ind,n,incdec)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-12-18  Time: 11:55:47

IMPLICIT NONE
!INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL*8, INTENT(IN)  :: x(n)
INTEGER, INTENT(OUT)   :: ind(n)
INTEGER, INTENT(IN)    :: n
INTEGER, INTENT(IN) :: incdec


!***************************************************************************

!                                                         ROBERT RENKA
!                                                 OAK RIDGE NATL. LAB.

!   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL (dp)
! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.

! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

!                      X - VECTOR OF LENGTH N TO BE SORTED.

!                    IND - VECTOR OF LENGTH >= N.
!                    incdec - integer, if positive then ind is returned so values in decreasing order 

! N AND X ARE NOT ALTERED BY THIS ROUTINE.

! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
!                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
!                          X IS DEFINED BY Y(I) = X(IND(I)).

!*********************************************************************

! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

!*********************************************************************

INTEGER   :: iu(21), il(21)
INTEGER   :: m, i, j, k, l, ij, it, itt, indx
REAL*8      :: r
REAL*8    :: t, temp

! LOCAL PARAMETERS -

! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
!            INDICES OF PORTIONS OF THE ARRAY X
! M =      INDEX FOR IU AND IL
! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =    INDICES IN THE RANGE I,...,J
! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =   TEMPORARY INDEX FOR X
! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =      CENTRAL ELEMENT OF X

IF (n <= 0) RETURN

! INITIALIZE IND, M, I, J, AND R

DO  i = 1, n
  ind(i) = i
END DO
m = 1
i = 1
j = n
r = .375

! TOP OF LOOP

20 IF (i >= j) GO TO 70
IF (r <= .5898437) THEN
  r = r + .0390625
ELSE
  r = r - .21875
END IF

! INITIALIZE K

30 k = i

! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

ij = i + r*(j-i)
it = ind(ij)
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) > t) THEN
  ind(ij) = indx
  ind(i) = it
  it = indx
  t = x(it)
END IF

! INITIALIZE L

l = j

! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!   INTERCHANGE IT WITH T

indx = ind(j)
IF (x(indx) >= t) GO TO 50
ind(ij) = indx
ind(j) = it
it = indx
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) <= t) GO TO 50
ind(ij) = indx
ind(i) = it
it = indx
t = x(it)
GO TO 50

! INTERCHANGE ELEMENTS K AND L

40 itt = ind(l)
ind(l) = ind(k)
ind(k) = itt

! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!   NOT LARGER THAN T

50 l = l - 1
indx = ind(l)
IF (x(indx) > t) GO TO 50

! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60 k = k + 1
indx = ind(k)
IF (x(indx) < t) GO TO 60

! IF K <= L, INTERCHANGE ELEMENTS K AND L

IF (k <= l) GO TO 40

! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!   ARRAY YET TO BE SORTED

IF (l-i > j-k) THEN
  il(m) = i
  iu(m) = l
  i = k
  m = m + 1
  GO TO 80
END IF

il(m) = k
iu(m) = j
j = l
m = m + 1
GO TO 80

! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70 m = m - 1
IF (m == 0) GO TO 110
i = il(m)
j = iu(m)

80 IF (j-i >= 11) GO TO 30
IF (i == 1) GO TO 20
i = i - 1

! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90 i = i + 1
IF (i == j) GO TO 70
indx = ind(i+1)
t = x(indx)
it = indx
indx = ind(i)
IF (x(indx) <= t) GO TO 90
k = i

100 ind(k+1) = ind(k)
k = k - 1
indx = ind(k)
IF (t < x(indx)) GO TO 100

ind(k+1) = it
GO TO 90


110 IF (incdec.gt.0) THEN 
  DO i=1,size(x)/2
  temp=ind(size(x)-i+1)
  ind(size(x)-i+1)=ind(i)
  ind(i)=temp
  END DO
END IF

RETURN
  

END SUBROUTINE qsortd
!endShift
SUBROUTINE endShift(arr,rownum,mn,n)
  !puts a row at rownum to the last row of the array arr
  IMPLICIT NONE
  
  INTEGER, DIMENSION(mn,n), INTENT(INOUT) :: arr
  INTEGER, INTENT(IN) :: rownum, mn, n
  
  IF(rownum.ne.mn) THEN
    arr(rownum:mn-1,:)= arr(rownum+1:mn,:)
  END IF
  
END SUBROUTINE endShift
!makeOrds
SUBROUTINE makeOrds(topog,ords, m,n)
  IMPLICIT NONE
  
  REAL*8, DIMENSION(m,n), INTENT(IN) :: topog
  INTEGER, DIMENSION(SIZE(topog),2), INTENT(OUT) :: ords
  INTEGER, INTENT(IN) :: m, n
  INTEGER :: i, j, k, mn, x,y
  INTEGER, DIMENSION(SIZE(topog)) :: ind2
  REAL*8, DIMENSION(SIZE(topog)) :: oneDtopog
  
  mn = m*n
  k=1
  DO i=1,m
    DO j=1,n
      oneDtopog(k)= topog(i,j)
      k=k+1
    END DO
  END DO
  
  ind2=(/(j,j=1,mn)/)
  CALL qsortd(oneDtopog, ind2, mn,1)
  DO i=1,mn
    If(Mod(ind2(i),n).eq.0) THEN
      x=0
      y=n
    ELSE
      x=1
      y=0
    END IF
    ords(i,:)=(/Floor(Real(ind2(i))/Real(n))+x , ind2(i) - Floor(Real(ind2(i))/Real(n))*n+y /)
    END DO
END SUBROUTINE makeOrds
!VegChange
SUBROUTINE VegChange(veg,m,n,vegParams, store, actualET,isEmerge,isGrow)
 IMPLICIT NONE
 
 INTEGER, INTENT(IN) :: m, n, isGrow,isEmerge
 REAL*8, DIMENSION(5), INTENT(IN) :: vegParams
 INTEGER, DIMENSION(m,n), INTENT(INOUT) :: veg, store, actualET
 
 REAL*8 :: pc !prob of collonisation of bare soil
 INTEGER :: vegmax, i, j, sEmerge, etPersist, usePCFlag
 !REAL*8 :: rnd
 
 vegmax = int(vegParams(1)) ! parameter set to define maximum biomass was 9
 sEmerge = int(vegParams(2))
 etPersist = int(vegParams(3))
 pc = vegParams(4)
 usePCFlag = int(vegParams(5))  !flag denotes whether to use random collonisation pc or storage based sEmerge
 
 Do i=1,m
  DO j=1,n
    If(veg(i, j).eq.0) THEN 
     If (isEmerge.gt.0) THEN
      IF (usePCFlag.gt.0) THEN
        If(store(i, j) > sEmerge) THEN
          store(i, j) = store(i, j) - sEmerge
          actualET(i,j) = actualET(i,j)+sEmerge
          veg(i, j) = 1
        !ELSE
          !store(i, j) = 0
        END IF
      ELSE
        CALL random_number(rnd)
        If(rnd < pc) THEN
          veg(i, j) = 1
        END IF
      END IF !end if (usePCFlag.lt.0)
     END IF
    ELSE
    If (isGrow.gt.0) THEN
      If(real(actualET(i, j)) > REAL(eTPersist) ) THEN
        veg(i, j) = min(veg(i, j) + 1,vegmax)
      ELSE
        veg(i, j) = veg(i, j) - 1  !was veg(i, j) - 1
      END IF
    END IF
    END IF
   END DO
 END DO
  
END SUBROUTINE VegChange

END PROGRAM Sensitivity


