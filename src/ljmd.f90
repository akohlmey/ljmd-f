! 
! simple lennard-jones potential MD code with velocity verlet.
! units: Length=Angstrom, Mass=amu, Energy=kcal
!
! optimized f95 version using cell lists
!

MODULE kinds
  IMPLICIT NONE
  INTEGER, PARAMETER :: dbl = selected_real_kind(14,200)  ! double precision floating point
  INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)    ! single precision floating point
  INTEGER, PARAMETER :: sln = 200                         ! length of I/O input line
  PRIVATE
  PUBLIC :: sgl, dbl, sln
END MODULE kinds

MODULE physconst
  USE kinds
  IMPLICIT NONE
  REAL(kind=dbl), PARAMETER :: kboltz =    0.0019872067_dbl   ! boltzman constant in kcal/mol/K
  REAL(kind=dbl), PARAMETER :: mvsq2e = 2390.05736153349_dbl  ! m*v^2 in kcal/mol
  PRIVATE
  PUBLIC :: kboltz, mvsq2e
END MODULE physconst

! module to hold the complete system information 
MODULE mdsys
  USE kinds
  IMPLICIT NONE
  INTEGER :: natoms,nfi,nsteps,nthreads
  REAL(kind=dbl) dt, mass, epsilon, sigma, box, rcut
  REAL(kind=dbl) ekin, epot, temp
  REAL(kind=dbl), POINTER, DIMENSION (:,:) :: pos, vel
  REAL(kind=dbl), POINTER, DIMENSION (:,:,:) :: frc
END MODULE mdsys

MODULE utils
  USE kinds
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pbc

CONTAINS
   
! helper function: apply minimum image convention 
  FUNCTION pbc(x, boxby2, box)
    REAL(kind=dbl), INTENT(IN)  :: x, boxby2, box
    REAL(kind=dbl) :: pbc

    pbc = x
    DO WHILE(pbc > boxby2)
       pbc = pbc - box
    END DO
    DO WHILE(pbc < -boxby2)
       pbc = pbc + box
    END DO
  END FUNCTION pbc
END MODULE utils

MODULE io
  USE kinds
  IMPLICIT NONE
  PRIVATE 
  INTEGER, PARAMETER :: stdin=5, stdout=6, log=30, xyz=31
  PUBLIC :: ioopen, ioclose, output, stdin, stdout, getline

CONTAINS
  SUBROUTINE getline(chan, line)
    INTEGER, INTENT(IN) :: chan
    CHARACTER(LEN=sln), INTENT(OUT) :: line
    INTEGER :: idx, i

    READ(CHAN, '(A)') line
    ! delete comment
    idx=INDEX(line,'#')
    IF (idx > 0) THEN
       DO i=idx,sln
          line(i:i) = ' '
       END DO
    END IF
  END SUBROUTINE getline

  SUBROUTINE ioopen(logname, xyzname)
    CHARACTER(LEN=sln) :: logname, xyzname
    OPEN(UNIT=log, FILE=TRIM(logname), STATUS='UNKNOWN', FORM='FORMATTED')
    OPEN(UNIT=xyz, FILE=TRIM(xyzname), STATUS='UNKNOWN', FORM='FORMATTED')
  END SUBROUTINE ioopen
  
  SUBROUTINE ioclose
    CLOSE(UNIT=log)
    CLOSE(UNIT=xyz)
  END SUBROUTINE ioclose
  
  ! append data to output.
  SUBROUTINE output
    USE mdsys
    IMPLICIT NONE
    INTEGER :: i
    WRITE(log, '(I8,1X,F20.8,1X,F20.8,1X,F20.8,1X,F20.8)') &
         nfi, temp, ekin, epot, ekin+epot
    WRITE(stdout, '(I8,1X,F20.8,1X,F20.8,1X,F20.8,1X,F20.8)') &
         nfi, temp, ekin, epot, ekin+epot
    WRITE(xyz, '(I8)') natoms
    WRITE(xyz, '(A,I8,1X,A,F20.8)') 'nfi=', nfi, 'etot=', ekin+epot
    DO i=1, natoms
       WRITE(xyz, '(A,1X,F20.8,1X,F20.8,1X,F20.8)') &
            'Ar ', pos(i,1), pos(i,2), pos(i,3)
    END DO
  END SUBROUTINE output
END MODULE io

MODULE cell
  USE kinds
  IMPLICIT NONE
  INTEGER :: npair, ncell, ngrid, nidx
  INTEGER, POINTER :: npercell(:), clist(:,:), plist(:)
  REAL(KIND=dbl), PARAMETER :: cellrat = 2.0_dbl
  INTEGER, PARAMETER :: cellfreq = 4
  REAL(KIND=dbl) :: dcell
  PRIVATE
  PUBLIC :: ncell, npair, clist, plist, npercell
  PUBLIC :: mkcell, updcell, rmcell
  
CONTAINS

  SUBROUTINE mkcell
    USE io
    USE mdsys
    USE utils
    REAL(kind=dbl) :: boxby2, boxoffs, x1, y1, z1, x2, y2, z2, rx, ry, rz
    INTEGER :: i, j, k
    
        
    ngrid   = FLOOR(cellrat * box / rcut)
    ncell   = ngrid * ngrid * ngrid
    dcell   = box / ngrid
    boxby2  = 0.5_dbl * box
    boxoffs = boxby2 - 0.5_dbl*dcell
    nidx = 2*natoms / ncell + 2
    nidx = ((nidx/2) + 1) * 2
        
    ! allocate cell list storage 
    ALLOCATE(npercell(ncell), clist(ncell,nidx), plist(2*ncell*ncell))

    ! build cell pair list, assuming newtons 3rd law. */
    npair = 0
    DO i=0, ncell-2
       k  = i/ngrid/ngrid
       x1 = k*dcell - boxoffs
       y1 = ((i-(k*ngrid*ngrid))/ngrid)*dcell - boxoffs
       z1 = MOD(i,ngrid)*dcell - boxoffs
       
       DO j=i+1, ncell-1
          k  = j/ngrid/ngrid
          x2 = k*dcell - boxoffs
          y2 = ((j-(k*ngrid*ngrid))/ngrid)*dcell - boxoffs
          z2 = MOD(j,ngrid)*dcell - boxoffs
          
          rx=pbc(x1-x2, boxby2, box)
          ry=pbc(y1-y2, boxby2, box)
          rz=pbc(z1-z2, boxby2, box)

          ! check for cells on a line that are too far apart
          IF (ABS(rx) > rcut+dcell) CYCLE
          IF (ABS(ry) > rcut+dcell) CYCLE
          if (ABS(rz) > rcut+dcell) CYCLE

          ! check for cell within a plane that are too far apart.
          IF (SQRT(rx*rx+ry*ry) > (rcut+SQRT(2.0_dbl)*dcell)) CYCLE
          IF (SQRT(rx*rx+rz*rz) > (rcut+SQRT(2.0_dbl)*dcell)) CYCLE
          IF (SQRT(ry*ry+rz*rz) > (rcut+SQRT(2.0_dbl)*dcell)) CYCLE

          ! other cells that are too far apart 
          IF (DSQRT(rx*rx + ry*ry + rz*rz) > (DSQRT(3.0_dbl)*dcell+rcut)) CYCLE
                
          ! cells are close enough. add to list.
          npair=npair+1
          plist(2*npair-1) = i+1
          plist(2*npair  ) = j+1
       END DO
    END DO

    WRITE(stdout,'(A,I2,"x",I2,"x",I2,"=",I6,A,I8,"/",I12,A,I4,A)') 'Cell list has ', &
         ngrid, ngrid, ngrid, ncell, ' cells with ', npair, ncell*(ncell-1)/2, &
         ' pairs and ', nidx, ' atoms/celllist'
  END SUBROUTINE mkcell

  ! update cell list contents
  SUBROUTINE updcell
    USE io
    USE mdsys
    USE utils

    INTEGER :: i, j, k, m, n, midx, idx
    REAL(kind=dbl) :: boxby2

    IF (MOD(nfi,cellfreq) > 0) RETURN

    npercell(:) = 0
    midx = 0

    DO i=1, natoms
       
       ! compute atom position in cell grid and its array index
       k=FLOOR((pbc(pos(i,1), boxby2, box)+boxby2)/dcell)
       m=FLOOR((pbc(pos(i,2), boxby2, box)+boxby2)/dcell)
       n=FLOOR((pbc(pos(i,3), boxby2, box)+boxby2)/dcell)
       j = ngrid*ngrid*k+ngrid*m+n+1
       idx = npercell(j) + 1
       clist(j,idx) = i
       npercell(j) = idx
       IF (idx > midx) midx=idx
    END DO
    IF (midx > nidx) THEN
       WRITE(stdout, '(A,I8,"/",I8,A)') 'Overflow in cell list: ', midx, nidx, ' atoms/cell.'
       STOP 'Fatal error'
    END IF
  END SUBROUTINE updcell

  SUBROUTINE rmcell
    ncell = 0
    nidx  = 0
    ngrid = 0
    npair = 0
    DEALLOCATE(npercell, clist, plist)
  END SUBROUTINE rmcell
END MODULE cell

! compute kinetic energy
SUBROUTINE getekin
  USE kinds
  USE mdsys, ONLY: natoms, mass, temp, ekin, vel
  USE physconst
  IMPLICIT NONE

  INTEGER :: i

  ekin = 0.0_dbl
  DO i=1,natoms
     ekin = ekin + 0.5_dbl * mvsq2e * mass * dot_product(vel(i,:),vel(i,:))
  END DO
  temp = 2.0_dbl * ekin/(3.0_dbl*DBLE(natoms-1))/kboltz
END SUBROUTINE getekin

! compute forces 
SUBROUTINE force
  USE kinds
  USE utils
  USE mdsys
  USE cell
  IMPLICIT NONE

  REAL(kind=dbl) :: rsq, rcutsq, rinv, pos1(3), delta(3)
  REAL(kind=dbl) :: boxby2, c12, c6, r6, ffac
  INTEGER :: i, j, k, n, m, ii, jj, kk, tid, fromidx, toidx
  INTEGER, EXTERNAL :: omp_get_thread_num
 
  epot=0.0_dbl

  !$OMP parallel default(SHARED) reduction(+:epot)              &
  !$OMP private(i,j,k,n,m,ii,jj,kk,tid,fromidx,toidx)           &
  !$OMP private(boxby2,c12,c6,r6,ffac,rsq,rcutsq,rinv,pos1,delta)
  tid = 0
  !$ tid = omp_get_thread_num() 
  tid = tid + 1
  frc(:,:,tid) = 0.0_dbl

  ! precompute some constants
  boxby2=0.5_dbl*box
  rcutsq=rcut*rcut
  c12 = 4.0_dbl*epsilon*sigma**12
  c6  = 4.0_dbl*epsilon*sigma**6

  ! first compute per cell self-interactions
  DO kk=0, ncell-1, nthreads
     i = kk + tid
     IF (i > ncell) EXIT

     DO j=1,npercell(i)-1
        ii = clist(i,j)
        pos1 = pos(ii,:)

        DO k=j+1,npercell(i)
           jj = clist(i,k)
           delta(1)=pbc(pos1(1)-pos(jj,1), boxby2, box)
           delta(2)=pbc(pos1(2)-pos(jj,2), boxby2, box)
           delta(3)=pbc(pos1(3)-pos(jj,3), boxby2, box)
           rsq = dot_product(delta,delta)
      
           ! compute force and energy if within cutoff */
           IF (rsq < rcutsq) THEN
              rinv = 1.0_dbl/rsq
              r6 = rinv*rinv*rinv
              ffac = (12.0_dbl*c12*r6 - 6.0_dbl*c6)*r6*rinv
              epot = epot + r6*(c12*r6 - c6)

              frc(ii,:,tid) = frc(ii,:,tid) + delta*ffac
              frc(jj,:,tid) = frc(jj,:,tid) - delta*ffac
           END IF
        END DO
     END DO
  END DO

  ! now compute per cell-cell interactions from pair list
  DO kk=0, npair-1, nthreads
     n = kk + tid
     IF (n > npair) EXIT
        
     i = plist(2*n-1)
     m = plist(2*n)
     DO j=1,npercell(i)
        ii = clist(i,j)
        pos1 = pos(ii,:)

        DO k=1,npercell(m)
           jj = clist(m,k)
           delta(1)=pbc(pos1(1)-pos(jj,1), boxby2, box)
           delta(2)=pbc(pos1(2)-pos(jj,2), boxby2, box)
           delta(3)=pbc(pos1(3)-pos(jj,3), boxby2, box)
           rsq = dot_product(delta,delta)
      
           ! compute force and energy if within cutoff */
           IF (rsq < rcutsq) THEN
              rinv = 1.0_dbl/rsq
              r6 = rinv*rinv*rinv
              ffac = (12.0_dbl*c12*r6 - 6.0_dbl*c6)*r6*rinv
              epot = epot + r6*(c12*r6 - c6)

              frc(ii,:,tid) = frc(ii,:,tid) + delta*ffac
              frc(jj,:,tid) = frc(jj,:,tid) - delta*ffac
           END IF
        END DO
     END DO
  END DO
  ! before reducing the forces, we have to make sure 
  ! that all threads are done adding to them.
  !$OMP barrier

  IF (nthreads > 1) THEN
     ! set equal chunks of index ranges
     i = 1 + (natoms/nthreads)
     fromidx = (tid-1)*i + 1
     toidx = fromidx + i - 1
     IF (toidx > natoms) toidx = natoms

     ! now reduce forces from threads with tid > 1 into
     ! the storage of the first thread. since we have
     ! threads already spawned, we do this in parallel.
     DO i=2,nthreads
        DO j=fromidx,toidx
           frc(j,:,1) = frc(j,:,1) + frc(j,:,i)
        END DO
     END DO
  END IF
  !$OMP END PARALLEL
END SUBROUTINE force


! velocity verlet
SUBROUTINE velverlet
  USE kinds
  USE mdsys
  USE physconst
  IMPLICIT NONE

  REAL(kind=dbl) :: vfac

  vfac = 0.5_dbl * dt / mvsq2e / mass

  ! first part: propagate velocities by half and positions by full step
  vel(:,:) = vel(:,:) + vfac*frc(:,:,1)
  pos(:,:) = pos(:,:) + dt*vel(:,:)

  ! compute forces and potential energy 
  CALL force

  ! second part: propagate velocities by another half step */
  vel(:,:) = vel(:,:) + vfac*frc(:,:,1) 
END SUBROUTINE velverlet


PROGRAM LJMD
  USE kinds
  USE io
  USE utils
  USE mdsys
  USE cell
  IMPLICIT NONE
  
  INTEGER :: nprint, i, j
  INTEGER, EXTERNAL :: omp_get_num_threads
  CHARACTER(len=sln) :: restfile, trajfile, ergfile

  nthreads = 1
  !$OMP parallel shared(nthreads)
  !$OMP master
  !$  nthreads = omp_get_num_threads()
  !$  WRITE(stdout,'(A,I2,A)') 'Running OpenMP version using ',nthreads,' thread(s).'
  !$OMP end master
  !$OMP end parallel

  READ(stdin,*) natoms
  READ(stdin,*) mass
  READ(stdin,*) epsilon
  READ(stdin,*) sigma
  READ(stdin,*) rcut
  READ(stdin,*) box
  CALL getline(stdin,restfile)
  CALL getline(stdin,trajfile)
  CALL getline(stdin,ergfile)
  READ(stdin,*) nsteps
  READ(stdin,*) dt
  READ(stdin,*) nprint

  ! allocate storage for simulation data.
  ALLOCATE(pos(natoms,3),vel(natoms,3),frc(natoms,3,nthreads))


  ! read restart 
  OPEN(UNIT=33, FILE=restfile, FORM='FORMATTED', STATUS='OLD')
  DO i=1,natoms
     READ(33,*) (pos(i,j),j=1,3)
  END DO
  DO i=1,natoms
     READ(33,*) (vel(i,j),j=1,3)
  END DO
  CLOSE(33)

  ! set up cell list
  CALL mkcell
  CALL updcell
  
  ! initialize forces and energies
  nfi=0
  frc(:,:,:) = 0.0_dbl
  CALL force
  CALL getekin
    
  CALL ioopen(ergfile, trajfile)

  WRITE(stdout, *) 'Starting simulation with ', natoms, ' atoms for', nsteps, ' steps'
  WRITE(stdout, *) '    NFI           TEMP                 EKIN                  EPOT&
       &                ETOT'
  CALL output

  ! main MD loop 
  DO nfi=1, nsteps
     ! write output, if requested
     IF (mod(nfi,nprint) == 0) THEN
        CALL output
     END IF

     ! propagate system and recompute energies
     CALL updcell
     CALL velverlet
     CALL getekin
  END DO

  ! clean up: close files, free memory
  WRITE(stdout,'(A)') 'Simulation Done.'
  CALL rmcell
  CALL ioclose

  DEALLOCATE(pos,vel,frc)
END PROGRAM LJMD
