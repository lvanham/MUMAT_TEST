PROGRAM MUMATERIAL_TEST
   USE MUMATERIAL_MOD
   IMPLICIT NONE
   INCLUDE "mpif.h"


   CHARACTER(LEN=256) :: filename
   CHARACTER(LEN=256) :: nearest
   CHARACTER(LEN=256) :: distance
   CHARACTER(LEN=256) :: grid

   DOUBLE PRECISION, DIMENSION(:), allocatable :: x, y, z, Hx, Hy, Hz, offset
   DOUBLE PRECISION, DIMENSION(:), allocatable :: Bx, By, Bz
   DOUBLE PRECISION, DIMENSION(:), allocatable :: Bx_local, By_local, Bz_local
   INTEGER :: i_int, nn
   DOUBLE PRECISION :: dist
   INTEGER :: start, finish, rate

   integer, parameter :: arg_len = 256

   INTEGER :: i, numargs, N_points
   CHARACTER*(arg_len) :: arg1
   CHARACTER*(arg_len), allocatable, dimension(:) :: args

   INTEGER :: istat, mysize, rank, comm
   INTEGER :: mystart, myend

   rank = 0

    !-----------------------------------------------------------------------
    !     Handle Input Arguments
    !-----------------------------------------------------------------------
    numargs = 0
    i = 0
    arg1 = ''
    nn = -1
    dist = -1

    ! First Handle the input arguments
    CALL GETCARG(1, arg1, numargs)
    ALLOCATE(args(numargs))
    ! Cycle through Arguments
    i = 1
    DO WHILE (i <= numargs)
       call GETCARG(i, args(i), numargs)
       select case (args(i))
          case ("-mumat")
                i = i + 1
                CALL GETCARG(i, filename, numargs)
          case ("-nearest")
                i = i + 1
                CALL GETCARG(i, nearest,  numargs)
		          read (nearest, '(I10)') nn
          case ("-distance")
                i = i + 1
                CALL GETCARG(i, distance, numargs)
                read (distance, '(F10.0)') dist
          case ("-grid")
                i = i + 1
                CALL GETCARG(i, grid, numargs)

       END SELECT
       i = i + 1
    END DO
    DEALLOCATE(args)

   CALL MPI_INIT(istat)
   comm = MPI_COMM_WORLD
   CALL MPI_COMM_SIZE(comm, mysize, istat)
   CALL MPI_COMM_RANK(comm, rank, istat)
   
   CALL MUMATERIAL_SETVERB(.FALSE.)
   IF (rank .eq. 0) THEN
      CALL SYSTEM_CLOCK(count_rate=rate)
      CALL SYSTEM_CLOCK(start)
      CALL MUMATERIAL_SETVERB(.TRUE.)
   END IF
   IF (rank .eq. 0) WRITE(*,*) dist
      allocate(offset(3))
      offset = [0.0, 0.0, 0.0]


      CALL MUMATERIAL_LOAD(TRIM(filename),istat, comm)
      ! if (istat/=0) EXIT(2) ! probably need to stop the program in this case?

      CALL MUMATERIAL_SETD(1.0d-5, 1000, 0.7d0, 0.75d0, nn, dist, comm) ! only set if values need to be changed
      
      IF (rank .eq. 0) THEN
         CALL MUMATERIAL_INFO(6)
      END IF
      CALL MUMATERIAL_INIT_NEW(BEXTERNAL, comm, offset)

      IF (rank .eq. 0) THEN
         CALL SYSTEM_CLOCK(finish)
         WRITE(*,*) "Time to finish loading: ", real(finish-start)/real(rate)
      END IF
      
      ! CALL gen_grid(x, y, z)
      IF (rank .eq. 0) WRITE(6,*) 'Reading grid'
      CALL read_grid(grid, x, y, z, istat)
      n_points = size(x)
      
      IF (rank .eq. 0) WRITE(6,*) 'Calculating B-field'

      allocate(Bx(n_points),By(n_points),Bz(n_points))
      allocate(Bx_local(n_points),By_local(n_points),Bz_local(n_points))

      CALL MPI_CALC_MYRANGE(comm, 1, n_points, mystart, myend)
      DO i = mystart, myend
         CALL mumaterial_getbmag_scalar(x(i), y(i), z(i), Bx_local(i), By_local(i), Bz_local(i))
      END DO

      CALL MPI_ALLREDUCE(Bx_local, Bx, n_points, MPI_DOUBLE_PRECISION, MPI_SUM, comm, istat)
      CALL MPI_ALLREDUCE(By_local, By, n_points, MPI_DOUBLE_PRECISION, MPI_SUM, comm, istat)
      CALL MPI_ALLREDUCE(Bz_local, Bz, n_points, MPI_DOUBLE_PRECISION, MPI_SUM, comm, istat)

      deallocate(Bx_local,By_local,Bz_local)

      IF (rank .eq. 0) THEN
         CALL SYSTEM_CLOCK(finish)
         WRITE(*,*) "Time to finish B-field calculations: ", real(finish-start)/real(rate)
      END IF

      IF (rank .eq. 0) THEN
         WRITE(6,*) "Outputting B-field"
         OPEN(14, file='./B.dat')
         DO i = 1, n_points
               WRITE(14, "(E15.7,A,E15.7,A,E15.7)") Bx(i), ',', By(i), ',', Bz(i)
         END DO
         CLOSE(14)
      END IF

      ! CALL MUMATERIAL_GETB(5.d0, 5.d0, 501.d0, Bx, By, Bz, BEXTERNAL)
      ! WRITE(*,*) "H:", Bx / (16 * atan(1.d0) * 1.d-7), By / (16 * atan(1.d0) * 1.d-7), Bz / (16 * atan(1.d0) * 1.d-7)

      !CALL MUMATERIAL_OUTPUT('./', x, y, z, BEXTERNAL, comm)




      CALL MUMATERIAL_FREE()

      IF (rank .eq. 0) THEN
         CALL SYSTEM_CLOCK(finish)
         WRITE(*,*) "Time to finish: ", real(finish-start)/real(rate)
      END IF

   CALL MPI_FINALIZE(istat)

   CONTAINS

   SUBROUTINE BEXTERNAL(x,y,z,bx,by,bz)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: x,y,z
      DOUBLE PRECISION, INTENT(OUT) :: bx,by,bz
      DOUBLE PRECISION :: theta, phi, mag

      ! python easy axis from spherical coordinate transformations
      theta = 0.0
      phi = 0.0
      mag = 1.0!16 * atan(1.d0) * 1.0E-7
      bx = mag*sin(theta)*cos(phi)
      by = mag*sin(theta)*sin(phi)
      bz = mag*cos(theta)

      !bx = 1.0; by = 0.0; bz = 0.0

      RETURN
   END SUBROUTINE BEXTERNAL

   subroutine read_grid(gridfile,x,y,z,istat)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: gridfile
      DOUBLE PRECISION, dimension(:), allocatable, intent(out) :: x, y, z
      INTEGER, INTENT(inout) :: istat
      INTEGER :: iunit, lines, il

      iunit = 380; istat = 0; lines = 0
      CALL safe_open(iunit,istat,TRIM(gridfile),'old','formatted')
      IF (istat /= 0) RETURN

      READ(iunit,*) lines

      allocate(x(lines))
      allocate(y(lines))
      allocate(z(lines))

      DO il = 1, lines
         READ(iunit,*) x(il),y(il),z(il)
      END DO

   end subroutine read_grid

END PROGRAM MUMATERIAL_TEST
