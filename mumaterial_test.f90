PROGRAM MUMATERIAL_TEST
   USE MUMATERIAL_MOD
   IMPLICIT NONE
   INCLUDE 'mpif.h'


   CHARACTER(LEN=256) :: filename
   CHARACTER(LEN=256) :: padding
   CHARACTER(LEN=256) :: lambdastart
   CHARACTER(LEN=256) :: lambdafactor
   CHARACTER(LEN=256) :: lambdacount
   CHARACTER(LEN=256) :: maxerror
   CHARACTER(LEN=256) :: maxiter
   CHARACTER(LEN=256) :: convcheck
   CHARACTER(LEN=256) :: Mfile
   CHARACTER(LEN=256) :: tree_depth
   CHARACTER(LEN=256) :: tree_leaf
   CHARACTER(LEN=256) :: tree_theta


   DOUBLE PRECISION :: pad, lambdaS, lambdaF, maxerr, cc, theta
   INTEGER :: lambdaC, maxi, depth, leaf

   INTEGER :: istat, comm_world, shar_comm, comm_master
   INTEGER :: shar_rank, master_rank
   LOGICAL :: lismaster, ldebug, lhasMfile, lnoiter, loutmag


   DOUBLE PRECISION, DIMENSION(:), allocatable :: x, y, z, Hx, Hy, Hz, offset
   INTEGER :: i_int
   DOUBLE PRECISION :: Bx, By, Bz
   INTEGER :: start, finish, rate

   integer, parameter :: arg_len = 256

   INTEGER :: i, numargs
   CHARACTER*(arg_len) :: arg1
   CHARACTER*(arg_len), allocatable, dimension(:) :: args

   shar_rank = 0
   master_rank = 0
   lismaster = .FALSE.
   !-----------------------------------------------------------------------
   !     Handle Input Arguments
   !-----------------------------------------------------------------------
   numargs = 0
   i = 0
   arg1 = ''
   pad = 1.0
   lambdaS = 0.9
   lambdaF = 0.99
   lambdaC = 9
   maxerr = 9.9d-3
   maxi = 999
   cc = 99.9
   lhasMfile = .FALSE.
   lnoiter = .FALSE.
   loutmag = .FALSE.
   theta = 0.20d0
   depth = 30
   leaf = 4

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
          case ("-padfactor")
            i = i + 1
            CALL GETCARG(i, padding, numargs)
            read (padding, '(F15.0)') pad
          case ("-lambdastart")
            i = i + 1
            CALL GETCARG(i, lambdastart, numargs)
            read (lambdastart, '(F15.0)') lambdaS
          case ("-lambdafactor")
            i = i + 1
            CALL GETCARG(i, lambdafactor, numargs)
            read (lambdafactor, '(F15.0)') lambdaF
          case ("-lambdacount")
            i = i + 1
            CALL GETCARG(i, lambdacount, numargs)
            read (lambdacount, '(I7)') lambdaC
          case ("-maxerror")
            i = i + 1
            CALL GETCARG(i, maxerror, numargs)
            read (maxerror, '(F15.0)') maxerr
          case ("-maxiter")
            i = i + 1
            CALL GETCARG(i, maxiter, numargs)
            read (maxiter, '(I7)') maxi
         case ("-convcheck")
            i = i + 1
            CALL GETCARG(i, convcheck, numargs)
            read (convcheck, '(F15.0)') cc
         case ("-magfile")
            i = i + 1
            CALL GETCARG(i, Mfile, numargs)
            lhasMfile = .TRUE.
         case ("-theta")
            i = i + 1
            CALL GETCARG(i, tree_theta, numargs)
            read (tree_theta, '(F15.0)') theta
         case ("-depth")
            i = i + 1
            CALL GETCARG(i, tree_depth, numargs)
            read (tree_depth, '(I7)') depth
         case ("-leaf")
            i = i + 1
            CALL GETCARG(i, tree_leaf, numargs)
            read (tree_leaf, '(I7)') leaf
         case ("-noiter")
            i = i + 1
            lnoiter = .TRUE.
         case ("-outmag")
            i = i + 1
            loutmag = .TRUE.
      END SELECT
      i = i + 1
   END DO
   DEALLOCATE(args)

   CALL MPI_INIT(istat)
   comm_world = MPI_COMM_WORLD
   CALL MUMATERIAL_SET_COMMS(comm_world, shar_comm, comm_master)
   CALL MPI_COMM_RANK( shar_comm, shar_rank, istat)
   ldebug = (shar_rank.EQ.0)
   IF (shar_rank.EQ.0) THEN
      CALL MPI_COMM_RANK( comm_master, master_rank, istat)
      lismaster = (master_rank.EQ.0)
   END IF

   CALL MUMATERIAL_SETVERB(lismaster)
 
   allocate(offset(3))
   offset = [0.0, 0.0, 0.0]

   CALL MUMATERIAL_LOAD(TRIM(filename),istat, shar_comm, comm_master,comm_world)
   CALL MUMATERIAL_SET_VARS(max_error=maxerr, max_iter = maxi, &
         max_depth = depth, max_leafsize = leaf, iter_theta=theta, eval_theta=theta, &
         min_conv_perc=cc, lambda_start = lambdaS, lambda_factor = lambdaF) 
   IF (lhasMfile) CALL MUMATERIAL_MAGFILE_READ(TRIM(Mfile))

   IF (lismaster) CALL MUMATERIAL_INFO(6, lnoiter)
   CALL MPI_BARRIER(comm_world, istat)

   IF (NOT(lnoiter)) THEN
      IF (lismaster) THEN
         CALL SYSTEM_CLOCK(count_rate=rate)
         CALL SYSTEM_CLOCK(start)
      END IF

      CALL MUMATERIAL_RUN(BEXTERNAL, offset, linitM = .NOT.lhasMfile)

      IF (lismaster) THEN
         CALL SYSTEM_CLOCK(finish)
         WRITE(*,*) "Time to finish loading: ", real(finish-start)/real(rate)
         
         OPEN(14, file='./time.dat')
         WRITE(14,"(E15.7)") real(finish-start)/real(rate)
         CLOSE(14)
      END IF
   END IF
   
   CALL gen_grid(x, y, z)
   
   CALL MUMATERIAL_OUTPUT('./', x, y, z, .TRUE.)
   
   IF (loutmag) CALL MUMATERIAL_MAGFILE_WRITE('test')
   
   CALL MUMATERIAL_FREE()

   IF (lismaster.AND.NOT(lnoiter)) THEN
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

   subroutine gen_grid(x, y, z)
      implicit none
      DOUBLE PRECISION, dimension(:), allocatable, intent(out) :: x, y, z
      integer, dimension(3) :: num_points
      integer :: n_temp, i, j, k, n_points
      DOUBLE PRECISION :: r, theta, phi, pi, ohpointfive, onepointfive, two
      DOUBLE PRECISION, dimension(3) :: min, max

      pi = 4.0 * atan(1.0)
      ohpointfive = 0.5
      onepointfive = 1.5
      two = 2.0

      min = [0.0, 0.0, 0.0]
      max = [two, 2*pi, 2.0*pi]
      num_points = [2001, 361, 1]
      
      n_temp = 1
      n_points = num_points(1)*num_points(2)*num_points(3)
      allocate(x(n_points))
      allocate(y(n_points))
      allocate(z(n_points))

      do i = 1, num_points(1)
         do j = 1, num_points(2)
               do k = 1, num_points(3)
                  if (num_points(1) .gt. 1) then
                     r = min(1) + 1.0*(i-1)*(max(1)-min(1))/(num_points(1)-1)
                     x(n_temp) = min(1) + 1.0*(i-1)*(max(1)-min(1))/(num_points(1)-1)
                  else
                     r = min(1)
                     x(n_temp) = min(1)
                  end if
                  if (num_points(2) .gt. 1) then
                     theta = min(2) + 1.0*(j-1)*(max(2)-min(2))/(num_points(2)-1)
                     y(n_temp) = min(2) + 1.0*(j-1)*(max(2)-min(2))/(num_points(2)-1)
                  else
                     theta = min(2)
                     y(n_temp) = min(2)
                  end if
                  if (num_points(3) .gt. 1) then
                     phi = min(3) + 1.0*(k-1)*(max(3)-min(3))/(num_points(3))
                     z(n_temp) = min(3) + 1.0*(k-1)*(max(3)-min(3))/(num_points(3)-1)
                  else
                     phi = min(3)
                     z(n_temp) = min(3)
                  end if
                  x(n_temp) = r*sin(theta)*cos(phi)
                  y(n_temp) = r*sin(theta)*sin(phi)
                  z(n_temp) = r*cos(theta)
                  n_temp = n_temp + 1
               enddo
         enddo
      enddo
   end subroutine gen_grid

END PROGRAM MUMATERIAL_TEST
