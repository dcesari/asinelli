MODULE asinelli
IMPLICIT NONE

INTEGER,PARAMETER :: wp=KIND(1.0), vers=1

TYPE layer_t
  PRIVATE
  TYPE(layer_t),POINTER,PUBLIC :: next=>NULL()
  INTEGER,PUBLIC :: ninp=0, nout=0, ninconn=0
  REAL(kind=wp),ALLOCATABLE :: w(:,:), b(:)
  INTEGER,ALLOCATABLE :: inconnstart(:), inconnend(:)
  CONTAINS
  GENERIC :: init=>layer_t_init
  PROCEDURE,PRIVATE :: layer_t_init
  GENERIC :: compute=>layer_t_compute
  PROCEDURE,PRIVATE :: layer_t_compute
  GENERIC :: display=>layer_t_display
  PROCEDURE,PRIVATE :: layer_t_display
  GENERIC :: write_to_file=>layer_t_write_to_file
  PROCEDURE,PRIVATE :: layer_t_write_to_file
  PROCEDURE,NOPASS :: trans_func
  
END TYPE layer_t

INTEGER :: ninp, nout, ninconn
NAMELIST/layer_def/ninp, nout, ninconn

TYPE nnet_t
  PRIVATE
  TYPE(layer_t),POINTER :: firstlayer=>NULL()
  INTEGER,PUBLIC :: nlayer=0, nmaxio=0, ninp=0, nout=0
  REAL(kind=wp),ALLOCATABLE ::  buff1(:), buff2(:)
  CONTAINS
  GENERIC :: init=>nnet_t_init
  PROCEDURE,PRIVATE :: nnet_t_init
!  GENERIC :: compute=>nnet_t_compute
!  PROCEDURE,PRIVATE :: nnet_t_compute
!  GENERIC :: display=>nnet_t_display
!  PROCEDURE,PRIVATE :: nnet_t_display
!  GENERIC :: write_to_file=>nnet_t_write_to_file
!  PROCEDURE,PRIVATE :: nnet_t_write_to_file
END TYPE nnet_t

INTEGER :: nlayer
NAMELIST/nnet_def/nlayer

CONTAINS

SUBROUTINE nnet_t_init(this, next, unit)
CLASS(nnet_t),INTENT(out) :: this
INTEGER,INTENT(in),OPTIONAL :: unit

TYPE(layer_t),POINTER :: curr, next
INTEGER :: i, ierr

IF (PRESENT(unit)) THEN
  READ(unit, nml=nnet_def, iostat=ierr)
  ALLOCATE(this%firstlayer)

  this%nlayer = nlayer
  curr => this%firstlayer
  DO i = 1, nlayer
    IF (i == 1) THEN
      thix%nin = curr%nin
    ENDIF
    IF (i < nlayer) THEN
      CALL curr%init(next, unit=unit)
    ELSE
      CALL curr%init(unit=unit)
      this%nout = curr%nout
    ENDIF
    IF (.NOT.ASSOCIATED(next)) THEN
      CALL pflow_error('nnet_init: end of file reading layer namelist')
    ENDIF
    this%nmaxio = MAX(this%nmaxio, curr%ninp, curr%nout)
    curr => next
  ENDDO
ENDIF
ALLOCATE(this%buff1(this%nmaxio), this%buff2(this%nmaxio))

END SUBROUTINE nnet_t_init


SUBROUTINE layer_t_init(this, next, unit, deffile)
CLASS(layer_t),INTENT(out) :: this
TYPE(layer_t),POINTER,OPTIONAL :: next
INTEGER,INTENT(in),OPTIONAL :: unit
CHARACTER(len=*),INTENT(in),OPTIONAL :: deffile

INTEGER :: lvers, ierr

IF (PRESENT(next)) THEN
  ALLOCATE(this%next)
  next => this%next
ENDIF
IF (PRESENT(unit)) THEN
  ninp=-1; nout=-1; ninconn=-1

  READ(unit, nml=layer_def, iostat=ierr)
  IF (nout <=0) THEN
    CALL pflow_error('nout must be >0 for each layer')
  ENDIF
  IF (ninp <=0) THEN
    CALL pflow_error('ninp must be >0 for each layer')
  ENDIF
  IF (ninconn <= 0) THEN
    ninconn = ninp
  ENDIF
  this%nout = nout
  this%ninp = ninp
  this%ninconn = MIN(ninconn, ninp)
  ALLOCATE(this%w(this%nout, this%ninconn), this%b(this%nout))
  this%w(:,:) = 1.0_wp/this%ninconn
  this%b(:) = 0.0_wp
  CALL def_limits()

ELSE IF (PRESENT(deffile)) THEN
  OPEN(10, file=deffile, form='formatted')
  READ(10,*) lvers,this%nout,this%ninp,this%ninconn
  ALLOCATE(this%w(this%nout, this%ninconn), this%b(this%nout))
  READ(10,*) this%w,this%b
  CLOSE(10)
  CALL def_limits()

ENDIF

CONTAINS

SUBROUTINE def_limits()
INTEGER :: i

IF (this%ninconn < this%ninp) THEN
  ALLOCATE(this%inconnstart(this%nout), this%inconnend(this%nout))
  DO i = 1, this%nout
    this%inconnstart(i) = &
     NINT((1.0_wp/(1-this%nout))*((i-this%nout)-(i-1)*(this%ninp-this%ninconn+1)))
    this%inconnend(i) = this%inconnstart(i) + this%ninconn-1
  ENDDO
ENDIF

END SUBROUTINE def_limits
END SUBROUTINE layer_t_init


SUBROUTINE nnet_t_compute(this, in, out)
CLASS(nnet_t),INTENT(in) :: this
REAL(kind=wp),INTENT(in) :: in(:)
REAL(kind=wp),INTENT(out) :: out(:)

TYPE(layer_t), POINTER :: curr
REAL(kind=wp),POINTER :: pin(:)
REAL(kind=wp),POINTER :: pout

curr => this%firstlayer

CALL curr%compute(in, this%buff1)

DO WHILE(ASSOCIATED(curr))
  
END DO
END SUBROUTINE nnet_t_compute

SUBROUTINE layer_t_compute(this, in, out)
CLASS(layer_t),INTENT(in) :: this
REAL(kind=wp),INTENT(in) :: in(:)
REAL(kind=wp),INTENT(out) :: out(:)

INTEGER :: i

#ifdef DEBUG
IF (SIZE(in) /= this%ninp) THEN
  CALL pflow_error('compute: in must be of size ninp')
ENDIF
IF (SIZE(out) /= this%nout) THEN
  CALL pflow_error('compute: out must be of size nout')
ENDIF
#endif

DO i = 1, this%nout
  out(i) = this%trans_func( &
   DOT_PRODUCT(this%w(i,this%inconnstart(i):this%inconnend(i)), in(:)) + &
   this%b(i))
ENDDO


!  out(:) = this%trans_func(MATMUL(this%w, in) + this%b(:))

END SUBROUTINE layer_t_compute


SUBROUTINE layer_t_display(this)
CLASS(layer_t),INTENT(in) :: this

PRINT*,this%ninp, this%nout, this%ninconn

END SUBROUTINE layer_t_display


SUBROUTINE layer_t_write_to_file(this, deffile)
CLASS(layer_t),INTENT(in) :: this
CHARACTER(len=*),INTENT(in) :: deffile

OPEN(10, file=deffile, form='formatted')
WRITE(10,*) vers,this%nout,this%ninp,this%ninconn
IF (ALLOCATED(this%w)) WRITE(10,*) this%w,this%b
IF (ALLOCATED(this%inconnstart)) WRITE(10,*) this%inconnstart,this%inconnend
CLOSE(10)

END SUBROUTINE layer_t_write_to_file


FUNCTION trans_func(x) RESULT(y)
REAL,INTENT(in) :: x
REAL :: y

y = 1.0_wp/(1.0_wp+EXP(-y))
END FUNCTION trans_func


SUBROUTINE pflow_error(message)
CHARACTER(len=*),INTENT(in) :: message
WRITE(*,'(A)') message
STOP 1
END SUBROUTINE pflow_error

END MODULE asinelli

