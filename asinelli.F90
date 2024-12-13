MODULE asinelli
IMPLICIT NONE

INTEGER,PARAMETER :: wp=KIND(1.0), vers=1

TYPE layer_t
  PRIVATE
  TYPE(layer_t),POINTER :: next=>NULL()
  INTEGER :: ninp=0, nout=0, ninconn=0
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
  TYPE(layer_t),POINTER :: firstlayer=>NULL()
END TYPE nnet_t


CONTAINS

SUBROUTINE layer_t_init(this, unit, namlfile, deffile)
CLASS(layer_t),INTENT(out) :: this
INTEGER,INTENT(in),OPTIONAL :: unit
CHARACTER(len=*),INTENT(in),OPTIONAL :: namlfile
CHARACTER(len=*),INTENT(in),OPTIONAL :: deffile

INTEGER :: lunit, lvers

IF (PRESENT(namlfile) .OR. PRESENT(unit)) THEN
  ninp=-1; nout=-1; ninconn=-1
  IF (PRESENT(namlfile)) THEN
    OPEN(10, file=namlfile)
    lunit = 10
  ELSE
    lunit = unit
  ENDIF

  READ(lunit, nml=layer_def)
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

  IF (.NOT.PRESENT(unit)) CLOSE(lunit)

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

