PROGRAM asino
USE asinelli
IMPLICIT NONE

TYPE(nnet_t) :: nnet

OPEN(10, file='nnet.naml')
CALL nnet%init(unit=10)
CLOSE(10)

OPEN(10, file='nnet.out')
CALL nnet%write_to_file(10)

END PROGRAM asino
