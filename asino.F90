PROGRAM asino
USE asinelli
IMPLICIT NONE

TYPE(layer_t) :: layer

OPEN(10, file='layer.naml')
CALL layer%init(unit=10)

CALL layer%write_to_file('layer.out')

END PROGRAM asino
