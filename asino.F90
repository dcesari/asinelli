PROGRAM asino
USE asinelli
IMPLICIT NONE

TYPE(layer_t) :: layer

CALL layer%init(namlfile='layer.naml')

CALL layer%write_to_file('layer.out')

END PROGRAM asino
