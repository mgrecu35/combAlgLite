
module keras_def
  use mod_kinds, only: ik, rk
  use mod_network, only: network_type

  implicit none

  type(network_type) :: net

  real(rk), allocatable :: result1(:), input(:)

end module keras_def
