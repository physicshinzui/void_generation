module Enum_Type
!

! calculate factor
!
  integer,parameter::TOTAL_POT     = 1
  integer,parameter::BOND_POT      = 2
  integer,parameter::ANGLE_POT     = 3
  integer,parameter::TORSION_POT   = 4
  integer,parameter::IMPROPER_POT  = 5

  integer,parameter::VDW14_POT     = 6
  integer,parameter::ELE14_POT     = 7
  integer,parameter::VDW15_POT     = 8
  integer,parameter::ELE15_POT     = 9
  integer,parameter::HYD15_POT     = 10

  character(9),dimension(HYD15_POT),parameter::Energy_Name =                   &
      (/'TOTAL ', 'BOND  ', 'ANGLE ', 'TORS. ', 'IMPRO.',                      &
        'VDW14 ', 'ELE14 ', 'VDW15 ', 'ELE15 ', 'HYD.  '/)

 end module Enum_Type
