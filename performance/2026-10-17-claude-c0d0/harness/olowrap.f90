! C-bound wrappers of OneLOop 3.6, identical to Ninja's avholo_wrapper.F90 (_cm entries)
subroutine olw_d0_cm(rslt,p1,p2,p3,p4,p12,p23,m1,m2,m3,m4,rmu) bind (c)
  use, intrinsic :: iso_c_binding
  use avh_olo
  implicit none
  complex(c_double), dimension(0:2), intent(out) :: rslt
  complex(c_double), intent(in)  :: p1,p2,p3,p4,p12,p23,m1,m2,m3,m4
  real(c_double), intent(in) :: rmu
  call olo_d0(rslt,p1,p2,p3,p4,p12,p23,m1,m2,m3,m4,rmu)
end subroutine
subroutine olw_c0_cm(rslt,p1,p2,p3,m1,m2,m3,rmu) bind (c)
  use, intrinsic :: iso_c_binding
  use avh_olo
  implicit none
  complex(c_double), dimension(0:2), intent(out) :: rslt
  complex(c_double), intent(in)  :: p1,p2,p3,m1,m2,m3
  real(c_double), intent(in) :: rmu
  call olo_c0(rslt,p1,p2,p3,m1,m2,m3,rmu)
end subroutine
subroutine olw_b0_cm(rslt,pp,m1,m2,rmu) bind (c)
  use, intrinsic :: iso_c_binding
  use avh_olo
  implicit none
  complex(c_double), dimension(0:2), intent(out) :: rslt
  complex(c_double), intent(in)  :: pp,m1,m2
  real(c_double), intent(in) :: rmu
  call olo_b0(rslt,pp,m1,m2,rmu)
end subroutine
subroutine olw_b11_cm(b11,b00,b1,b0,pp,m1,m2,rmu) bind (c)
  use, intrinsic :: iso_c_binding
  use avh_olo
  implicit none
  complex(c_double), dimension(0:2), intent(out) :: b11,b00,b1,b0
  complex(c_double), intent(in)  :: pp,m1,m2
  real(c_double), intent(in) :: rmu
  call olo_b11(b11,b00,b1,b0,pp,m1,m2,rmu)
end subroutine
subroutine olw_a0_cm(rslt,mm,rmu) bind (c)
  use, intrinsic :: iso_c_binding
  use avh_olo
  implicit none
  complex(c_double), dimension(0:2), intent(out) :: rslt
  complex(c_double), intent(in)  :: mm
  real(c_double), intent(in) :: rmu
  call olo_a0(rslt,mm,rmu)
end subroutine
subroutine olw_setup(thrs) bind (c)
  use, intrinsic :: iso_c_binding
  use avh_olo
  implicit none
  real(c_double), intent(in) :: thrs
  call olo_unit(-1,'warning')
  call olo_onshell(thrs)
end subroutine
