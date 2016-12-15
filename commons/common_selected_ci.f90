module common_selected_ci
  use types, only : rk,ik,ik_vec,i8b

  integer :: det_sel_iters ! number of times find connected dets and truncation are performed
  integer,allocatable :: norb_det_sel(:),n_sym_uniq_det_det_sel(:)
  integer :: ndet_det_sel
  real(rk), allocatable :: cdet_det_sel(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable :: dets_up_det_sel(:), dets_dn_det_sel(:)
#else
  integer(ik), allocatable :: dets_up_det_sel(:), dets_dn_det_sel(:)
#endif
  integer :: lanczos_iters,lanczos_initiators,lanczos_truncate
  integer(i8b),parameter :: max_nonzero_elements=int(5e9,i8b) ! Maximum number of nonzero matrix elements that can be stored in a sparse matrix for Lanczos
  logical :: too_big_to_store
  real(rk) :: log_num_nonzero_elements
  ! These are for Norm's IACI algorithm which is no longer used:
  integer :: cdets,tdets
  ! These are for our HCI algorithm:
  real(rk) :: eps_var ! Cutoff for terms in the sum for 1PT coefficients (typically 1e-3 - 2e-5)
  real(rk) :: eps_pt  ! Cutoff for terms in the sum for 2PT energy correction (typically 1e-6 - 1e-8)
  real(rk) :: target_error ! Maximum standard deviation in 2PT energy correction when computed stochastically (typically 5e-4)
  logical :: dump_wf_var
  type sparse_mat
    integer :: ndet = 0
    integer(i8b),allocatable :: indices(:),nonzero_elements(:)
    real(rk),allocatable :: values(:)
  end type sparse_mat
  type(sparse_mat), save :: sparse_ham
  integer :: n_states ! number of ground and excited states to compute with HCI
end module common_selected_ci
