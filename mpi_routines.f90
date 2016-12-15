module mpi_routines
  use types
  use common_walk, only: walk_dets_up, walk_dets_dn, walk_wt, matrix_elements, imp_distance
  !You have to use ISO_C_BINDING for the CPTR stuff
  use, intrinsic :: ISO_C_BINDING


implicit none
#ifdef MPI
  INCLUDE 'mpif.h'
#endif
public :: master_core, master_core_node
public :: cluster_init, cluster_finalize
public :: mpi_bsend, mpi_bsend_diffroot
public :: mpi_allred, mpi_gath, mpi_agathv, mpi_gathv, mpi_red_max, mpi_red_min, mpi_red, mpi_stop
public :: mpi_barr, t_sndlist, mpi_scattv, mpi_redscatt_real_dparray
public :: mpi_sendwalks, mpi_sendnewwalks, init_snd_buffers
public :: snd_table, snd_displs, mpi_push_nwalk, snd_cnt
public :: get_owner, init_hash_owners, get_det_owner, mpi_snd
public :: shmem_allocate

private

type t_sndlist
    integer :: to,wlkid
end type t_sndlist

type t_walk
    real(rk)     :: wt
    integer(i8b) :: det_u(2*num_words),det_d(2*num_words) ! 2 i8b's instead of 1 i16b because MPI can only handle size i8b
    integer(i1b) :: imp_ini(2) !contains both imp_distance and initiator
    real(rk) :: e_num,e_den,diag_elems
end type t_walk

!whoaminode is the rank within a node. master_core_node defines who is the master within a node
!sharedComm is the communicator within a node
integer,public :: whoami,mpierr,ncores,whoaminode,sharedComm
!declare shmemWindow here so it can be freed later
!A "window" is MPI's name for the memory that is open to other cores to see
integer :: shmemWindow(7) ! integrals + q,J,prob for each of same and opposite spin double excitations in heatbath
integer :: win_num
logical :: master_core,master_core_node
integer(i4b) :: nwbuff,loc_mwalk
integer,allocatable :: snd_table(:),snd_displs(:),snd_dsp(:),snd_cnt(:),rcv_cnt(:),recv_displs(:),list_cowner(:)
type(t_sndlist),allocatable :: snd_list(:)
integer istat
integer,parameter :: mpi_send_limit = 15000

!***buffers for send process [WARNING: put this in common_walk]
  real(rk), allocatable ::  snd_matrix_elements(:)

  type(t_walk), allocatable :: snd_buff(:),recv_buff(:)

  integer :: MPI_WALK_TYPE !** needed for communications

  integer :: snd_count,recv_count

interface shmem_allocate
  module procedure shmem_allocate_rk, shmem_allocate_real, shmem_allocate_int
end interface shmem_allocate

interface mpi_bsend
  module procedure mpi_bsend_int, mpi_bsend_real_sp,mpi_bsend_real_dp, mpi_bsend_complex,mpi_bsend_string,mpi_bsend_iarray, &!
                   mpi_bsend_logical,mpi_bsend_int64,mpi_bsend_irand_seed,mpi_bsend_real_dparray,mpi_bsend_i64arr,mpi_bsend_i128arr,mpi_bsend_ikvec,mpi_bsend_ikvec_arr
end interface

interface mpi_snd
  module procedure mpi_snd_int,mpi_snd_real_sp,mpi_snd_real_dp,mpi_snd_real_dparray,mpi_snd_iarray,mpi_snd_i64array,&
                   mpi_snd_int64,mpi_snd_i128arr,mpi_snd_ikvec,mpi_snd_ikvec_array
end interface

interface mpi_bsend_diffroot
  module procedure mpi_bsend_diffroot_int, mpi_bsend_diffroot_real_sp,mpi_bsend_diffroot_real_dp,mpi_bsend_diffroot_logical,&
                    mpi_bsend_diffroot_int64,mpi_bsend_diffroot_ik_vec,mpi_bsend_diffroot_int128,mpi_bsend_diffroot_int8,&
                    mpi_bsend_diffroot_sndlist,mpi_bsend_diffroot_iarray
end interface

interface mpi_gath
  module procedure mpi_gath_int, mpi_gath_real_sp, mpi_gath_real_dp, mpi_gath_int8 !, mpi_gath_ik_vec, mpi_gath_i128
end interface

interface mpi_gathv
  module procedure mpi_gathv_int, mpi_gathv_real_sp, mpi_gathv_real_dp, mpi_gathv_int8, mpi_gathv_ik_vec,mpi_gathv_i128
end interface

interface mpi_scattv
  module procedure mpi_scattv_int, mpi_scattv_real_sp, mpi_scattv_real_dp, mpi_scattv_int8, mpi_scattv_ik_vec,mpi_scattv_i128,mpi_scattv_iwalk
end interface

interface mpi_agathv
  module procedure mpi_agathv_int, mpi_agathv_real_sp, mpi_agathv_real_dp, mpi_agathv_int8, mpi_agathv_ik_vec, &
                   mpi_agathv_i128
end interface

interface mpi_allred
  module procedure mpi_allred_int, mpi_allred_real_sp,mpi_allred_real_dp,mpi_allred_real_dparray,mpi_allred_iarray!,mpi_bsend_complex,mpi_bsend_string, &
                 !  mpi_bsend_logical,mpi_bsend_int64,mpi_bsend_irand_seed!,mpi_bsend_3string
end interface

interface mpi_red_max
  module procedure mpi_red_max_real_sp, mpi_red_max_real_dp, mpi_red_max_real_dparray
end interface

interface mpi_red_min
  module procedure mpi_red_min_real_sp, mpi_red_min_real_dp, mpi_red_min_real_dparray
end interface

interface mpi_red
  module procedure mpi_red_real_sp,mpi_red_real_dp,mpi_red_real_dparray
end interface

interface conv_128_to_64
  module procedure conv_128_to_64,conv_064_to_64
end interface

interface conv_64_to_128
  module procedure conv_64_to_128,conv_64_to_064
end interface


!---- HASH TABLE part
integer :: RandomHash1(0:255) !,RandomHash2(0:255)
!---- HASH TABLE part
!--- for 128-bit integers send
integer*8,allocatable  :: i128_high(:),i128_low(:),i128_buff(:)
!--- for 128-bit integers send

contains

subroutine init_snd_buffers(mwalk)
    integer(i8b), intent(in) :: mwalk
    integer :: icore
#ifdef MPI
    loc_mwalk=int(mwalk,i4b)
    allocate(snd_list(mwalk),stat=istat)
    if(istat.ne.0) stop 'failed to allocate snd_list'
    allocate(list_cowner(mwalk),stat=istat)
    if(istat.ne.0) stop 'failed to allocate list_cowner'
    allocate(snd_matrix_elements(mwalk),stat=istat)
    if(istat.ne.0) stop 'failed to allocate snd_matrix_elements'

    allocate(snd_buff(mwalk),stat=istat)
    if(istat.ne.0) stop 'failed to allocate snd_buff'
    allocate(recv_buff(mwalk),stat=istat)
    if(istat.ne.0) stop 'failed to allocate recv_buff'

    nwbuff = int(mwalk/ncores,i4b)

    do icore=0,(ncores-1)
      snd_dsp(icore+1)=icore*nwbuff
    enddo

#endif

end subroutine init_snd_buffers
!------------------------------------------------------

!---- HASH TABLE part [based on George Booth's code]
subroutine init_hash_owners()
    implicit none
    integer :: i,j,map
    logical :: found
    real(rk) rannyu
#ifdef MPI
    if(ncores>1) then
      if(master_core) then
        do i=0,255
            RandomHash1(i) = i
        enddo

        !Find mapping integer function from range 0 -> 255 to unique integers in range 0 -> 255*20,000
        do i=0,255
            found=.false.
            do while(.not.found)
                map=int(255*rannyu()*20000)

                do j=0,i-1
                    if(RandomHash1(j).eq.map) exit
                enddo
                if(j.eq.i) found=.true.
            enddo
            RandomHash1(i)=map
        enddo
      endif
      call mpi_bsend(RandomHash1)
    endif
#endif
end subroutine
!------------------------------------------------------

!Return a hash between the range 0 -> 'range-1' from the bit-representation 'det'
subroutine hash(det,range,hashx,RandomHash)
    use types, only : ik
    implicit none
    integer(kind=ik), intent(in) :: det(2*num_words)
    integer, intent(in) :: range,RandomHash(0:255)
    integer, intent(out) :: hashx
    integer(kind=ik) :: acc
    integer(kind=ik) :: test_nk
    integer :: test_int,i,j
    integer :: val
    integer, parameter :: n_bits_nk = bit_size(test_nk)
    integer, parameter :: n_bits = bit_size(test_int)

    acc = 0
    do i=1,2*num_words    !run over integers defining the bit-string
        do j = 0, n_bits_nk -1, 8
            !val = int(iand(ishft(det(i),-j), int(255,ik)),kind(test_int))
            val = int(iand(ishft(det(i)*368296850_ik,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
            !val is now a number between 0 -> 255
            !Ensure that RandomHash has a mapping for the 0th element too
            !1099511628211_ik = ibset(0_ik,27)
            acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
        enddo
    enddo
    hashx = int(abs(mod(acc,int(range,ik))),kind(test_int))

end subroutine hash
!------------------------------------------------------

function get_det_owner(det_up,det_dn) result(coreid)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_up,det_dn
#else
    integer(kind=ik), intent(in) :: det_up,det_dn
#endif
    integer :: coreid
    integer(kind=ik) :: det(2*num_words)
#ifdef MPI
    if(ncores==1) then
      coreid=0
    else
#ifdef NUM_ORBITALS_GT_127
      det(1:num_words)=det_up%v(1:num_words)
      det(num_words+1:2*num_words)=det_dn%v(1:num_words)
#else
      det(1)=det_up
      det(2)=det_dn
#endif
      call hash(det,ncores,coreid,RandomHash1)
    endif
#else
    coreid = 0
#endif
end function get_det_owner
!------------------------------------------------------

function get_owner(windex,wtot) result(coreid)
integer, intent(in) :: windex,wtot
integer :: coreid
integer :: cscale,ic

    cscale=wtot/ncores
    coreid=-1
    do ic=0,(ncores-2)
        if((windex .ge. (ic*cscale+1)).and.(windex .le. ((ic+1)*cscale))) then
            coreid=ic
            return
        endif
    enddo
    if(coreid < 0) coreid=ncores-1 !*just for now that the last one has different upper-bound

!coreid=MOD(windex,2)
end function get_owner
!------------------------------------------------------
!---- End of HASH TABLE part --------------------------

subroutine init_mpitype_walker

  integer :: block_len(0:6)
  integer :: oldtypes(0:6),offsets(0:6)
  integer :: i8extent,r8extent,i1extent

#ifdef MPI
  call MPI_TYPE_EXTENT(MPI_INTEGER8, i8extent, mpierr)
  call MPI_TYPE_EXTENT(MPI_REAL8, r8extent, mpierr)
  call MPI_TYPE_EXTENT(MPI_INTEGER1, i1extent, mpierr)

  oldtypes(0) = MPI_REAL8    ! wt
  oldtypes(1) = MPI_INTEGER8 ! det_up
  oldtypes(2) = MPI_INTEGER8 ! det_dn
  oldtypes(3) = MPI_INTEGER1 ! combination of imp_distance and initiator (?)
  oldtypes(4) = MPI_REAL8    ! e_num
  oldtypes(5) = MPI_REAL8    ! e_den
  oldtypes(6) = MPI_REAL8    ! diagonal Hamiltonian element

  block_len(0) = 1
#ifdef NUM_ORBITALS_GT_127
  block_len(1) = 2*num_words
  block_len(2) = 2*num_words
#else
  block_len(1) = 2
  block_len(2) = 2
#endif
  block_len(3) = 2
  block_len(4) = 1
  block_len(5) = 1
  block_len(6) = 1

  offsets(0) = 0
  offsets(1) = r8extent*block_len(0)
  offsets(2) = offsets(1) + i8extent*block_len(1)
  offsets(3) = offsets(2) + i8extent*block_len(2)
  offsets(4) = offsets(3) + i1extent*block_len(3)
  offsets(5) = offsets(4) + r8extent*block_len(4)
  offsets(6) = offsets(5) + r8extent*block_len(5)

  call MPI_TYPE_STRUCT(7,block_len,offsets,oldtypes,MPI_WALK_TYPE,mpierr)
  call MPI_TYPE_COMMIT(MPI_WALK_TYPE,mpierr)
#endif

end subroutine init_mpitype_walker
!------------------------------------------------------

!*** these routines shall be moved elsewhere
subroutine conv_128_to_64(i128,i64_1,i64_2)
integer*16, intent(in) :: i128
integer*8, intent(out) :: i64_1,i64_2
! integer :: ibit

integer*8 :: arr(2)

arr = TRANSFER(i128, arr)
i64_1=arr(1)
i64_2=arr(2)
!----
!i64_1=0
!i64_2=0
!do ibit = 0,63
!    if(btest(i128,ibit)) i64_1=ibset(i64_1,ibit)
!enddo
!do ibit = 64,127
!    if(btest(i128,ibit)) i64_2=ibset(i64_2,ibit-64)
!enddo
end subroutine conv_128_to_64
!------------------------------------------------------

subroutine conv_64_to_128(i128,i64_1,i64_2)
integer*16, intent(out) :: i128
integer*8, intent(in) :: i64_1,i64_2
integer :: ibit

! integer*16 :: lowbits

!i128=0
!lowbits=0
!i128 = i64_2
!i128=ISHFT(i128,63)
!lowbits = i64_1
!i128=IOR(i128,lowbits)
!---
i128=0
do ibit = 0,63
    if(btest(i64_1,ibit)) i128=ibset(i128,ibit)
enddo
do ibit = 0,63
    if(btest(i64_2,ibit)) i128=ibset(i128,ibit+64)
enddo
end subroutine conv_64_to_128
!------------------------------------------------------

! These are same as above, but for when ik = i8b
subroutine conv_064_to_64(i128,i64_1,i64_2)
integer*8, intent(in) :: i128
integer*8, intent(out) :: i64_1,i64_2
! integer :: ibit

integer*8 :: arr(2)

arr = TRANSFER(i128, arr)
i64_1=arr(1)
i64_2=arr(2)
!----
!i64_1=0
!i64_2=0
!do ibit = 0,63
!    if(btest(i128,ibit)) i64_1=ibset(i64_1,ibit)
!enddo
!do ibit = 64,127
!    if(btest(i128,ibit)) i64_2=ibset(i64_2,ibit-64)
!enddo
end subroutine conv_064_to_64
!------------------------------------------------------

subroutine conv_64_to_064(i128,i64_1,i64_2)
integer*8, intent(out) :: i128
integer*8, intent(in) :: i64_1,i64_2
integer :: ibit

! integer*16 :: lowbits

!i128=0
!lowbits=0
!i128 = i64_2
!i128=ISHFT(i128,63)
!lowbits = i64_1
!i128=IOR(i128,lowbits)
!---
i128=0
do ibit = 0,63
    if(btest(i64_1,ibit)) i128=ibset(i128,ibit)
enddo
do ibit = 0,63
    if(btest(i64_2,ibit)) i128=ibset(i128,ibit+64)
enddo
end subroutine conv_64_to_064
!------------------------------------------------------

subroutine cluster_init

  character(len=16) filename

#ifdef MPI
  call MPI_INIT(mpierr)

  if(mpierr /= 0) then
    write(6,*) "MPI init failed!! errorcode:",mpierr
    stop "MPI init failed!! errorcode:"
  endif

  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, ncores, mpierr )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, whoami, mpierr )

  !Set up a communicator for each node
  CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,sharedComm,mpierr)
  !Figure out the rank for each processor, within the node
  CALL MPI_COMM_RANK(sharedComm,whoaminode,mpierr)


!**setup sending process
  allocate(snd_table(ncores*ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate snd_table'
  allocate(snd_displs(ncores*ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate snd_displs'
  allocate(snd_dsp(ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate snd_dsp'
  allocate(rcv_cnt(ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate rcv_cnt'

  allocate(recv_displs(ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate recv_displs'

  call init_mpitype_walker !** to define MPI_WALK_TYPE

  master_core=.false.
  if(whoami==0) master_core=.true.

  !Set master_core_node to true if you are the mastercore for the node
  master_core_node=.false.
  if(whoaminode==0) master_core_node=.true.
!** added redirection of stdout for all but master_core [09/12/13]
  if(.not.master_core) then
    close(6)
! Uncomment the following lines and comment the /dev/null line if you want output from slave processes
    if(whoami.le.9) then
      write(filename,'(''slave.'',i1)') whoami
    elseif(whoami.le.99) then
      write(filename,'(''slave.'',i2)') whoami
    elseif(whoami.le.999) then
      write(filename,'(''slave.'',i3)') whoami
    elseif(whoami.le.9999) then
      write(filename,'(''slave.'',i4)') whoami
    elseif(whoami.le.99999) then
      write(filename,'(''slave.'',i4)') whoami
    endif
   !open(6,file=filename)
    open(6,file='/dev/null')
  endif

  win_num = 0

  write(6,*) "MPI is ENABLED [nc=",ncores,"]!"
#else
  whoami = 0
  master_core = .true.
  whoaminode = 0
  master_core_node = .true.
  ncores = 1
  write(6,*) "MPI is DISABLED!"
#endif
  allocate(snd_cnt(ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate snd_cnt'
end subroutine cluster_init
!------------------------------------------------------

subroutine cluster_finalize
  integer :: i
#ifdef MPI
  !free the window
  do i=1,win_num
    call MPI_WIN_FREE(shmemWindow(i),mpierr)
  enddo
  call MPI_FINALIZE(mpierr)
  if(mpierr /= 0) then
    write(6,*) "MPI finalize failed!! errorcode:",mpierr
    stop "MPI finalize failed!! errorcode:"
  endif

!  deallocate(snd_table,snd_displs,recv_displs)
#endif
end subroutine cluster_finalize
!------------------------------------------------------
subroutine mpi_barr()
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_barr
!--------------------------------------------------------------------------
!--- MPI_ALLTOALLV ---------------------------------------------------------

subroutine mpi_alltoallv_iwalk(sndbuf,sndcounts,sdispls,recvbuf,recvcounts,rdispls)
type(t_walk), intent(in) :: sndbuf(:)
integer, intent(in) :: sndcounts(:)
integer, intent(in) :: sdispls(:),rdispls(:)
type(t_walk), intent(inout) :: recvbuf(:)
integer, intent(inout) :: recvcounts(:)
#ifdef MPI
call MPI_ALLTOALLV(sndbuf,sndcounts,sdispls,MPI_WALK_TYPE,recvbuf,recvcounts,rdispls,MPI_WALK_TYPE,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_alltoallv_iwalk
!--------------------------------------------------------------------------

!--- MPI_AGATH ---------------------------------------------------------
subroutine mpi_agath_int(sndbuf,sndcount,recvbuf,recvcount)
integer, intent(in) :: sndbuf(:)
integer, intent(in) :: sndcount
integer, intent(inout) :: recvbuf(:)
integer, intent(inout) :: recvcount
#ifdef MPI
call MPI_ALLGATHER(sndbuf,sndcount,MPI_INTEGER,recvbuf,recvcount,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_agath_int
!--------------------------------------------------------------------------

!--- MPI_AGATHV ---------------------------------------------------------
subroutine mpi_agathv_int(sndbuf,sndcount,recvbuf,recvcounts,displs)
    integer, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    integer, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)
#ifdef MPI
    call MPI_ALLGATHERV(sndbuf,sndcount,MPI_INTEGER,recvbuf,recvcounts,displs,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_agathv_int
!------------------------------------------------------

subroutine mpi_agathv_real_sp(sndbuf,sndcount,recvbuf,recvcounts,displs)
    real, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    real, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)
#ifdef MPI
    call MPI_ALLGATHERV(sndbuf,sndcount,MPI_REAL,recvbuf,recvcounts,displs,MPI_REAL,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_agathv_real_sp
!------------------------------------------------------

subroutine mpi_agathv_real_dp(sndbuf,sndcount,recvbuf,recvcounts,displs)
    real*8, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    real*8, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)
#ifdef MPI
    call MPI_ALLGATHERV(sndbuf,sndcount,MPI_REAL8,recvbuf,recvcounts,displs,MPI_REAL8,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_agathv_real_dp
!------------------------------------------------------

subroutine mpi_agathv_int8(sndbuf,sndcount,recvbuf,recvcounts,displs)
    integer*1, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    integer*1, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)
#ifdef MPI
    call MPI_ALLGATHERV(sndbuf,sndcount,MPI_INTEGER1,recvbuf,recvcounts,displs,MPI_INTEGER1,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_agathv_int8
!------------------------------------------------------

subroutine mpi_agathv_ik_vec(sndbuf,sndcount,recvbuf,recvcounts,displs)
    type(ik_vec), intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    type(ik_vec), intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(recvcounts(1:ncores))
    allocate(high(tot_snd))
    allocate(low(tot_snd))
    allocate(buff(tot_snd))

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii)%v(1),low(ii),high(ii))
    enddo
    call MPI_ALLGATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_ALLGATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii)%v(1),low(ii),high(ii))
    enddo

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii)%v(2),low(ii),high(ii))
    enddo
    call MPI_ALLGATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_ALLGATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii)%v(2),low(ii),high(ii))
    enddo

    deallocate(high,low,buff)
#endif
end subroutine mpi_agathv_ik_vec
!------------------------------------------------------

subroutine mpi_agathv_i128(sndbuf,sndcount,recvbuf,recvcounts,displs)
    integer*16, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    integer*16, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(recvcounts(1:ncores))
    allocate(high(tot_snd))
    allocate(low(tot_snd))
    allocate(buff(tot_snd))

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii),low(ii),high(ii))
    enddo
    call MPI_ALLGATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_ALLGATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii),low(ii),high(ii))
    enddo

    deallocate(high,low,buff)
#endif
end subroutine mpi_agathv_i128
!--------------------------------------------------------------------------

!--- MPI_GATH ---------------------------------------------------------
subroutine mpi_gath_int(sndbuf,sndcount,recvbuf,recvcount,root)
    integer, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    integer, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcount
#ifdef MPI
    call MPI_GATHER(sndbuf,sndcount,MPI_INTEGER,recvbuf,recvcount,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
#else
    recvbuf(1:recvcount)=sndbuf(1:recvcount)
#endif
end subroutine mpi_gath_int
!------------------------------------------------------

subroutine mpi_gath_real_sp(sndbuf,sndcount,recvbuf,recvcount,root)
    real, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    real, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcount
#ifdef MPI
    call MPI_GATHER(sndbuf,sndcount,MPI_REAL,recvbuf,recvcount,MPI_REAL,root,MPI_COMM_WORLD,mpierr)
#else
    recvbuf(1:recvcount)=sndbuf(1:recvcount)
#endif
end subroutine mpi_gath_real_sp
!------------------------------------------------------

subroutine mpi_gath_real_dp(sndbuf,sndcount,recvbuf,recvcount,root)
    real*8, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    real*8, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcount
#ifdef MPI
    call MPI_GATHER(sndbuf,sndcount,MPI_REAL8,recvbuf,recvcount,MPI_REAL8,root,MPI_COMM_WORLD,mpierr)
#else
    recvbuf(1:recvcount)=sndbuf(1:recvcount)
#endif
end subroutine mpi_gath_real_dp
!------------------------------------------------------

subroutine mpi_gath_int8(sndbuf,sndcount,recvbuf,recvcount,root)
    integer*1, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    integer*1, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcount
#ifdef MPI
    call MPI_GATHER(sndbuf,sndcount,MPI_INTEGER1,recvbuf,recvcount,MPI_INTEGER1,root,MPI_COMM_WORLD,mpierr)
#else
    recvbuf(1:recvcount)=sndbuf(1:recvcount)
#endif
end subroutine mpi_gath_int8
!------------------------------------------------------

!subroutine mpi_gath_ik_vec(sndbuf,sndcount,recvbuf,recvcount,root)
!    type(ik_vec), intent(in) :: sndbuf(:)
!    integer, intent(in) :: sndcount,root
!    type(ik_vec), intent(out) :: recvbuf(:)
!    integer, intent(in) :: recvcount
!
!    integer*8,allocatable  :: high(:),low(:),buff(:)
!    integer :: ii,tot_snd
!#ifdef MPI
!    tot_snd=sum(recvcount(1:ncores))
!    allocate(high(tot_snd))
!    allocate(low(tot_snd))
!    allocate(buff(tot_snd))
!
!    do ii=1,sndcount
!        call conv_128_to_64(sndbuf(ii)%v(1),low(ii),high(ii))
!    enddo
!    call MPI_GATHER(low,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    low=buff !***now buff contains the full vector
!    call MPI_GATHER(high,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    high=buff
!    do ii=1,tot_snd
!        call conv_64_to_128(recvbuf(ii)%v(1),low(ii),high(ii))
!    enddo
!
!    do ii=1,sndcount
!        call conv_128_to_64(sndbuf(ii)%v(2),low(ii),high(ii))
!    enddo
!    call MPI_GATHER(low,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    low=buff !***now buff contains the full vector
!    call MPI_GATHER(high,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    high=buff
!    do ii=1,tot_snd
!        call conv_64_to_128(recvbuf(ii)%v(2),low(ii),high(ii))
!    enddo
!    deallocate(high,low,buff)
!#endif
!end subroutine mpi_gath_ik_vec
!!------------------------------------------------------
!
!subroutine mpi_gath_i128(sndbuf,sndcount,recvbuf,recvcount,root)
!    integer*16, intent(in) :: sndbuf(:)
!    integer, intent(in) :: sndcount,root
!    integer*16, intent(out) :: recvbuf(:)
!    integer, intent(in) :: recvcount
!
!    integer*8,allocatable  :: high(:),low(:),buff(:)
!    integer :: ii,tot_snd
!#ifdef MPI
!    tot_snd=sum(recvcount(1:ncores))
!    allocate(high(tot_snd))
!    allocate(low(tot_snd))
!    allocate(buff(tot_snd))
!
!    do ii=1,sndcount
!        call conv_128_to_64(sndbuf(ii),low(ii),high(ii))
!    enddo
!    call MPI_GATHER(low,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    low=buff !***now buff contains the full vector
!    call MPI_GATHER(high,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    high=buff
!    do ii=1,tot_snd
!        call conv_64_to_128(recvbuf(ii),low(ii),high(ii))
!    enddo
!
!    deallocate(high,low,buff)
!#endif
!end subroutine mpi_gath_i128
!!--------------------------------------------------------------------------

!--- MPI_GATHV ---------------------------------------------------------
subroutine mpi_gathv_int(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    integer, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    integer, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)
#ifdef MPI
    call MPI_GATHERV(sndbuf,sndcount,MPI_INTEGER,recvbuf,recvcounts,displs,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_gathv_int
!------------------------------------------------------

subroutine mpi_gathv_real_sp(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    real, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    real, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)
#ifdef MPI
    call MPI_GATHERV(sndbuf,sndcount,MPI_REAL,recvbuf,recvcounts,displs,MPI_REAL,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_gathv_real_sp
!------------------------------------------------------

subroutine mpi_gathv_real_dp(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    real*8, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    real*8, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)
#ifdef MPI
    call MPI_GATHERV(sndbuf,sndcount,MPI_REAL8,recvbuf,recvcounts,displs,MPI_REAL8,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_gathv_real_dp
!------------------------------------------------------

subroutine mpi_gathv_int8(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    integer*1, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    integer*1, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)
#ifdef MPI
    call MPI_GATHERV(sndbuf,sndcount,MPI_INTEGER1,recvbuf,recvcounts,displs,MPI_INTEGER1,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_gathv_int8
!------------------------------------------------------

subroutine mpi_gathv_ik_vec(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    type(ik_vec), intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    type(ik_vec), intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(recvcounts(1:ncores))
    allocate(high(tot_snd))
    allocate(low(tot_snd))
    allocate(buff(tot_snd))

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii)%v(1),low(ii),high(ii))
    enddo
    call MPI_GATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_GATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii)%v(1),low(ii),high(ii))
    enddo

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii)%v(2),low(ii),high(ii))
    enddo
    call MPI_GATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_GATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii)%v(2),low(ii),high(ii))
    enddo
    deallocate(high,low,buff)
#endif
end subroutine mpi_gathv_ik_vec
!------------------------------------------------------

subroutine mpi_gathv_i128(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    integer*16, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    integer*16, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(recvcounts(1:ncores))
    allocate(high(tot_snd))
    allocate(low(tot_snd))
    allocate(buff(tot_snd))

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii),low(ii),high(ii))
    enddo
    call MPI_GATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_GATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii),low(ii),high(ii))
    enddo

    deallocate(high,low,buff)
#endif
end subroutine mpi_gathv_i128
!--------------------------------------------------------------------------

!--- MPI_SCATTV ---------------------------------------------------------
subroutine mpi_scattv_int(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    integer, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    integer, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount
#ifdef MPI
    call MPI_SCATTERV(sndbuf,sndcounts,displs,MPI_INTEGER,recvbuf,recvcount,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_scattv_int
!------------------------------------------------------

subroutine mpi_scattv_real_dp(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    real*8, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    real*8, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount
#ifdef MPI
    call MPI_SCATTERV(sndbuf,sndcounts,displs,MPI_REAL8,recvbuf,recvcount,MPI_REAL8,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_scattv_real_dp
!------------------------------------------------------

subroutine mpi_scattv_real_sp(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    real, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    real, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount
#ifdef MPI
    call MPI_SCATTERV(sndbuf,sndcounts,displs,MPI_REAL,recvbuf,recvcount,MPI_REAL,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_scattv_real_sp
!------------------------------------------------------

subroutine mpi_scattv_int8(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    integer*1, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    integer*1, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount
#ifdef MPI
    call MPI_SCATTERV(sndbuf,sndcounts,displs,MPI_INTEGER1,recvbuf,recvcount,MPI_INTEGER1,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_scattv_int8
!------------------------------------------------------

subroutine mpi_scattv_iwalk(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    type(t_walk), intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    type(t_walk), intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount
#ifdef MPI
    call MPI_SCATTERV(sndbuf,sndcounts,displs,MPI_WALK_TYPE,recvbuf,recvcount,MPI_WALK_TYPE,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_scattv_iwalk
!------------------------------------------------------

subroutine mpi_scattv_ik_vec(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    type(ik_vec), intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    type(ik_vec), intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(sndcounts(1:ncores))
    allocate(high(tot_snd))
    allocate(low(tot_snd))
    allocate(buff(tot_snd))

    do ii=1,tot_snd
        call conv_128_to_64(sndbuf(ii)%v(1),low(ii),high(ii))
    enddo
    call MPI_SCATTERV(low,sndcounts,displs,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_SCATTERV(high,sndcounts,displs,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,recvcount
        call conv_64_to_128(recvbuf(ii)%v(1),low(ii),high(ii))
    enddo

    do ii=1,tot_snd
        call conv_128_to_64(sndbuf(ii)%v(2),low(ii),high(ii))
    enddo
    call MPI_SCATTERV(low,sndcounts,displs,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_SCATTERV(high,sndcounts,displs,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,recvcount
        call conv_64_to_128(recvbuf(ii)%v(2),low(ii),high(ii))
    enddo
    deallocate(high,low,buff)
#endif
end subroutine mpi_scattv_ik_vec
!------------------------------------------------------

subroutine mpi_scattv_i128(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    integer*16, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    integer*16, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount

    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(sndcounts(1:ncores))

    if(.not.allocated(i128_low)) then
      allocate(i128_low(loc_mwalk),stat=istat)
      if(istat.ne.0) stop 'failed to allocate i128_low'
    endif
    if(.not.allocated(i128_high)) then
      allocate(i128_high(loc_mwalk),stat=istat)
      if(istat.ne.0) stop 'failed to allocate i128_high'
    endif
    if(.not.allocated(i128_buff)) then
      allocate(i128_buff(loc_mwalk),stat=istat)
      if(istat.ne.0) stop 'failed to allocate i128_buff'
    endif
    if(whoami==root) then
        do ii=1,tot_snd
            call conv_128_to_64(sndbuf(ii),i128_low(ii),i128_high(ii))
        enddo
    endif
    call MPI_SCATTERV(i128_low,sndcounts,displs,MPI_INTEGER8,i128_buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    i128_low(1:recvcount)=i128_buff(1:recvcount)
    call MPI_SCATTERV(i128_high,sndcounts,displs,MPI_INTEGER8,i128_buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    do ii=1,recvcount
        call conv_64_to_128(recvbuf(ii),i128_low(ii),i128_buff(ii))
    enddo
#endif
end subroutine mpi_scattv_i128
!--------------------------------------------------------------------------
!--- MPI_REDUCE for max --------------------------------------------------------
subroutine mpi_red_max_real_sp(rspvar_in,rspvar_out,root)
real, intent(in) :: rspvar_in
real, intent(out) :: rspvar_out
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rspvar_in,rspvar_out,1,MPI_REAL,MPI_MAX,root,MPI_COMM_WORLD,mpierr)
#else
  rspvar_out=rspvar_in
#endif
end subroutine mpi_red_max_real_sp
!------------------------------------------------------

subroutine mpi_red_max_real_dp(rdpvar_in,rdpvar_out,root)
real*8, intent(in) :: rdpvar_in
real*8, intent(out) :: rdpvar_out
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rdpvar_in,rdpvar_out,1,MPI_REAL8,MPI_MAX,root,MPI_COMM_WORLD,mpierr)
#else
  rdpvar_out=rdpvar_in
#endif
end subroutine mpi_red_max_real_dp
!------------------------------------------------------

subroutine mpi_red_max_real_dparray(rdparr_in,rdparr_out,root)
real*8, intent(in) :: rdparr_in(:)
real*8, intent(out) :: rdparr_out(:)
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rdparr_in,rdparr_out,size(rdparr_in),MPI_REAL8,MPI_MAX,root,MPI_COMM_WORLD,mpierr)
#else
  rdparr_out=rdparr_in
#endif
end subroutine mpi_red_max_real_dparray
!--------------------------------------------------------------------------

!--- MPI_REDUCE for min --------------------------------------------------------
subroutine mpi_red_min_real_sp(rspvar_in,rspvar_out,root)
real, intent(in) :: rspvar_in
real, intent(out) :: rspvar_out
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rspvar_in,rspvar_out,1,MPI_REAL,MPI_MIN,root,MPI_COMM_WORLD,mpierr)
#else
  rspvar_out=rspvar_in
#endif
end subroutine mpi_red_min_real_sp
!------------------------------------------------------

subroutine mpi_red_min_real_dp(rdpvar_in,rdpvar_out,root)
real*8, intent(in) :: rdpvar_in
real*8, intent(out) :: rdpvar_out
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rdpvar_in,rdpvar_out,1,MPI_REAL8,MPI_MIN,root,MPI_COMM_WORLD,mpierr)
#else
  rdpvar_out=rdpvar_in
#endif
end subroutine mpi_red_min_real_dp
!------------------------------------------------------

subroutine mpi_red_min_real_dparray(rdparr_in,rdparr_out,root)
real*8, intent(in) :: rdparr_in(:)
real*8, intent(out) :: rdparr_out(:)
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rdparr_in,rdparr_out,size(rdparr_in),MPI_REAL8,MPI_MIN,root,MPI_COMM_WORLD,mpierr)
#else
  rdparr_out=rdparr_in
#endif
end subroutine mpi_red_min_real_dparray
!--------------------------------------------------------------------------

!--- MPI_REDUCE for sum --------------------------------------------------------
subroutine mpi_red_real_sp(rspvar,root)
real, intent(inout) :: rspvar
integer, intent(in) :: root
real :: rspsum
#ifdef MPI
  call MPI_REDUCE(rspvar,rspsum,1,MPI_REAL,MPI_SUM,root,MPI_COMM_WORLD,mpierr)
  if(whoami==root) rspvar=rspsum
#endif
end subroutine mpi_red_real_sp
!------------------------------------------------------

subroutine mpi_red_real_dp(rdpvar,root)
real*8, intent(inout) :: rdpvar
integer, intent(in) :: root
real*8 :: rdpsum
#ifdef MPI
  call MPI_REDUCE(rdpvar,rdpsum,1,MPI_REAL8,MPI_SUM,root,MPI_COMM_WORLD,mpierr)
  if(whoami==root) rdpvar=rdpsum
#endif
end subroutine mpi_red_real_dp
!------------------------------------------------------

subroutine mpi_red_real_dparray(rdparr,root)
real*8, intent(inout) :: rdparr(:)
integer, intent(in) :: root
real*8 :: rdpsum(size(rdparr))
#ifdef MPI
  call MPI_REDUCE(rdparr,rdpsum,size(rdparr),MPI_REAL8,MPI_SUM,root,MPI_COMM_WORLD,mpierr)
  if(whoami==root) rdparr=rdpsum
#endif
end subroutine mpi_red_real_dparray
!--------------------------------------------------------------------------

!--- MPI_REDUCE_SCATTER for sum ---------------------------------------------------
subroutine mpi_redscatt_real_dparray(rdparr,counts)
real*8, intent(inout) :: rdparr(:)
integer, intent(in) :: counts(:)
real*8  :: rdpsum(counts(whoami+1))
#ifdef MPI
  call MPI_REDUCE_SCATTER(rdparr,rdpsum,counts, MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
  rdparr(1:counts(whoami+1))=rdpsum(1:counts(whoami+1))
#endif
end subroutine mpi_redscatt_real_dparray
!--------------------------------------------------------------------------

!--- MPI_ALLREDUCE --------------------------------------------------------
subroutine mpi_allred_int(ivar)
integer, intent(inout) :: ivar
integer :: isum
#ifdef MPI
  call MPI_ALLREDUCE(ivar,isum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
  ivar=isum
#endif
end subroutine mpi_allred_int
!------------------------------------------------------

subroutine mpi_allred_real_sp(rspvar)
real, intent(inout) :: rspvar
real :: rspsum
#ifdef MPI
  call MPI_ALLREDUCE(rspvar,rspsum,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,mpierr)
  rspvar=rspsum
#endif
end subroutine mpi_allred_real_sp
!------------------------------------------------------

subroutine mpi_allred_real_dp(rdpvar)
real*8, intent(inout) :: rdpvar
real*8 :: rdpsum
#ifdef MPI
  call MPI_ALLREDUCE(rdpvar,rdpsum,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
  rdpvar=rdpsum
#endif
end subroutine mpi_allred_real_dp

subroutine mpi_allred_real_dparray(rdparr)
real*8, intent(inout) :: rdparr(:)
real*8 :: rdpsum(size(rdparr))
#ifdef MPI
  call MPI_ALLREDUCE(rdparr,rdpsum,size(rdparr),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
  rdparr=rdpsum
#endif
end subroutine mpi_allred_real_dparray
!------------------------------------------------------

subroutine mpi_allred_iarray(iarr)
integer, intent(inout) :: iarr(:)
integer :: isum(size(iarr))
#ifdef MPI
  call MPI_ALLREDUCE(iarr,isum,size(iarr),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
  iarr=isum
#endif
end subroutine mpi_allred_iarray
!--------------------------------------------------------------------------

!--- MPI_BSEND_DIFFROOT ------------------------------------------------------
subroutine mpi_bsend_diffroot_int(ivar,id_root)
integer, intent(inout) :: ivar
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(ivar,1,MPI_INTEGER,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_int
!------------------------------------------------------

subroutine mpi_bsend_diffroot_iarray(iarr,id_root)
integer,intent(inout) :: iarr(:)
integer, intent(in) :: id_root
#ifdef MPI
call MPI_Bcast(iarr,size(iarr),MPI_INTEGER,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_iarray
!------------------------------------------------------

subroutine mpi_bsend_diffroot_real_sp(rspvar,id_root)
real, intent(inout) :: rspvar
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(rspvar,1,MPI_REAL,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_real_sp
!------------------------------------------------------

subroutine mpi_bsend_diffroot_real_dp(rdpvar,id_root)
real*8, intent(inout) :: rdpvar
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(rdpvar,1,MPI_REAL8,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_real_dp
!------------------------------------------------------

subroutine mpi_bsend_diffroot_logical(lvar,id_root)
logical, intent(inout) :: lvar
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(lvar,1,MPI_LOGICAL,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_logical
!------------------------------------------------------

subroutine mpi_bsend_diffroot_int64(ivar64, id_root)
integer*8, intent(inout) :: ivar64
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(ivar64,1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_int64
!------------------------------------------------------

subroutine mpi_bsend_diffroot_sndlist(isndlist, id_root)
type(t_sndlist), intent(inout) :: isndlist
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(isndlist%to,1,MPI_INTEGER,id_root,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(isndlist%wlkid,1,MPI_INTEGER,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_sndlist
!------------------------------------------------------

subroutine mpi_bsend_diffroot_int128(ivar128, id_root)
integer*16, intent(inout) :: ivar128
integer, intent(in) :: id_root
integer*8  :: high,low
#ifdef MPI

  call conv_128_to_64(ivar128,low,high)
  call MPI_Bcast(low,1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(high,1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
  call conv_64_to_128(ivar128,low,high)

#endif
end subroutine mpi_bsend_diffroot_int128
!------------------------------------------------------

subroutine mpi_bsend_diffroot_int8(ivar8, id_root)
integer*1, intent(inout) :: ivar8
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(ivar8,1,MPI_INTEGER1,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_int8
!------------------------------------------------------

subroutine mpi_bsend_diffroot_ik_vec(ik_vecvar,id_root)
type(ik_vec), intent(inout) :: ik_vecvar
integer, intent(in) :: id_root
integer*8  :: high,low
integer :: i
#ifdef MPI

  if (ik==selected_int_kind(38)) then
    do i=1,num_words
      call conv_128_to_64(ik_vecvar%v(i),low,high)
      call MPI_Bcast(low,1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(high,1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
      call conv_64_to_128(ik_vecvar%v(i),low,high)
    enddo
  else
    do i=1,num_words
      call MPI_Bcast(ik_vecvar%v(i),1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
    enddo
  endif

#endif
end subroutine mpi_bsend_diffroot_ik_vec
!--------------------------------------------------------------------------

!--- MPI_BSEND ------------------------------------------------------------
subroutine mpi_bsend_int(ivar)
integer, intent(inout) :: ivar
#ifdef MPI
  call MPI_Bcast(ivar,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_int
!------------------------------------------------------

subroutine mpi_bsend_real_sp(rspvar)
real, intent(inout) :: rspvar
#ifdef MPI
  call MPI_Bcast(rspvar,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_real_sp
!------------------------------------------------------

subroutine mpi_bsend_real_dp(rdpvar)
real*8, intent(inout) :: rdpvar
#ifdef MPI
  call MPI_Bcast(rdpvar,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_real_dp
!------------------------------------------------------

subroutine mpi_bsend_real_dparray(rdparr)
real*8, intent(inout) :: rdparr(:)
#ifdef MPI
  call MPI_Bcast(rdparr,size(rdparr),MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_real_dparray
!------------------------------------------------------

subroutine mpi_bsend_complex(cvar)
complex, intent(inout) :: cvar
#ifdef MPI
  call MPI_Bcast(cvar,1,MPI_COMPLEX,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_complex
!------------------------------------------------------

subroutine mpi_bsend_string(svar)
character, intent(inout) :: svar*16
#ifdef MPI
 call MPI_Bcast(svar,len(svar),MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_string
!------------------------------------------------------

subroutine mpi_bsend_iarray(iarr)
integer,intent(inout) :: iarr(:)
#ifdef MPI
call MPI_Bcast(iarr,size(iarr),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_iarray
!------------------------------------------------------

!ad-hoc for rand seed, should be substitute with one for matrices
subroutine mpi_bsend_irand_seed(irand)
integer,intent(inout),dimension(4,2) :: irand
#ifdef MPI
call MPI_Bcast(irand(:,1),4,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
call MPI_Bcast(irand(:,2),4,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_irand_seed
!------------------------------------------------------

subroutine mpi_bsend_logical(lvar)
logical, intent(inout) :: lvar
#ifdef MPI
  call MPI_Bcast(lvar,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_logical
!------------------------------------------------------

subroutine mpi_bsend_int64(ivar64)
integer*8, intent(inout) :: ivar64
#ifdef MPI
  call MPI_Bcast(ivar64,1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_int64
!------------------------------------------------------

subroutine mpi_bsend_i64arr(i64arr)
integer*8, intent(inout) :: i64arr(:)
#ifdef MPI
    call MPI_Bcast(i64arr,size(i64arr),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_i64arr
!------------------------------------------------------

subroutine mpi_bsend_i128arr(i128arr)
integer*16, intent(inout) :: i128arr(:)

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii
#ifdef MPI
    allocate(high(size(i128arr)))
    allocate(low(size(i128arr)))
    allocate(buff(size(i128arr)))
    if(master_core) then
        do ii=1,size(i128arr)
            call conv_128_to_64(i128arr(ii),low(ii),high(ii))
        enddo
    endif
    call MPI_Bcast(low,size(low),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(high,size(high),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    do ii=1,size(i128arr)
        call conv_64_to_128(i128arr(ii),low(ii),high(ii))
    enddo
    deallocate(high,low,buff)
#endif
end subroutine mpi_bsend_i128arr
!------------------------------------------------------

subroutine mpi_bsend_ikvec(ikv)
    type(ik_vec), intent(inout) :: ikv
    integer*8,allocatable  :: high(:),low(:)
    integer :: ii
#ifdef MPI
  if (ik==selected_int_kind(38)) then
    allocate(high(num_words))
    allocate(low(num_words))
    if(master_core) then
        do ii=1,num_words
            call conv_128_to_64(ikv%v(ii),low(ii),high(ii))
        enddo
    endif
    call MPI_Bcast(low,size(low),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(high,size(high),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    do ii=1,num_words
        call conv_64_to_128(ikv%v(ii),low(ii),high(ii))
    enddo
    deallocate(high,low)
  else
    do ii=1,num_words
      call MPI_Bcast(ikv%v(ii),1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    enddo
  endif

#endif
end subroutine mpi_bsend_ikvec
!------------------------------------------------------

subroutine mpi_bsend_ikvec_arr(ikvarr)
  type(ik_vec), intent(inout) :: ikvarr(:)
  integer*8,allocatable  :: high(:),low(:)
  integer :: ii,jj,kk
#ifdef MPI
  if (ik==selected_int_kind(38)) then
    allocate(high(num_words*size(ikvarr)))
    allocate(low(num_words*size(ikvarr)))
    if(master_core) then
        kk=1
        do jj=1,size(ikvarr)
            do ii=1,num_words
                call conv_128_to_64(ikvarr(jj)%v(ii),low(kk),high(kk))
                kk=kk+1
            enddo
        enddo
    endif
    call MPI_Bcast(low,size(low),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(high,size(high),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    kk=1
    do jj=1,size(ikvarr)
        do ii=1,num_words
            call conv_64_to_128(ikvarr(jj)%v(ii),low(kk),high(kk))
            kk=kk+1
        enddo
    enddo
    deallocate(high,low)
  else
    do jj=1,size(ikvarr)
      do ii=1,num_words
        call MPI_Bcast(ikvarr(jj)%v(ii),1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
      enddo
    enddo
  endif

#endif
end subroutine mpi_bsend_ikvec_arr
!-------------------------------------------------------------------------

!--- MPI_SND ------------------------------------------------------------
!-- assuming master_core is the sender!!
!-------------------------------------------------------------------------
subroutine mpi_snd_int(ivar,target_core)
integer, intent(inout) :: ivar
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
  if(master_core) call mpi_send(ivar,1,MPI_INTEGER,target_core,69,MPI_COMM_WORLD,mpierr)
  if(whoami==target_core) call mpi_recv(ivar,1,MPI_INTEGER,0,69,MPI_COMM_WORLD,cstatus,mpierr)
#endif
end subroutine mpi_snd_int
!------------------------------------------------------

subroutine mpi_snd_real_sp(rspvar,target_core)
real, intent(inout) :: rspvar
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
  if(master_core) call mpi_send(rspvar,1,MPI_REAL,target_core,69,MPI_COMM_WORLD,mpierr)
  if(whoami==target_core) call mpi_recv(rspvar,1,MPI_REAL,0,69,MPI_COMM_WORLD,cstatus,mpierr)
#endif
end subroutine mpi_snd_real_sp
!------------------------------------------------------

subroutine mpi_snd_real_dp(rdpvar,target_core)
real*8, intent(inout) :: rdpvar
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
  if(master_core) call mpi_send(rdpvar,1,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
  if(whoami==target_core) call mpi_recv(rdpvar,1,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
#endif
end subroutine mpi_snd_real_dp
!------------------------------------------------------

subroutine mpi_snd_real_dparray(rdparr,target_core)
real*8, intent(inout) :: rdparr(:)
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
integer :: total_size,istart,iend,i
  total_size=size(rdparr)
  if (total_size<=mpi_send_limit) then
    if(master_core) call mpi_send(rdparr,total_size,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
    if(whoami==target_core) call mpi_recv(rdparr,total_size,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
  else
    do i=1,(total_size-1)/mpi_send_limit+1
      istart=(i-1)*mpi_send_limit+1
      iend=min(total_size,i*mpi_send_limit)
     !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
      if(master_core) call mpi_send(rdparr(istart:iend),iend-istart+1,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
      if(whoami==target_core) call mpi_recv(rdparr(istart:iend),iend-istart+1,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
    enddo
  endif
#endif
end subroutine mpi_snd_real_dparray
!------------------------------------------------------

subroutine mpi_snd_iarray(iarr,target_core)
integer,intent(inout) :: iarr(:)
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
integer :: total_size,istart,iend,i
  total_size=size(iarr)
  if (total_size<=mpi_send_limit) then
    if(master_core) call mpi_send(iarr,total_size,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
    if(whoami==target_core) call mpi_recv(iarr,total_size,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
  else
    do i=1,(total_size-1)/mpi_send_limit+1
      istart=(i-1)*mpi_send_limit+1
      iend=min(total_size,i*mpi_send_limit)
     !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
      if(master_core) call mpi_send(iarr(istart:iend),iend-istart+1,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
      if(whoami==target_core) call mpi_recv(iarr(istart:iend),iend-istart+1,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
    enddo
  endif
#endif
end subroutine mpi_snd_iarray
!------------------------------------------------------

subroutine mpi_snd_i64array(i64arr,target_core)
integer*8,intent(inout) :: i64arr(:)
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
integer :: total_size,istart,iend,i
  total_size=size(i64arr)
  if (total_size<=mpi_send_limit) then
    if(master_core) call mpi_send(i64arr,total_size,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
    if(whoami==target_core) call mpi_recv(i64arr,total_size,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
  else
    do i=1,(total_size-1)/mpi_send_limit+1
      istart=(i-1)*mpi_send_limit+1
      iend=min(total_size,i*mpi_send_limit)
     !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
      if(master_core) call mpi_send(i64arr(istart:iend),iend-istart+1,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
      if(whoami==target_core) call mpi_recv(i64arr(istart:iend),iend-istart+1,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
    enddo
  endif
#endif
end subroutine mpi_snd_i64array
!------------------------------------------------------

subroutine mpi_snd_int64(int64,target_core)
integer*8,intent(inout) :: int64
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
  if(master_core) call mpi_send(int64,1,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
  if(whoami==target_core) call mpi_recv(int64,1,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
#endif
end subroutine mpi_snd_int64

subroutine mpi_snd_i128arr(i128arr,target_core)
    integer*16, intent(inout) :: i128arr(:)
    integer, intent(in) :: target_core

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
integer :: total_size,istart,iend,i
  total_size=size(i128arr)
  if(master_core .or. (whoami==target_core)) then
    if (total_size<=mpi_send_limit) then
      allocate(high(total_size))
      allocate(low(total_size))
      allocate(buff(total_size))
      if(master_core) then
          do ii=1,total_size
              call conv_128_to_64(i128arr(ii),low(ii),high(ii))
          enddo
          call mpi_send(low,total_size,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
          call mpi_send(high,total_size,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
      endif
      if(whoami == target_core) then
          call mpi_recv(low,total_size,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
          call mpi_recv(high,total_size,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
          do ii=1,total_size
              call conv_64_to_128(i128arr(ii),low(ii),high(ii))
          enddo
      endif
    else
      allocate(high(mpi_send_limit))
      allocate(low(mpi_send_limit))
      allocate(buff(mpi_send_limit))
      do i=1,(total_size-1)/mpi_send_limit+1
        istart=(i-1)*mpi_send_limit+1
        iend=min(total_size,i*mpi_send_limit)
       !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
        if(master_core) then
            do ii=istart,iend
                call conv_128_to_64(i128arr(ii),low(ii-istart+1),high(ii-istart+1))
            enddo
            call mpi_send(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
            call mpi_send(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
        endif
        if(whoami == target_core) then
            call mpi_recv(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
            call mpi_recv(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
            do ii=istart,iend
                call conv_64_to_128(i128arr(ii),low(ii-istart+1),high(ii-istart+1))
            enddo
        endif
      enddo
    endif
    deallocate(high,low,buff)
  endif
#endif
end subroutine mpi_snd_i128arr
!------------------------------------------------------

subroutine mpi_snd_ikvec(ikv,target_core)
    type(ik_vec), intent(inout) :: ikv
    integer, intent(in) :: target_core

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii
#ifdef MPI
  integer :: cstatus(MPI_STATUS_SIZE)
  integer :: total_size,istart,iend,i
  total_size=num_words
  if(master_core .or. (whoami==target_core)) then
    if (total_size<=mpi_send_limit) then
      allocate(high(total_size))
      allocate(low(total_size))
      allocate(buff(total_size))
      if(master_core) then
          do ii=1,total_size
              call conv_128_to_64(ikv%v(ii),low(ii),high(ii))
          enddo
          call mpi_send(low,total_size,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
          call mpi_send(high,total_size,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
      endif
      if(whoami == target_core) then
          call mpi_recv(low,total_size,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
          call mpi_recv(high,total_size,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
          do ii=1,total_size
              call conv_64_to_128(ikv%v(ii),low(ii),high(ii))
          enddo
      endif
    else
      allocate(high(mpi_send_limit))
      allocate(low(mpi_send_limit))
      allocate(buff(mpi_send_limit))
      do i=1,(total_size-1)/mpi_send_limit+1
        istart=(i-1)*mpi_send_limit+1
        iend=min(total_size,i*mpi_send_limit)
       !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
        if(master_core) then
            do ii=istart,iend
                call conv_128_to_64(ikv%v(ii),low(ii-istart+1),high(ii-istart+1))
            enddo
            call mpi_send(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
            call mpi_send(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
        endif
        if(whoami == target_core) then
            call mpi_recv(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
            call mpi_recv(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
            do ii=istart,iend
                call conv_64_to_128(ikv%v(ii),low(ii-istart+1),high(ii-istart+1))
            enddo
        endif
      enddo
    endif
    deallocate(high,low,buff)
  endif
#endif
end subroutine mpi_snd_ikvec
!------------------------------------------------------

subroutine mpi_snd_ikvec_array(ikvarr,target_core)
    type(ik_vec), intent(inout) :: ikvarr(:)
    integer, intent(in) :: target_core

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,jj,kk
#ifdef MPI
  integer :: cstatus(MPI_STATUS_SIZE)
  integer :: total_size,istart,iend,i
  total_size=num_words*size(ikvarr)
  if(master_core .or. (whoami==target_core)) then
    if (total_size<=mpi_send_limit) then
      allocate(high(total_size))
      allocate(low(total_size))
      allocate(buff(total_size))
      if(master_core) then
          kk=1
          do jj=1,size(ikvarr)
              do ii=1,num_words
                  call conv_128_to_64(ikvarr(jj)%v(ii),low(kk),high(kk))
                  kk=kk+1
              enddo
          enddo
          call mpi_send(low,total_size,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
          call mpi_send(high,total_size,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
      endif
      if(whoami == target_core) then
          call mpi_recv(low,total_size,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
          call mpi_recv(high,total_size,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
          kk=1
          do jj=1,size(ikvarr)
              do ii=1,num_words
                  call conv_64_to_128(ikvarr(jj)%v(ii),low(kk),high(kk))
                  kk=kk+1
              enddo
          enddo
      endif
    else
      allocate(high(mpi_send_limit))
      allocate(low(mpi_send_limit))
      allocate(buff(mpi_send_limit))
      do i=1,(total_size-1)/mpi_send_limit+1
        istart=(i-1)*mpi_send_limit+1
        iend=min(total_size,i*mpi_send_limit)
       !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
        if(master_core) then
            do ii=istart,iend
                call conv_128_to_64(ikvarr(int((ii-1)/num_words)+1)%v(mod(ii-1,num_words)+1),low(ii-istart+1),high(ii-istart+1))
            enddo
            call mpi_send(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
            call mpi_send(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
        endif
        if(whoami == target_core) then
            call mpi_recv(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
            call mpi_recv(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
            do ii=istart,iend
                call conv_64_to_128(ikvarr(int((ii-1)/num_words)+1)%v(mod(ii-1,num_words)+1),low(ii-istart+1),high(ii-istart+1))
            enddo
        endif
      enddo
    endif
    deallocate(high,low,buff)
  endif
#endif
end subroutine mpi_snd_ikvec_array
!------------------------------------------------------------------------------

!**WARNING: may break consistency if we change the kind of some of the walker arrays, should work with a generic object walker [TODO]
!my_cnt: number of walks to be sent
!my_nwalk: current population on core, used as offset for addresses [buf_snd=buf_snd(my_nwalk+1:my_nwalk+my_cnt)]
!nwalk: total population, used for current get_owner [CHANGE WHEN USING HASH_TABLE]
!my_nwalk_1 : first element of my new population [shall be set to zero when everyone works with its own walks, just use my_nwalk instead]
! Only used once at beginning of run

subroutine mpi_sendwalks(isnd_start,my_nwalk,my_cnt,my_nwalk_1,initiator,e_num,e_den,diag_elems)
  integer,intent(in) :: isnd_start,my_nwalk_1,my_cnt
  integer(i1b), intent(inout) :: initiator(:)
  real(rk),intent(inout) :: e_num(:),e_den(:),diag_elems(:)

  integer,intent(inout) :: my_nwalk

  integer :: isnd,iown,tot_snd,icore,my_ic,my_own,iw
#ifdef NUM_ORBITALS_GT_127
  integer :: islice
#endif

#ifdef MPI
  snd_table=0
  do isnd=1,my_cnt
    iown=get_det_owner(walk_dets_up(isnd_start+isnd),walk_dets_dn(isnd_start+isnd))
    snd_list(isnd)%to=iown
    snd_list(isnd)%wlkid=isnd_start+isnd
    snd_table(whoami*ncores + iown + 1) = snd_table(whoami*ncores + iown + 1) + 1
  enddo
    call mpi_allred(snd_table)

  recv_displs=0
  do icore=1,(ncores-1)
    recv_displs(icore+1)=snd_table((icore-1)*ncores+whoami+1)
  enddo
  do icore=2,ncores
    recv_displs(icore)=recv_displs(icore)+recv_displs(icore-1)
  enddo
  recv_displs(1)=0

  recv_displs=recv_displs+my_nwalk_1-1

  tot_snd=sum(snd_table(whoami*ncores + 1:(whoami+1)*ncores))

  snd_displs=0
  do icore=2,ncores
    snd_displs(whoami*ncores + icore)=sum(snd_table(whoami*ncores + 1:whoami*ncores + icore-1))
  enddo

  if(tot_snd > 0) then
    do icore=0,(ncores-1)
        iown=0
        my_ic=whoami*ncores + icore + 1
        do isnd=1,(tot_snd+my_own)
          if(snd_list(isnd)%to == icore) then
            iown=iown+1
            list_cowner(iown)=snd_list(isnd)%wlkid
          endif
          if(iown.ge.snd_table(my_ic)) exit
        enddo
        if(iown /= snd_table(my_ic)) stop "something's WRONG!!!"

        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%wt=walk_wt(list_cowner(1:snd_table(my_ic)))
        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%imp_ini(1) = imp_distance(list_cowner(1:snd_table(my_ic)))
        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%imp_ini(2) = initiator(list_cowner(1:snd_table(my_ic)))
        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%e_num = e_num(list_cowner(1:snd_table(my_ic)))
        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%e_den = e_den(list_cowner(1:snd_table(my_ic)))
        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%diag_elems = diag_elems(list_cowner(1:snd_table(my_ic)))
        do iw=1,snd_table(my_ic)
#ifdef NUM_ORBITALS_GT_127
            do islice=1,num_words
                call conv_128_to_64(walk_dets_up(list_cowner(iw))%v(islice),snd_buff(snd_displs(my_ic)+iw)%det_u(2*(islice-1)+1),snd_buff(snd_displs(my_ic)+iw)%det_u(2*(islice-1)+2))
                call conv_128_to_64(walk_dets_dn(list_cowner(iw))%v(islice),snd_buff(snd_displs(my_ic)+iw)%det_d(2*(islice-1)+1),snd_buff(snd_displs(my_ic)+iw)%det_d(2*(islice-1)+2))
            enddo
#else
            call conv_128_to_64(walk_dets_up(list_cowner(iw)),snd_buff(snd_displs(my_ic)+iw)%det_u(1),snd_buff(snd_displs(my_ic)+iw)%det_u(2))
            call conv_128_to_64(walk_dets_dn(list_cowner(iw)),snd_buff(snd_displs(my_ic)+iw)%det_d(1),snd_buff(snd_displs(my_ic)+iw)%det_d(2))
#endif
        enddo
        snd_matrix_elements(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic)) = matrix_elements(list_cowner(1:snd_table(my_ic)))
    enddo
  endif

  do icore=0,(ncores-1)
    if(sum(snd_table(icore*ncores + 1:(icore+1)*ncores)) > 0) then !**do the following only if icore has something for the others
      call mpi_scattv(snd_buff(1:tot_snd),snd_table(icore*ncores + 1:(icore+1)*ncores),snd_displs(icore*ncores + 1:(icore+1)*ncores),recv_buff(recv_displs(icore+1)+1:recv_displs(icore+1)+snd_table(icore*ncores+ whoami +1)),snd_table(icore*ncores+ whoami +1),icore)
      call mpi_scattv(snd_matrix_elements(1:tot_snd),snd_table(icore*ncores + 1:(icore+1)*ncores),snd_displs(icore*ncores + 1:(icore+1)*ncores),matrix_elements(recv_displs(icore+1)+1:recv_displs(icore+1)+snd_table(icore*ncores+ whoami +1)),snd_table(icore*ncores+ whoami +1),icore)
    endif
  enddo

  do icore=0,(ncores-1)
    my_nwalk = my_nwalk +snd_table(icore*ncores+whoami+1)
  enddo

  walk_wt(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%wt
  do iw=recv_displs(1)+1,my_nwalk
#ifdef NUM_ORBITALS_GT_127
    do islice=1,num_words
      call conv_64_to_128(walk_dets_up(iw)%v(islice),recv_buff(iw)%det_u(2*(islice-1)+1),recv_buff(iw)%det_u(2*(islice-1)+2))
      call conv_64_to_128(walk_dets_dn(iw)%v(islice),recv_buff(iw)%det_d(2*(islice-1)+1),recv_buff(iw)%det_d(2*(islice-1)+2))
    enddo
#else
    call conv_64_to_128(walk_dets_up(iw),recv_buff(iw)%det_u(1),recv_buff(iw)%det_u(2))
    call conv_64_to_128(walk_dets_dn(iw),recv_buff(iw)%det_d(1),recv_buff(iw)%det_d(2))
#endif
  enddo
  imp_distance(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%imp_ini(1)
  initiator(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%imp_ini(2)
  e_num(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%e_num
  e_den(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%e_den
  diag_elems(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%diag_elems

  snd_table=0
#endif
end subroutine mpi_sendwalks
!------------------------------------------------------------------------------

subroutine mpi_push_nwalk(iw,iinitiator,e_num,e_den,diag_elem)
integer,intent(in) :: iw
integer(i1b), intent(in) :: iinitiator
real(rk),intent(in) :: e_num,e_den,diag_elem

integer :: iown,sdispls
#ifdef NUM_ORBITALS_GT_127
integer :: islice
#endif

#ifdef MPI
iown=get_det_owner(walk_dets_up(iw),walk_dets_dn(iw))
snd_cnt(iown+1) = snd_cnt(iown+1) + 1

if(snd_cnt(iown+1) > nwbuff) stop "Need to enlarge dimension of SND_BUFFs!"

sdispls=nwbuff*iown+snd_cnt(iown+1)

snd_buff(sdispls)%wt=walk_wt(iw)
snd_buff(sdispls)%imp_ini(1)=imp_distance(iw)
snd_buff(sdispls)%imp_ini(2)=iinitiator
snd_buff(sdispls)%e_num=e_num
snd_buff(sdispls)%e_den=e_den
snd_buff(sdispls)%diag_elems=diag_elem

#ifdef NUM_ORBITALS_GT_127
do islice=1,num_words
    call conv_128_to_64(walk_dets_up(iw)%v(islice),snd_buff(sdispls)%det_u(2*(islice-1)+1),snd_buff(sdispls)%det_u(2*(islice-1)+2))
    call conv_128_to_64(walk_dets_dn(iw)%v(islice),snd_buff(sdispls)%det_d(2*(islice-1)+1),snd_buff(sdispls)%det_d(2*(islice-1)+2))
enddo
#else
call conv_128_to_64(walk_dets_up(iw),snd_buff(sdispls)%det_u(1),snd_buff(sdispls)%det_u(2))
call conv_128_to_64(walk_dets_dn(iw),snd_buff(sdispls)%det_d(1),snd_buff(sdispls)%det_d(2))
#endif

#endif
end subroutine mpi_push_nwalk
!------------------------------------------------------------------------------

!**improved version to use when we know that the walkers are newly generated:
!** - don't sends matrix_elements but set them to 1e51_rk
!** - makes just one scatterv with a derived type
subroutine mpi_sendnewwalks(my_nwalk,initiator,e_num,e_den,diag_elems)
    integer(i1b), intent(inout) :: initiator(:)
    integer,intent(inout) :: my_nwalk
    real(rk),intent(inout) :: e_num(:),e_den(:),diag_elems(:)

    integer :: tot_snd,icore,my_own,iw
#ifdef NUM_ORBITALS_GT_127
    integer :: islice
#endif
#ifdef MPI
    call mpi_agath_int(snd_cnt,ncores,snd_table,ncores)

    recv_displs=0
    do icore=1,(ncores-1)
        recv_displs(icore+1)=snd_table((icore-1)*ncores+whoami+1)
    enddo
    do icore=2,ncores
        recv_displs(icore)=recv_displs(icore)+recv_displs(icore-1)
    enddo
    recv_displs(1)=0

    recv_displs=recv_displs+my_nwalk

    my_own=snd_table(whoami*ncores+whoami+1)
    snd_table(whoami*ncores+whoami+1) = 0 !**don't send back my own walkers

    tot_snd=sum(snd_table(whoami*ncores + 1:(whoami+1)*ncores))

    do icore=0,(ncores-1)
      rcv_cnt(1+icore)=snd_table(icore*ncores+ whoami +1)
    enddo

    call mpi_alltoallv_iwalk(snd_buff(1:tot_snd),snd_table(whoami*ncores + 1:(whoami+1)*ncores),snd_dsp,recv_buff,rcv_cnt,recv_displs)

    my_nwalk = my_nwalk + sum(rcv_cnt(1:ncores))

    my_nwalk = my_nwalk+my_own

    if(.not.master_core) then
      walk_wt(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%wt
      do iw=recv_displs(1)+1,recv_displs(whoami+1)
#ifdef NUM_ORBITALS_GT_127
        do islice=1,num_words
          call conv_64_to_128(walk_dets_up(iw)%v(islice),recv_buff(iw)%det_u(2*(islice-1)+1),recv_buff(iw)%det_u(2*(islice-1)+2))
          call conv_64_to_128(walk_dets_dn(iw)%v(islice),recv_buff(iw)%det_d(2*(islice-1)+1),recv_buff(iw)%det_d(2*(islice-1)+2))
        enddo
#else
        call conv_64_to_128(walk_dets_up(iw),recv_buff(iw)%det_u(1),recv_buff(iw)%det_u(2))
        call conv_64_to_128(walk_dets_dn(iw),recv_buff(iw)%det_d(1),recv_buff(iw)%det_d(2))
#endif
      enddo
      imp_distance(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%imp_ini(1)
      initiator(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%imp_ini(2)
      e_num(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%e_num
      e_den(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%e_den
      diag_elems(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%diag_elems
    endif

    walk_wt(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%wt
    do iw=1,my_own
#ifdef NUM_ORBITALS_GT_127
        do islice=1,num_words
          call conv_64_to_128(walk_dets_up(recv_displs(whoami+1)+iw)%v(islice),snd_buff(snd_dsp(whoami+1)+iw)%det_u(2*(islice-1)+1),snd_buff(snd_dsp(whoami+1)+iw)%det_u(2*(islice-1)+2))
          call conv_64_to_128(walk_dets_dn(recv_displs(whoami+1)+iw)%v(islice),snd_buff(snd_dsp(whoami+1)+iw)%det_d(2*(islice-1)+1),snd_buff(snd_dsp(whoami+1)+iw)%det_d(2*(islice-1)+2))
        enddo
#else
        call conv_64_to_128(walk_dets_up(recv_displs(whoami+1)+iw),snd_buff(snd_dsp(whoami+1)+iw)%det_u(1),snd_buff(snd_dsp(whoami+1)+iw)%det_u(2))
        call conv_64_to_128(walk_dets_dn(recv_displs(whoami+1)+iw),snd_buff(snd_dsp(whoami+1)+iw)%det_d(1),snd_buff(snd_dsp(whoami+1)+iw)%det_d(2))
#endif
    enddo
    imp_distance(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%imp_ini(1)
    initiator(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%imp_ini(2)
    e_num(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%e_num
    e_den(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%e_den
    diag_elems(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%diag_elems

    if(whoami+2 .le. ncores) then
      walk_wt(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%wt
      do iw=recv_displs(whoami+2)+1,my_nwalk
#ifdef NUM_ORBITALS_GT_127
        do islice=1,num_words
          call conv_64_to_128(walk_dets_up(iw)%v(islice),recv_buff(iw)%det_u(2*(islice-1)+1),recv_buff(iw)%det_u(2*(islice-1)+2))
          call conv_64_to_128(walk_dets_dn(iw)%v(islice),recv_buff(iw)%det_d(2*(islice-1)+1),recv_buff(iw)%det_d(2*(islice-1)+2))
        enddo
#else
        call conv_64_to_128(walk_dets_up(iw),recv_buff(iw)%det_u(1),recv_buff(iw)%det_u(2))
        call conv_64_to_128(walk_dets_dn(iw),recv_buff(iw)%det_d(1),recv_buff(iw)%det_d(2))
#endif
      enddo
      imp_distance(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%imp_ini(1)
      initiator(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%imp_ini(2)
      e_num(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%e_num
      e_den(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%e_den
      diag_elems(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%diag_elems
    endif

    matrix_elements(recv_displs(1)+1:my_nwalk) = 1e51_rk
#else
    my_nwalk = my_nwalk + sum(snd_cnt)
#endif
end subroutine mpi_sendnewwalks
!---------------------------------------------------------------

!------------Shared Memory Allocation---------------------------
!-- This routine uses MPI3 routines the allocate shared memory--
!-- MPI routines used:
!-- MPI_WIN_ALLOCATE_SHARED, MPI_WIN_SHARED_QUERY
!-- Additionally, C_F_POINTER (from ISO_C_BINDING) is used
!-- to set up the pointers needed.
!--
!--Use: call shmem_allocate(array,numberOfElements)
!-- in: numberOfElements: the number of elements of the array
!--                       to be allocated
!-- out: array:           pointer to the sharedMemory array
!--                       allocated
!--------------------------------------------------------------
!--Added by: Matthew Otten, April 26, 2014
!---------------------------------------------------------------
!---------------------------------------------------------------
!-------------How to Implement----------------------------------
!-- for any array for which you want shared, declare it as
!-- real(rk),dimension(:),pointer :: array
!-- and then call shmem_allocate(array,numberOfElements)
!--
!-- To change from arrays that are currently allocated:
!-- change real(rk),allocatable :: array(:) to
!-- real(rk),dimension(:),pointer :: array
!-- and change allocate(array(numberOfElements))
!-- to call shmem_allocate(array,numberOfElements)
!---------------------------------------------------------------

subroutine shmem_allocate_rk(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                        my_ncores_on_node
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_local_cptr,array_cptr
  !This array is a f90 pointer
  real(rk),dimension(:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements
  integer :: info

 !if MPI is defined, use the MPI sharedmenmory things
#ifdef MPI
  !figure out how many nodes on a core
  call MPI_COMM_SIZE(sharedComm,my_ncores_on_node,mpierr)
  !set the buffer size; it should be the total size/ncores_on_node
  !1._rk is just chosen as an arbitrary number of the type we
  !want the sizeof
  write (6,*) "number of elements, ncores on node=",numberOfElements,my_ncores_on_node; call flush(6)
  shmem_buf_size = ceiling(real(numberOfElements)/real(my_ncores_on_node))*sizeof(1._rk)
  write (6,*) "Allocating shared memory of size",shmem_buf_size; call flush(6)
  !Allocate shared memory. Arguments: bufferSize, elementSize, MPI_INFO
  !communicator,c_ptr for array,window,mpierr
  win_num = win_num + 1
  call MPI_INFO_CREATE(info,mpierr)
  call MPI_INFO_SET(info,'alloc_shared_noncontig','true',mpierr)
  call MPI_WIN_ALLOCATE_SHARED(shmem_buf_size,sizeof(1._rk),info,sharedComm,array_local_cptr,shmemWindow(win_num),mpierr)
  !query the master processor (of the window [so, of the sharedComm, or of the node]) for
  !the location of the start of the memory.
  call MPI_WIN_SHARED_QUERY(shmemWindow(win_num),0,shmem_buf_size,shmem_buf_size,array_cptr,mpierr)
  !remap that c_ptr to the array
  call C_F_POINTER(array_cptr,array,[numberOfElements])
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements))
#endif
end subroutine shmem_allocate_rk
!------------------------------------------------------------------------------

subroutine shmem_allocate_real(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  ! Added type real to Matt's shmem_allocate routine
  ! A Holmes, 5 Apr 2015
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                        my_ncores_on_node
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_local_cptr,array_cptr
  !This array is a f90 pointer
  real,dimension(:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements
  real :: test
  integer :: info

 !if MPI is defined, use the MPI sharedmenmory things
#ifdef MPI
  !figure out how many nodes on a core
  call MPI_COMM_SIZE(sharedComm,my_ncores_on_node,mpierr)
  !set the buffer size; it should be the total size/ncores_on_node
  !1._rk is just chosen as an arbitrary number of the type we
  !want the sizeof
  write (6,*) "number of elements, ncores on node=",numberOfElements,my_ncores_on_node; call flush(6)
  shmem_buf_size = ceiling(real(numberOfElements)/real(my_ncores_on_node))*sizeof(test)
  write (6,*) "Allocating shared memory of size",shmem_buf_size; call flush(6)
  !Allocate shared memory. Arguments: bufferSize, elementSize, MPI_INFO
  !communicator,c_ptr for array,window,mpierr
  win_num = win_num + 1
  call MPI_INFO_CREATE(info,mpierr)
  call MPI_INFO_SET(info,'alloc_shared_noncontig','true',mpierr)
  call MPI_WIN_ALLOCATE_SHARED(shmem_buf_size,sizeof(test),info,sharedComm,array_local_cptr,shmemWindow(win_num),mpierr)
 !call MPI_WIN_ALLOCATE_SHARED(shmem_buf_size,sizeof(test),MPI_INFO_NULL,sharedComm,array_local_cptr,shmemWindow(win_num),mpierr)
  !query the master processor (of the window [so, of the sharedComm, or of the node]) for 
  !the location of the start of the memory.
  call MPI_WIN_SHARED_QUERY(shmemWindow(win_num),0,shmem_buf_size,shmem_buf_size,array_cptr,mpierr)
  !remap that c_ptr to the array
  call C_F_POINTER(array_cptr,array,[numberOfElements])
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements))
#endif 
end subroutine shmem_allocate_real
!------------------------------------------------------------------------------

subroutine shmem_allocate_int(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  ! Added type int to Matt's shmem_allocate routine
  ! A Holmes, 5 Apr 2015
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                        my_ncores_on_node
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_local_cptr,array_cptr
  !This array is a f90 pointer
  integer,dimension(:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements
  integer :: test
  integer :: info

 !if MPI is defined, use the MPI sharedmenmory things
#ifdef MPI
  !figure out how many nodes on a core
  call MPI_COMM_SIZE(sharedComm,my_ncores_on_node,mpierr)
  !set the buffer size; it should be the total size/ncores_on_node
  !1._rk is just chosen as an arbitrary number of the type we
  !want the sizeof
  write (6,*) "number of elements, ncores on node=",numberOfElements,my_ncores_on_node; call flush(6)
  shmem_buf_size = ceiling(real(numberOfElements)/real(my_ncores_on_node))*sizeof(test)
  write (6,*) "Allocating shared memory of size",shmem_buf_size; call flush(6)
  !Allocate shared memory. Arguments: bufferSize, elementSize, MPI_INFO
  !communicator,c_ptr for array,window,mpierr
  win_num = win_num + 1
  call MPI_INFO_CREATE(info,mpierr)
  call MPI_INFO_SET(info,'alloc_shared_noncontig','true',mpierr)
  call MPI_WIN_ALLOCATE_SHARED(shmem_buf_size,sizeof(test),info,sharedComm,array_local_cptr,shmemWindow(win_num),mpierr)
  !query the master processor (of the window [so, of the sharedComm, or of the node]) for 
  !the location of the start of the memory.
  call MPI_WIN_SHARED_QUERY(shmemWindow(win_num),0,shmem_buf_size,shmem_buf_size,array_cptr,mpierr)
  !remap that c_ptr to the array
  call C_F_POINTER(array_cptr,array,[numberOfElements])
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements))
#endif 
end subroutine shmem_allocate_int
!------------------------------------------------------------------------------

subroutine mpi_stop(msg)
  USE ISO_FORTRAN_ENV, only : error_unit ! access computing environment
  character(len=*) msg
  integer mpi_err

  write(6,'(a)') msg ; call flush(6)
  write(error_unit,'(a)') msg ! write message to standard error
#ifdef MPI
  call MPI_ABORT(MPI_COMM_WORLD,1,mpi_err)
#else
  stop 'mpi_stop'
#endif

end subroutine mpi_stop

end module mpi_routines
