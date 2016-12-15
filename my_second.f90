  subroutine my_second (n,title)
! Written by Cyrus Umrigar
! Prints out cpu and wall-clock time,
! both since beginning of run and since last call.

! I used to use etime for cpu time and time for wall-clock time.
! I now use cpu_time and system_clock because they are standard.
! I tested out time and cpu_time for cpu time and
! system_clock and secnds for wall-clock time (system_clock is standard)
! to see if they are more accurate.  All four gave variations of
! a few tenths of a second for a job that ran less than 10 secs.
! cpu_time and time were the same using g77 in Linux

!**** Edited by AR[7/23/13]: for now the job is done just by master_core, may be extended to act globally later 

  use types, only: rk, i8b
  
  use mpi_routines, only: master_core
  
  implicit none
  integer, intent(in) :: n
  character(len=*), intent(in) :: title

  real(rk) etim,etim1,etim2,etimtot,etimlast
! integer itim,itim1, itim2, itimtot, itimlast, icall, icount_rate, icountmax
  integer(i8b) itim,itim1, itim2, itimtot, itimlast, icall, icount_rate, icountmax

  save icall,itim1,itim2,etim1,etim2
  data icall/0/

! etim=etime(i)        ! unix
! etim=etime(tarray)   ! linux
! etim=mclock()*1.d-2  ! aix
! etim=mclock()*1.d-6

  call cpu_time(etim)  ! standard but use pgf90, not pgf77

! itim=time()
! itim=time(0)       ! aix

  call system_clock(itim,icount_rate,icountmax) ! standard but but use pgf90, not pgf77

  if(icall.eq.0) then
    icall=1
    itim1=itim
    itim2=itim
    etim1=etim
    etim2=etim
  endif

  if(itim.gt.itim1) then
    itimtot=nint(real(itim-itim1,rk)/icount_rate)
  else
    itimtot=nint(real(itim-itim1+icountmax,rk)/icount_rate)
  endif
  if(itim.gt.itim2) then
    itimlast=nint(real(itim-itim2,rk)/icount_rate)
  else
    itimlast=nint(real(itim-itim2+icountmax,rk)/icount_rate)
  endif
  itim2=itim

  etimtot=etim-etim1
  etimlast=etim-etim2
  etim2=etim

  if(n.eq.1) write (6,'(/,''BEGINNING OF '',a,'' CP, REAL TIME IS '', 2f11.2,2i9)') title,etimtot,etimlast,itimtot,itimlast
  if(n.eq.2) write (6,'(''END       OF '',a,'' CP, REAL TIME IS '', 2f11.2,2i9)') title,etimtot,etimlast,itimtot,itimlast

#ifdef DEBUG
  call my_memory() ! Also print out memory usage when in debug mode.
#endif

  call flush(6)

  return
  end
