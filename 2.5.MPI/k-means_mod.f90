module kmeans_omp_module
  !$ use omp_lib
  implicit none
  include "mpif.h"

contains
  
  function mpiInit() result(ierr)
    integer :: ierr
    !F2PY INTENT(OUT) :: ierr
    call MPI_INIT(ierr)
    return 
  end function mpiInit
  ! Need MPI wall time
  function mpiwtime() result(seconds)
      DOUBLE PRECISION :: seconds
      !F2PY INTENT(OUT) :: seconds
      seconds = MPI_WTIME()
      return
  end function
  ! Get rank
  function mpiGetRank() result(rank) 
    integer :: ierr, rank
    !F2PY INTENT(OUT) :: rank
    call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
    return
  end function mpiGetRank

  ! Get size
  function mpiGetSize() result(tprocs)
    integer :: ierr, tprocs
    !F2PY INTENT(OUT) :: tprocs
    call MPI_COMM_SIZE(MPI_COMM_WORLD, tprocs, ierr)
    return
  end function mpiGetSize
  ! Comm = MPI_COMM_WORLD to prevent any bad datatype transformations.

  function set_num_threads(nn) result(nt)
    integer :: nn,nt
    integer :: rank
    !F2PY INTENT(IN) :: nn
    !F2PY INTENT(OUT) :: nt
    
    !$ nt = omp_get_max_threads()
    rank = mpiGetRank()
    if (rank == 0) then
        print*,"Max # of Threads=",nt
        !$ call omp_set_num_threads(nn)
        !!! nt = omp_get_max_threads()
        print*, "Set Max # of Threads=",nt
    endif
  end function set_num_threads

  real(8) function calc_sqdist(aa,bb,nn)
  !!! Calculate squared sum of distance
    integer :: nn,ii
    real(8), dimension(nn),intent(in) :: aa,bb
    !F2PY INTENT(HIDE) :: nn
    !F2PY INTENT(OUT) :: calc_sqdist

    calc_sqdist=0.d0
    do ii=1,nn
       calc_sqdist=calc_sqdist+(aa(ii)-bb(ii))**2
    enddo

    return
  end function calc_sqdist

  real(8) function calc_dist(aa,bb,nn) 
    integer :: nn,ii
    real(8),intent(in) :: aa(nn),bb(nn) 
!    real(8) :: calc_sqdist
    !F2PY INTENT(HIDE) :: nn
    !F2PY INTENT(OUT) :: calc_dist

    calc_dist=sqrt(calc_sqdist(aa,bb,nn))

    return
  end function calc_dist

  subroutine &
      assign_and_get_newsum(indata,ctd,nk,cl,outsum,ncl,nelem,nrec,rank,tprocs)
  !!! Calculate sum of data by clusters
  !!! Need to get new centroid
    integer :: nelem,nrec,nk,ncl
    integer :: ii,jj,kk,idx
    integer, intent(out) :: cl(nrec)
    real(8), intent(in) :: indata(nelem,nrec),ctd(nelem,ncl)
    real(8), intent(out) :: outsum(nelem,ncl)
    real(8) :: mindd,tmpdd
    !!! MPI VARIABLES
    integer, intent(in) :: rank, tprocs
    integer :: startRec, stopRec, ierr
    integer :: tcl(nrec)
    real(8) :: toutsum(nelem, ncl)
    !F2PY INTENT(HIDE) :: ncl,nelem,nrec
    !F2PY INTENT(OUT) :: cl,outsum
    !F2PY INTENT(IN) :: nk

    cl = 0
    tcl = 0

    !rank = mpigetrank()
    !tprocs = mpigetsize()
    call get_record_spans(nrec, rank, tprocs, startRec, stopRec)
    !print*,"rank",rank,"tprocs",tprocs,"startRec",startRec,"stopRec",stopRec,"nrec",nrec

    outsum=0.d0
    toutsum=0.d0

    !!!--- Assigning Clusters

    ! If OpenMP can split this, so can MPI!

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(indata,ctd,cl,tcl,outsum,ncl,nk,nrec,nelem,startRec,stopRec)
    do ii=startRec,stopRec,nk
    !do ii=1,nrec,nk
       mindd=10000.d0;idx=0
       do kk=1,ncl
          tmpdd=calc_sqdist(indata(:,ii),ctd(:,kk),nelem)
          if (tmpdd.lt.mindd) then
             mindd=tmpdd; idx=kk
          endif
       enddo
       if (idx.eq.0) then
          print*,"Not assigned",idx,ii,indata(:,ii)
          stop
       endif
       !cl(ii)=idx
       tcl(ii)=idx
    enddo
    !$OMP END PARALLEL DO
    !!!--- Sum for New Centroid
    !cl = tcl
    call MPI_ALLREDUCE(tcl, cl, nrec, MPI_INTEGER, &
      MPI_SUM, MPI_COMM_WORLD, ierr)
    ! Loop order needs to be swapped for MPI
    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(indata,ctd,cl,outsum,toutsum,ncl,nk,nrec,nelem, startRec, stopRec)
    do ii=1,nelem
       do jj=startRec,stopRec,nk
       !do jj=1,nrec,nk
          toutsum(ii,cl(jj))=toutsum(ii,cl(jj))+indata(ii,jj)
       enddo
    enddo
    !$OMP END PARALLEL DO

    !outsum = toutsum
    do ii=1,ncl
      call MPI_ALLREDUCE(toutsum(:,ii), outsum(:,ii), nelem, MPI_REAL8, &
        MPI_SUM, MPI_COMM_WORLD, ierr)
    enddo
    !print*, outsum,cluster(1000:1100)

  end subroutine assign_and_get_newsum

  subroutine get_wcv_sum(indata,ctd,cl,outsum,ncl,nelem,nrec,rank,tprocs)
  !!! Calculate sum of data by clusters
  !!! Need to get new centroid
    integer :: nelem,nrec,ncl
    integer :: mm,ii
    integer, intent(in) :: cl(nrec)
    real(8), intent(in) :: indata(nelem,nrec),ctd(nelem,ncl)
    real(8), intent(out) :: outsum(nelem,ncl)
    !!! MPI VARIABLES
    real(8) :: toutsum(nelem, ncl)
    integer :: startRec, stopRec, ierr
    integer, intent(in) :: rank, tprocs
    !F2PY INTENT(HIDE) :: ncl,nelem,nrec
    !F2PY INTENT(OUT) :: outsum

    outsum=0.d0
    toutsum=0.d0
    !rank = mpigetrank()
    !tprocs = mpigetsize()
    call get_record_spans(nrec, rank, tprocs, startRec, stopRec)
    ! Needs to be adapted for MPI
    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(toutsum,cl,indata,ctd,nrec,nelem,startRec,stopRec)
    do mm=1,nelem
        do ii=startRec,stopRec
          toutsum(mm,cl(ii))=toutsum(mm,cl(ii))+(indata(mm,ii)-ctd(mm,cl(ii)))**2
       enddo
    enddo
    !$OMP END PARALLEL DO
    do ii=1,ncl
      call MPI_ALLREDUCE(toutsum(:,ii), outsum(:,ii), nelem, MPI_REAL8, &
        MPI_SUM, MPI_COMM_WORLD, ierr)
    enddo

    ! An MPI sum would have to happen here
  end subroutine get_wcv_sum

  subroutine get_record_spans(nrec,rank,tprocs,startRec,stopRec)
    integer :: l_nrec, rem
    integer, intent(out) :: startRec, stopRec
    integer, intent(in) :: nrec, rank, tprocs
    !F2PY INTENT(OUT) :: startRec,stopRec
    !F2PY INTENT(IN)  :: nrec,rank,tprocs
    !!! Calculate the total number of records for each process
    l_nrec = nrec / tprocs
    rem = MOD(nrec,  tprocs)
    if (rem == 0) then
      startRec = l_nrec*rank 
    else 
      if (rank < rem) then
        !!! Pick up an extra record
        l_nrec = l_nrec + 1
        startRec = rank*l_nrec
      else
        !!! Accounts for additional records
        startRec = l_nrec*rank + rem
      endif
    endif
    !!! Fortran indexes at 1
    startRec = startRec + 1
    stopRec = startRec + l_nrec - 1
  end subroutine

end module kmeans_omp_module
