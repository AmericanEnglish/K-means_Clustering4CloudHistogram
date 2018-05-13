module kmeans_omp_module
  !$ use omp_lib
  implicit none
  include 'mpif.h'

contains
  !!! Naiive MPI Wrappers
  !!! Comm = MPI_COMM_WORLD to prevent any bad datatype transformations.
  function mpiInit() result(ierr)
    integer :: ierr
    !F2PY INTENT(OUT) :: ierr
    call MPI_INIT(ierr)
    return 
  end function mpiInit

  !!! Get rank
  function mpiGetRank() result(rank) 
    integer :: ierr, rank
    !F2PY INTENT(OUT) :: rank
    call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
    return
  end function mpiGetRank

  !!! Get size
  function mpiGetSize() result(tprocs)
    integer :: ierr, tprocs
    !F2PY INTENT(OUT) :: tprocs
    call MPI_COMM_SIZE(MPI_COMM_WORLD, tprocs, ierr)
    return
  end function mpiGetSize

  function set_num_threads(nn) result(nt)
    integer :: nn,nt
    !F2PY INTENT(IN) :: nn
    !F2PY INTENT(OUT) :: nt
    
    !$ nt = omp_get_max_threads()
    print*,"Max # of Threads=",nt
    !$ call omp_set_num_threads(nn)
    !$ nt = omp_get_max_threads()
    print*, "Set Max # of Threads=",nt
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

  subroutine
      assign_and_get_newsum(indata,ctd,nk,cl,outsum,ncl,nelem,nrec,rank,tprocs)
  !!! Calculate sum of data by clusters
  !!! Need to get new centroid
  !!! ncl = number of clusters
  !!! nk = limited number of records to use
  !!! ii,jj,kk = indexes used during iteration
  !!! nelem =  number of bins in the histogram
  !!! nrec = the number of total records in the dataset
  !!! l_nrec = the number of records that each process should handle
  !!! rem = remainder
  !!! startRec = each process' record to start mathing at
  !!! stopRec = each process' last record to process
  !!! tcl = Temporary CL which MPI can use for reduce operations
    integer :: nelem,nrec,nk,ncl,rank,tprocs
    integer :: ii,jj,kk,idx,l_nrec,rem,startRec,stopRec,ierr
    integer, intent(out) :: cl(nrec), tcl(nrec)
    real(8), intent(in) :: indata(nelem,nrec),ctd(nelem,ncl)
    real(8), intent(out) :: outsum(nelem,ncl)
    real(8) :: mindd,tmpdd
    !F2PY INTENT(HIDE) :: ncl,nelem,nrec
    !F2PY INTENT(OUT) :: cl,outsum
    !F2PY INTENT(IN) :: nk
    
    outsum=0.d0
    !!! 0 out cl, so MPI doesn't cause artifacts
    cl = 0
    tcl = 0
    !!! Calculate the total number of records for each process
    l_nrec = nrec / tprocs
    rem = nrec % tprocs
    if (rem == 0)
      startRec = l_nrec*rank 
    else 
      if (rank < rem)
        !!! Pick up an extra record
        startRec = l_nrec*rank + 1
      else
        !!! Accounts for additional records
        startRec = l_nrec*rank + rem
    
    !!! Fortran indexes at 1
    startRec = startRec + 1
    stopRec = startRec + l_nrec


    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(indata,ctd,cl,tcl,outsum,ncl,nk,nrec,nelem,startRec,stopRec)

    !!!--- Assigning Clusters

    !!! From 1->nrec as 1+n*nk

    !$OMP DO 
    !!!do ii=1,nrec,nk
    do ii=startRec,stopRec,nk
       !!! Min distance
       mindd=10000.d0;idx=0
       do kk=1,ncl
          tmpdd=calc_sqdist(indata(:,ii), ctd(:,kk), nelem)
          !!! If calculated distance is less than min distance
          if (tmpdd.lt.mindd) then
             mindd=tmpdd; idx=kk
          endif
       enddo
       if (idx.eq.0) then
          print*,"Not assigned",idx,ii,indata(:,ii)
          stop
       endif
       tcl(ii)=idx
    enddo

    
    !$OMP END DO
    !$OMP BARRIER

    
    !!! Do this with only 1 thread for face the consequences of the ultimate
    !!! race condition
    !$OMP SINGLE
    !!! Merge tcl's across all processes into the complete cl
    ! send data, recieve data, entries, datatype, operation, communicator, ierr
    MPI_ALLREDUCE(tcl, cl, nrec, MPI_DOUBLE_PRECISION, MPI_SUM, 
    &             MPI_COMM_WORLD, ierr)
    !$OMP END SINGLE
    !!! There is an implied barrier using the END clause 


    !!!--- Sum for New Centroid

    !!! This loop order keeps OpenMP thread safe but will require an MPI reduce
    !!! to merge all of the histograms back into the correct data structure
    !!! if the loop is parallelized via MPI as well

    
    !$OMP DO
    do ii=1,nelem
       do jj=1,nrec,nk
          outsum(ii,cl(jj))=outsum(ii,cl(jj))+indata(ii,jj)
       enddo
    enddo
    !$OMP END DO
    !print*, outsum,cluster(1000:1100)
    !$OMP END PARALLEL

  end subroutine assign_and_get_newsum

  subroutine get_wcv_sum(indata,ctd,cl,outsum,ncl,nelem,nrec)
  !!! Calculate sum of data by clusters
  !!! Need to get new centroid
    integer :: nelem,nrec,ncl
    integer :: mm,ii
    integer, intent(in) :: cl(nrec)
    real(8), intent(in) :: indata(nelem,nrec),ctd(nelem,ncl)
    real(8), intent(out) :: outsum(nelem,ncl)
    !F2PY INTENT(HIDE) :: ncl,nelem,nrec
    !F2PY INTENT(OUT) :: outsum

    outsum=0.d0
    !!! Similar story to the previous outsum loop in assign_and_get
    !!! can be parallelized but will require a Reduce call in order
    !!! have the correct numbers inside outsum
    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(outsum,cl,indata,ctd,nrec,nelem)
    do mm=1,nelem
       do ii=1,nrec
          outsum(mm,cl(ii))=outsum(mm,cl(ii))+(indata(mm,ii)-ctd(mm,cl(ii)))**2
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine get_wcv_sum


end module kmeans_omp_module
