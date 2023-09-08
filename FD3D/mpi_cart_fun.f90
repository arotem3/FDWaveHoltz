!c-----------------------------------------------------------------------

subroutine PROCDIST3D( nig, njg, nkg, nproc, pmax )

  !***********************************************************************
  !***
  !*** Distribute processors in three space dimensions, such that each
  !*** patch becomes as close to a square as possible.
  !***
  !*** Input: nig, njg - Global array sizes.
  !***        nprocs   - Total number of processors.
  !***
  !*** Output: pmax  - Number of processors in each direction.
  !***
  !***********************************************************************

  implicit none
  integer nig, njg, nkg, nproc, pmax(3), p1, p2, p3
  real*8  f, fmin
  logical first
  data first/.true./

  fmin = nig+njg+nkg

  do p3=1,nproc
     do p1=1,nproc
        if( MOD(nproc,p1*p3).eq.0 )then
           p2 = nproc/(p1*p3)
           f = ABS( DBLE(nig)/p1 - DBLE(njg)/p2 ) +&
                ABS( DBLE(njg)/p2 - DBLE(nkg)/p3 )
           if( f .lt. fmin .or. first )then
              fmin = f
              pmax(1) = p1
              pmax(2) = p2
              pmax(3) = p3
              first = .false.
           endif
        endif
     enddo
  enddo
end subroutine PROCDIST3D


!c-----------------------------------------------------------------------
subroutine BLOCKDEF( n, b, igtr, nloc, olap )

  !***********************************************************************
  !***
  !*** Distribute blocks as evenly as possible, using overlap olap.
  !***
  !***  Input: n - Number of points
  !***         b - Number of blocks to split the points into
  !***         olap - Size of overlap between blocks ( in number of points )
  !***
  !***  Output: igtr(*) - Transform such that iglobal = ilocal + igtr(k)
  !***             iglobal - Global index, ilocal - starting from 1 in block k
  !***          nloc(*) - Number of points in each block
  !***
  !***********************************************************************

  implicit none

  integer n, b, igtr(b), nloc(b), olap
  integer s0, r, k

  s0 = ( n + (b-1)*olap )/b
  r  = MOD( n + (b-1)*olap, b )

  do k=1,r
     nloc(k) = s0+1
  enddo
  do k=r+1,b
     nloc(k) = s0
  enddo
  do k=1,b
     igtr(k) = (k-1)*(s0-olap) + MIN( k-1, r )
  enddo

  return
end subroutine BLOCKDEF

!c-----------------------------------------------------------------------
subroutine MYMESH( myname, pmax, proc )

  !***********************************************************************
  !***
  !*** My ID in three dimensional processor mesh.
  !***
  !*** Mapping:
  !***       myname = proc(1)-1 + pmax(1)*(proc(2)-1)+
  !***                                      pmax(2)*pmax(1)*(proc(3)-1)
  !***
  !*** Input:   myname - My ID in parallel machine. 0,1,2,..,P-1
  !***          pmax   - Number of procs in each coordinate direction
  !***                    P = pmax(1)*pmax(2)*pmax(3)
  !***
  !*** Output: proc - My ID in three dimensional processor mesh
  !***                1,..,pmax(1) x 1,..,pmax(2) x 1,..,pmax(3)
  !***
  !***********************************************************************

  implicit none

  integer myname, pmax(3), proc(3), remain

  proc(1) = MOD( myname, pmax(1) ) + 1
  remain  = ( myname - proc(1) + 1 )/pmax(1)
  proc(2) = MOD( remain, pmax(2) ) + 1
  proc(3) = ( remain - proc(2) + 1 )/pmax(2) + 1

  return
end subroutine MYMESH
