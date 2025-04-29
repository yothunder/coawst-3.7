      MODULE mod_coupler_iounits
!
!svn $Id: mod_coupler_iounits.F jcwarner $
!=======================================================================
!                                                                      !
!  stdinp      Unit number for standard input (often 5).               !
!  stdout      Unit number for standard output (often 6).              !
!  Aname       Atmosphere model stadard input file name.               !
!  IWBNDname   Input boundary data file name for InWave model          !
!  IWSWNname   Input spectral SWAN data file name for InWave model     !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE swan_iounits
      USE mod_param
      implicit none
!      integer, parameter :: IOnamesize = 160
!      character (len=IOnamesize) :: Wname
! this flag is temporary for SCRIP option.
      integer :: scrip_opt
      character (len=80) :: SCRIPname
      character (len=80), dimension(:,:), pointer :: O2Wname
      character (len=80), dimension(:,:), pointer :: W2Oname
      character (len=80), dimension(:,:), pointer :: O2Aname
      character (len=80), dimension(:,:), pointer :: A2Oname
      character (len=80), dimension(:,:), pointer :: W2Aname
      character (len=80), dimension(:,:), pointer :: A2Wname
      character (len=80) :: Aname
      CONTAINS
      SUBROUTINE allocate_coupler_iounits
!=======================================================================
!                                                                      !
!  This routine initialize all the coupler io vars.                    !
!                                                                      !
!=======================================================================
      character (len=1 ), parameter :: blank = ' '
      integer :: i, io, ia, iw, ih
      allocate (O2Wname(Nocn_grids,Nwav_grids))
      allocate (W2Oname(Nwav_grids,Nocn_grids))
      DO io=1,Nocn_grids
        DO iw=1,Nwav_grids
          DO i=1,LEN(O2Wname(io,iw))
            O2Wname(io,iw)(i:i)=blank
            W2Oname(iw,io)(i:i)=blank
          END DO
        END DO
      END DO
      allocate (O2Aname(Nocn_grids,Natm_grids))
      allocate (A2Oname(Natm_grids,Nocn_grids))
      DO ia=1,Natm_grids
        DO io=1,Nocn_grids
          DO i=1,LEN(O2Aname(io,ia))
            O2Aname(io,ia)(i:i)=blank
            A2Oname(ia,io)(i:i)=blank
          END DO
        END DO
      END DO
      allocate (W2Aname(Nwav_grids,Natm_grids))
      allocate (A2Wname(Natm_grids,Nwav_grids))
      DO ia=1,Natm_grids
        DO iw=1,Nwav_grids
          DO i=1,LEN(A2Wname(iw,ia))
            A2Wname(ia,iw)(i:i)=blank
            W2Aname(iw,ia)(i:i)=blank
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE allocate_coupler_iounits
      END MODULE mod_coupler_iounits
