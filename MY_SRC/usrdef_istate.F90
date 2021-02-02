MODULE usrdef_istate
   !!======================================================================
   !!                   ***  MODULE  usrdef_istate   ***
   !!
   !!                     ===  Constant TS configuration  ===    ! VAL CHANGED, BEFORE IT WAS WRITEN: GYRE configuration
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2016-03  (S. Flavoni, G. Madec) Original code     ! VAL CHANGED, BEFORE IT WAS WRITTEN :  History :  4.0 ! 2016-03  (S. Flavoni) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce , ONLY : glamt  ! VAL LINE ADDED
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called in istate.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Set a constant T&S, to use for testing when ln_2d=False
      !!                (when it is True, T&S are automatically set to 0)       ! VAL CHANGED, BEFORE IT WAS WRITTEN : Here GYRE configuration example : (double gyre with rotated domain)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity   field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      INTEGER  ::   jk     ! dummy loop indices            ! VAL BEFORE IT WAS : INTEGER :: ji, jj, jk  ! dummy loop indices
      REAL(wp) ::   zdam   ! location of dam [Km]          ! VAL LINE ADDED
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : ' ! VAL CHANGED, I REMOVED THIS : analytical definition of initial state '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest, with a constant temperature and salinity' ! VAL CHANGED, I REMOVED THIS:  horizontally uniform T and S profiles'
      !
      ! VAL LINES ADDED
      !  rn_a0 =  0.2   !  thermal expension coefficient (nn_eos= 1)
      !  rho = rau0 - rn_a0 * (T-10)
      !  delta_T = 25 degrees  ==>>  delta_rho = 25 * rn_a0 = 5 kg/m3
      ! VAL END OF ADDED LINES
      !
      pu  (:,:,:) = 0._wp        ! ocean at rest
      pv  (:,:,:) = 0._wp
      pssh(:,:)   = 0._wp
      !
      ! VAL COMMENT LINE DO jk = 1, jpk             ! horizontally uniform T & S profiles
      ! VAL COMMENT LINE    DO jj = 1, jpj
      ! VAL COMMENT LINE       DO ji = 1, jpi
      ! VAL COMMENT LINE          pts(ji,jj,jk,jp_tem) =  (  (  16. - 12. * TANH( (pdept(ji,jj,jk) - 400) / 700 ) )   &
      ! VAL COMMENT LINE               &           * (-TANH( (500. - pdept(ji,jj,jk)) / 150. ) + 1.) / 2.             &
      ! VAL COMMENT LINE               &           + ( 15. * ( 1. - TANH( (pdept(ji,jj,jk)-50.) / 1500.) )            &
      ! VAL COMMENT LINE               &           - 1.4 * TANH((pdept(ji,jj,jk)-100.) / 100.)                        &
      ! VAL COMMENT LINE               &           + 7.  * (1500. - pdept(ji,jj,jk) ) / 1500.)                        &
      ! VAL COMMENT LINE               &           * (-TANH( (pdept(ji,jj,jk) - 500.) / 150.) + 1.) / 2.  ) * ptmask(ji,jj,jk)

      ! VAL COMMENT LINE          pts(ji,jj,jk,jp_sal) =  (  (  36.25 - 1.13 * TANH( (pdept(ji,jj,jk) - 305) / 460 ) )  &
      ! VAL COMMENT LINE               &         * (-TANH((500. - pdept(ji,jj,jk)) / 150.) + 1.) / 2                  &
      ! VAL COMMENT LINE               &         + ( 35.55 + 1.25 * (5000. - pdept(ji,jj,jk)) / 5000.                 &
      ! VAL COMMENT LINE               &         - 1.62 * TANH( (pdept(ji,jj,jk) - 60.  ) / 650. )                    &
      ! VAL COMMENT LINE               &         + 0.2  * TANH( (pdept(ji,jj,jk) - 35.  ) / 100. )                    &
      ! VAL COMMENT LINE               &         + 0.2  * TANH( (pdept(ji,jj,jk) - 1000.) / 5000.) )                  &
      ! VAL COMMENT LINE               &         * (-TANH( (pdept(ji,jj,jk) - 500.) / 150.) + 1.) / 2  ) * ptmask(ji,jj,jk)
      ! VAL COMMENT LINE       END DO
      ! VAL COMMENT LINE    END DO
      ! VAL COMMENT LINE   END DO
      !                          ! T & S profiles ! VAL LINE ADDED
      zdam = 32.                      ! density front position in kilometers ! VAL LINE ADDED
      pts(:,:,:,jp_tem) = 10._wp * ptmask(:,:,:) ! VAL LINE ADDED
      pts(:,:,:,jp_sal) = 35._wp * ptmask(:,:,:) ! VAL LINE ADDED
      !   
   END SUBROUTINE usr_def_istate

   !!======================================================================
END MODULE usrdef_istate
