#define SPONGE

Module agrif_opa_sponge
#if defined key_agrif  && ! defined key_offline
   USE par_oce
   USE oce
   USE dom_oce
   USE in_out_manager
   USE agrif_oce
   USE wrk_nemo  

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_Sponge_Tra, Agrif_Sponge_Dyn, interptsn, interpun, interpvn

   !!----------------------------------------------------------------------
   !! NEMO/NST 3.3 , NEMO Consortium (2010)
   !! $Id: agrif_opa_sponge.F90 3186 2011-11-27 08:16:19Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   CONTAINS

   SUBROUTINE Agrif_Sponge_Tra
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Sponge_Tra ***
      !!---------------------------------------------
#include "domzgr_substitute.h90"
      !!
      INTEGER :: ji,jj,jk,jn
      INTEGER :: spongearea
      REAL(wp) :: timecoeff
      REAL(wp) :: ztsa, zabe1, zabe2, zbtr
      REAL(wp), POINTER, DIMENSION(:,:    ) :: localviscsponge
      REAL(wp), POINTER, DIMENSION(:,:    ) :: ztu, ztv
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: ztab
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: tsbdiff

#if defined SPONGE
      CALL wrk_alloc( jpi, jpj, localviscsponge, ztu, ztv )
      CALL wrk_alloc( jpi, jpj, jpk, jpts, ztab, tsbdiff  )

      timecoeff = REAL(Agrif_NbStepint(),wp)/Agrif_rhot()

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = .TRUE.
      ztab = 0.e0
      CALL Agrif_Bc_Variable(ztab, tsa_id,calledweight=timecoeff,procname=interptsn)
      Agrif_UseSpecialValue = .FALSE.

      tsbdiff(:,:,:,:) = tsb(:,:,:,:) - ztab(:,:,:,:)

      spongearea = 2 + 2 * Agrif_irhox()

      localviscsponge = 0.
      
      IF (.NOT. spongedoneT) THEN
         spe1ur(:,:) = 0.
         spe2vr(:,:) = 0.

      IF ((nbondi == -1).OR.(nbondi == 2)) THEN
         DO ji = 2, spongearea
            localviscsponge(ji,:) = visc_tra * (spongearea-ji)/real(spongearea-2)
         ENDDO
	 
	 spe1ur(2:spongearea-1,:)=0.5 * (localviscsponge(2:spongearea-1,:) + localviscsponge(3:spongearea,:)) &
	       * e2u(2:spongearea-1,:) / e1u(2:spongearea-1,:)

         spe2vr(2:spongearea,1:jpjm1) = 0.5 * (localviscsponge(2:spongearea,1:jpjm1) + &
	          localviscsponge(2:spongearea,2:jpj)) &
	        * e1v(2:spongearea,1:jpjm1) / e2v(2:spongearea,1:jpjm1)
      ENDIF

      IF ((nbondi == 1).OR.(nbondi == 2)) THEN
         DO ji = nlci-spongearea + 1,nlci-1
            localviscsponge(ji,:) = visc_tra * (ji - (nlci-spongearea+1))/real(spongearea-2)
         ENDDO
	 
	 spe1ur(nlci-spongearea + 1:nlci-2,:)=0.5 * (localviscsponge(nlci-spongearea + 1:nlci-2,:) + &
	        localviscsponge(nlci-spongearea + 2:nlci-1,:)) &
	       * e2u(nlci-spongearea + 1:nlci-2,:) / e1u(nlci-spongearea + 1:nlci-2,:)

         spe2vr(nlci-spongearea + 1:nlci-1,1:jpjm1) = 0.5 * (localviscsponge(nlci-spongearea + 1:nlci-1,1:jpjm1) &
	           + localviscsponge(nlci-spongearea + 1:nlci-1,2:jpj)) &
	        * e1v(nlci-spongearea + 1:nlci-1,1:jpjm1) / e2v(nlci-spongearea + 1:nlci-1,1:jpjm1)
      ENDIF


      IF ((nbondj == -1).OR.(nbondj == 2)) THEN
         DO jj = 2, spongearea
            localviscsponge(:,jj) = visc_tra * (spongearea-jj)/real(spongearea-2)
         ENDDO
	 
	 spe1ur(1:jpim1,2:spongearea)=0.5 * (localviscsponge(1:jpim1,2:spongearea) + &
	        localviscsponge(2:jpi,2:spongearea)) &
	       * e2u(1:jpim1,2:spongearea) / e1u(1:jpim1,2:spongearea)

         spe2vr(:,2:spongearea-1) = 0.5 * (localviscsponge(:,2:spongearea-1) + &
	          localviscsponge(:,3:spongearea)) &
	        * e1v(:,2:spongearea-1) / e2v(:,2:spongearea-1)
      ENDIF

      IF ((nbondj == 1).OR.(nbondj == 2)) THEN
         DO jj = nlcj-spongearea + 1,nlcj-1
            localviscsponge(:,jj) = visc_tra * (jj - (nlcj-spongearea+1))/real(spongearea-2)
         ENDDO
	 
	 spe1ur(1:jpim1,nlcj-spongearea + 1:nlcj-1)=0.5 * (localviscsponge(1:jpim1,nlcj-spongearea + 1:nlcj-1) + &
	         localviscsponge(2:jpi,nlcj-spongearea + 1:nlcj-1)) &
	       * e2u(1:jpim1,nlcj-spongearea + 1:nlcj-1) / e1u(1:jpim1,nlcj-spongearea + 1:nlcj-1)

         spe2vr(:,nlcj-spongearea + 1:nlcj-2) = 0.5 * (localviscsponge(:,nlcj-spongearea + 1:nlcj-2) + &
	         localviscsponge(:,nlcj-spongearea + 2:nlcj-1)) &
	        * e1v(:,nlcj-spongearea + 1:nlcj-2) / e2v(:,nlcj-spongearea + 1:nlcj-2)
      ENDIF
      
         spbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:))

         spongedoneT = .TRUE.
      ENDIF

      DO jn = 1, jpts
         DO jk = 1, jpkm1
            !
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zabe1 = umask(ji,jj,jk) * spe1ur(ji,jj) * fse3u(ji,jj,jk)
                  zabe2 = vmask(ji,jj,jk) * spe2vr(ji,jj) * fse3v(ji,jj,jk)
                  ztu(ji,jj) = zabe1 * ( tsbdiff(ji+1,jj  ,jk,jn) - tsbdiff(ji,jj,jk,jn) )
                  ztv(ji,jj) = zabe2 * ( tsbdiff(ji  ,jj+1,jk,jn) - tsbdiff(ji,jj,jk,jn) )
               ENDDO
            ENDDO

            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zbtr = spbtr2(ji,jj) / fse3t(ji,jj,jk)
                  ! horizontal diffusive trends
                  ztsa = zbtr * (  ztu(ji,jj) - ztu(ji-1,jj  )   &
                  &              + ztv(ji,jj) - ztv(ji  ,jj-1)  )
                  ! add it to the general tracer trends
                  tsa(ji,jj,jk,jn) = tsa(ji,jj,jk,jn) + ztsa
               END DO
            END DO
            !
         ENDDO
      ENDDO

      CALL wrk_dealloc( jpi, jpj, localviscsponge, ztu, ztv )
      CALL wrk_dealloc( jpi, jpj, jpk, jpts, ztab, tsbdiff  )
#endif

   END SUBROUTINE Agrif_Sponge_Tra

   SUBROUTINE Agrif_Sponge_dyn
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Sponge_dyn ***
      !!---------------------------------------------
#include "domzgr_substitute.h90"
      !!
      INTEGER :: ji,jj,jk
      INTEGER :: spongearea
      REAL(wp) :: timecoeff
      REAL(wp) :: ze2u, ze1v, zua, zva, zbtr
      REAL(wp), POINTER, DIMENSION(:,:) :: localviscsponge
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ubdiff, vbdiff
      REAL(wp), POINTER, DIMENSION(:,:,:) :: rotdiff, hdivdiff
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztab

#if defined SPONGE
      CALL wrk_alloc( jpi, jpj, localviscsponge )
      CALL wrk_alloc( jpi, jpj, jpk, ztab, ubdiff, vbdiff, rotdiff, hdivdiff )

      timecoeff = REAL(Agrif_NbStepint(),wp)/Agrif_rhot()

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn
      ztab = 0.e0
      CALL Agrif_Bc_Variable(ztab, ua_id,calledweight=timecoeff,procname=interpun)
      Agrif_UseSpecialValue = .FALSE.

      ubdiff(:,:,:) = (ub(:,:,:) - ztab(:,:,:))*umask(:,:,:)

      ztab = 0.e0
      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn
      CALL Agrif_Bc_Variable(ztab, va_id,calledweight=timecoeff,procname=interpvn)
      Agrif_UseSpecialValue = .FALSE.

      vbdiff(:,:,:) = (vb(:,:,:) - ztab(:,:,:))*vmask(:,:,:)

      spongearea = 2 + 2 * Agrif_irhox()

      localviscsponge = 0.

      IF (.NOT. spongedoneU) THEN
         spe1ur2(:,:) = 0.
         spe2vr2(:,:) = 0.

      IF ((nbondi == -1).OR.(nbondi == 2)) THEN
         DO ji = 2, spongearea
            localviscsponge(ji,:) = visc_dyn * (spongearea-ji)/real(spongearea-2)
         ENDDO
	 
	 spe1ur2(2:spongearea-1,:)=0.5 * (localviscsponge(2:spongearea-1,:) + localviscsponge(3:spongearea,:))

         spe2vr2(2:spongearea,1:jpjm1) = 0.5 * (localviscsponge(2:spongearea,1:jpjm1) + &
	          localviscsponge(2:spongearea,2:jpj))
      ENDIF

      IF ((nbondi == 1).OR.(nbondi == 2)) THEN
         DO ji = nlci-spongearea + 1,nlci-1
            localviscsponge(ji,:) = visc_dyn * (ji - (nlci-spongearea+1))/real(spongearea-2)
         ENDDO
	 
	 spe1ur2(nlci-spongearea + 1:nlci-2,:)=0.5 * (localviscsponge(nlci-spongearea + 1:nlci-2,:) + &
	        localviscsponge(nlci-spongearea + 2:nlci-1,:))

         spe2vr2(nlci-spongearea + 1:nlci-1,1:jpjm1) = 0.5 * (localviscsponge(nlci-spongearea + 1:nlci-1,1:jpjm1) &
	           + localviscsponge(nlci-spongearea + 1:nlci-1,2:jpj))
      ENDIF


      IF ((nbondj == -1).OR.(nbondj == 2)) THEN
         DO jj = 2, spongearea
            localviscsponge(:,jj) = visc_dyn * (spongearea-jj)/real(spongearea-2)
         ENDDO
	 
	 spe1ur2(1:jpim1,2:spongearea)=0.5 * (localviscsponge(1:jpim1,2:spongearea) + &
	        localviscsponge(2:jpi,2:spongearea))

         spe2vr2(:,2:spongearea-1) = 0.5 * (localviscsponge(:,2:spongearea-1) + &
	          localviscsponge(:,3:spongearea))
      ENDIF

      IF ((nbondj == 1).OR.(nbondj == 2)) THEN
         DO jj = nlcj-spongearea + 1,nlcj-1
            localviscsponge(:,jj) = visc_dyn * (jj - (nlcj-spongearea+1))/real(spongearea-2)
         ENDDO
	 
	 spe1ur2(1:jpim1,nlcj-spongearea + 1:nlcj-1)=0.5 * (localviscsponge(1:jpim1,nlcj-spongearea + 1:nlcj-1) + &
	         localviscsponge(2:jpi,nlcj-spongearea + 1:nlcj-1))

         spe2vr2(:,nlcj-spongearea + 1:nlcj-2) = 0.5 * (localviscsponge(:,nlcj-spongearea + 1:nlcj-2) + &
	         localviscsponge(:,nlcj-spongearea + 2:nlcj-1))
      ENDIF

         spongedoneU = .TRUE.
	 
	  spbtr3(:,:) = 1./( e1f(:,:) * e2f(:,:))
      ENDIF
      
      IF (.NOT. spongedoneT) THEN
        spbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:))      
      ENDIF
      
      DO jk=1,jpkm1
      ubdiff(:,:,jk) = ubdiff(:,:,jk) * spe1ur2(:,:)
      vbdiff(:,:,jk) = vbdiff(:,:,jk) * spe2vr2(:,:)
      ENDDO
      
      hdivdiff = 0.
      rotdiff = 0.

      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         !                                             ! --------
         ! Horizontal divergence                       !   div
         !                                             ! --------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zbtr = spbtr2(ji,jj) / fse3t(ji,jj,jk)
               hdivdiff(ji,jj,jk) =   &
                  (  e2u(ji,jj)*fse3u(ji,jj,jk) * & 
                  ubdiff(ji,jj,jk) - e2u(ji-1,jj  )* &
                  fse3u(ji-1,jj  ,jk)  * ubdiff(ji-1,jj  ,jk)       &
                  + e1v(ji,jj)*fse3v(ji,jj,jk) * &
                  vbdiff(ji,jj,jk) - e1v(ji  ,jj-1)* &
                  fse3v(ji  ,jj-1,jk)  * vbdiff(ji  ,jj-1,jk)  ) * zbtr
            END DO
         END DO

         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               zbtr = spbtr3(ji,jj) * fse3f(ji,jj,jk)
               rotdiff(ji,jj,jk) = (  e2v(ji+1,jj  ) * vbdiff(ji+1,jj  ,jk) - e2v(ji,jj) * vbdiff(ji,jj,jk)    &
                  &              - e1u(ji  ,jj+1) * ubdiff(ji  ,jj+1,jk) + e1u(ji,jj) * ubdiff(ji,jj,jk)  ) &
                  &           * fmask(ji,jj,jk) * zbtr
            END DO
         END DO

      ENDDO

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ze2u = rotdiff (ji,jj,jk)
               ze1v = hdivdiff(ji,jj,jk)
               ! horizontal diffusive trends
               zua = - ( ze2u - rotdiff (ji,jj-1,jk)) / ( e2u(ji,jj) * fse3u(ji,jj,jk) )   &
                  + ( hdivdiff(ji+1,jj,jk) - ze1v      &
                  ) / e1u(ji,jj)

               zva = + ( ze2u - rotdiff (ji-1,jj,jk)) / ( e1v(ji,jj) * fse3v(ji,jj,jk) )   &
                  + ( hdivdiff(ji,jj+1,jk) - ze1v    &
                  ) / e2v(ji,jj)

               ! add it to the general momentum trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, localviscsponge )
      CALL wrk_dealloc( jpi, jpj, jpk, ztab, ubdiff, vbdiff, rotdiff, hdivdiff )

#endif

   END SUBROUTINE Agrif_Sponge_dyn

   SUBROUTINE interptsn(tabres,i1,i2,j1,j2,k1,k2,n1,n2)
      !!---------------------------------------------
      !!   *** ROUTINE interptsn ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"       
      
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,n1,n2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) :: tabres

      tabres(i1:i2,j1:j2,k1:k2,n1:n2) = tsn(i1:i2,j1:j2,k1:k2,n1:n2)

   END SUBROUTINE interptsn

   SUBROUTINE interpun(tabres,i1,i2,j1,j2,k1,k2)
      !!---------------------------------------------
      !!   *** ROUTINE interpun ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"       
      
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres

      tabres(i1:i2,j1:j2,k1:k2) = un(i1:i2,j1:j2,k1:k2)

   END SUBROUTINE interpun

   SUBROUTINE interpvn(tabres,i1,i2,j1,j2,k1,k2)
      !!---------------------------------------------
      !!   *** ROUTINE interpvn ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"       
      
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres

      tabres(i1:i2,j1:j2,k1:k2) = vn(i1:i2,j1:j2,k1:k2)

   END SUBROUTINE interpvn

#else
CONTAINS

   SUBROUTINE agrif_opa_sponge_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_OPA_sponge_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_opa_sponge : You should not have seen this print! error?'
   END SUBROUTINE agrif_opa_sponge_empty
#endif

END MODULE agrif_opa_sponge
