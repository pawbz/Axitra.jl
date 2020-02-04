module Axitra

using ProgressMeter


function Axitra.update!()


nc=1; # number of components?
nfreq=512; # number of freqs
tl=70.; # total length in time?
aw=2.; #??
nr=5; #number of receivers #
ns=1; # number of source
xl=1000000.; # ??
ikmax=100000; # ??

dfreq = inv(tl)
# latlon=.true.,freesurface=.true.,sourcefile="source",statfile="station"



	# loop over frequencies
@showprogress	for jf=1:nfreq
      freq = (jf - 1)*dfreq
#      println("freq", jf, " ", nfreq)
rw = π * π * freq 
      omega = complex(rw, aw)
      omega2 = omega*omega
      a1 = 0.5 *inv(omega2)*inv(xl) # ??
      zom = sqrt(rw*rw + aw*aw)

      if (jf == 1) then
         phi = -π/2
      else
         phi = atan(aw/rw) # phase?
      end
      #=

      do ir = 1, nr
      do is = 1, ns
         tconv(ir, is) = .false.
      enddo
      enddo

      ttconv = .false.
#define OLDEF
#ifdef OLDEF
      xlnf = (ai*phi + dlog(zom))/pi
#else
! Futterman
      xlnf = (ai*phi + dlog(zom/(pi2*fref)))
! Kjartansson
      xlnf = zom/(pi2*fref)
#endif

!            ******************************************
!            ******************************************
!            **  RESOLUTION PAR BOUCLE EXTERNE EN Kr **
!            ******************************************
!            ******************************************

      do ik = 0, ikmax

	      kr = (ik + .258)*pil #  radial wavenumber, Herrmann & Mandal (1986)
         kr2 = kr*kr

!+++++++++++++
!              Calcul de nombreux coefficients et des fonctions de Bessel
!+++++++++++++

         call reflect0(ik + 1, nc, nr, ns, nrs)

!+++++++++++++
!              Calcul des coefficients de reflexion/transmission
!               Matrice de Reflection/Transmission et Dephasage
!+++++++++++++

         call reflect1(freesurface, nc)

!+++++++++++++
!              Calcul des matrices de reflectivite : mt(),mb(),nt(),nb()
!              (rapport des differents potentiels montant/descendant
!                        en haut et en bas de chaque couche)
!+++++++++++++

         call reflect2(nc)

!+++++++++++++
!               Calcul des matrices de passage des vecteurs potentiel
!                source, aux vecteurs potentiel PHI, PSI et KHI au sommet
!                de chaque couche
!+++++++++++++

         call reflect3(ncs)

!+++++++++++++
!               Calcul des potentiels et des deplacement dus aux sources du
!                tenseur, en chaque recepteur (termes en kr, r, z)
!+++++++++++++

         call reflect4(jf, ik, ik .gt. ikmin, tconv, nc, nr, ns, ncs, ncr)

         if (ttconv) exit

      enddo ! wavenumber loop

!+++++++++++++
!               Calcul des deplacements aux recepteurs
!
!                Sortie des resultats
!+++++++++++++
      lastik = ik - 1
      write (out, *) 'freq =', freq, 'iter =', lastik
      write (6,"(1a1,'freq ',I5,'/',I5,' iter=',I10,$)") char(13),jf, nfreq,lastik

      if (jf .eq. 1) lastik = 0

      call reflect5(jf, nr, ns)

      if (ik .ge. ikmax) then
         write (6, *) 'Depassement du nombre d iteration maximum'
         stop
      endif

   enddo ! frequency loop
!$OMP END PARALLEL
   write(6,*) 'Done'
   stop

=#

	end

end


end # module
