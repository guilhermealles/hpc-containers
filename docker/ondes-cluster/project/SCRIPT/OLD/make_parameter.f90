!       Sept 2016
!       Fabrice Dupros / Hideo Aochi
!       
!       Input
!       -        magnitude.in : xx
!        -       epicentre.in : xx xx xx 
!        -       bbox.in : x1 y1
!                          x2 y2
!        -       directory.in : input_dir
!                               output_dir



!       Output (cf. notice Ondes3D)
!       -       nice2.prm
!       -       nice2.map
!       -       nice2_station.map
!       -       nice2.hist


!       Box initial (Xmin / Xman -- Ymin / Ymax)
!       -        976 200 / 1 076 200
!       -      1 836 200 / 1 906 000


!       Taille initial du domaine
!       -        100 x 70 x 26 KM

!       Zero inital
!       -       976 000 // 1 836 000 // -25 000
!       -       x1 / y1 devient le nouveau ZERO de la simulation


!       Hypocentre initial -- posiiton source
!       -       xhypo = 1 056 330 //  yhypo = 1 877 363 //   zhypo = -7000.



!       Magnitude initiale par defaut = 6.2

!       Fichier recepteur --> nice2_station.map

!	        dsmax -- Valeur en DURE --- pour la sortie au niveau des recepteurs
!               Pour le moment on sort l'information tous les 4 KM (4 000m)

!	x0,y0,z0 Zero initial de la simulation de Nice (976 000 // 1 836 000 )
!	x1,y1,z1 coin bas gauche issu du fichier bbox.in 
!       x1,y1,z1 --> nouveau zero de la simulation dans le fichier .prm
!	x2,y2,z2 coin haut droit



!       Dans le nouveau fichier prm.
!       On considere z1 = z
!       On part de imin =1 / jmin = 1 / kmin = 1
!       On part de imax = x2-x1 / jmax = y2-y1 / kmax = z0
!       Duree - Arbitraire - de la simulation =  10 secondes
!       Pas d espace arbitraire = 400 m
!       On considere une seule source --> Dans le fichier nice.prm initial on
!       decrit 16471 sources ....
!
!       Garder un ratio entier entre le pas de temps du modele initial (200) et le pas
!       d'espace du nouveau modele (ds en parametre)
!       
!       PArametre pour la sauvegarde en surface en dur -- dsmax=1000m
!---------------------------------------------------------------------------------------------



        DOUBLE PRECISION stf(100000)
        DOUBLE PRECISION spline3, sum
        CHARACTER resp*3
 	character(len=100)::dir1
 	character(len=100)::dir2
	
	integer i1,i2,j1,j2,k,k1,k2,dsmax,cpt
        parameter (ds=400,x0=976000.,y0=1836000,z0=-25000.,xmin_ini=976000.,xmax_ini=1076000.)
        parameter (ymin_ini=1836000.,ymax_ini=1906000.,dt_ori=0.004d0,wd=300)
        parameter (tmax=10.,dsmax=4000)


        double precision ratio,dt,dtmax


        EXTERNAL spline3




!       Open des fichiers de parametres





!       Lecture des repertoires de lecture et d ecriture
        OPEN(24, FILE='directory.in')
	read(24,"(a)")	dir1
	read(24,"(a)") 	dir2
	write(*,*) 'Input directory:  ',trim(dir1)
	write(*,*) 'Output directory: ',trim(dir2)
	close(24)




!       Lecture de la nouvelle bbox

	OPEN(21, FILE='bbox.in') 
        read(21,*) x1, y1
        read(21,*) x2, y2
        close(21)
        
        zmin = z0
        z1 = z0



!       Controle des valeurs de la bbox
        if (x2.le.x1.or.y2.le.y1) then
          write(*,*) "wrong position of model frame"
          stop
        endif
        if (x1.lt.xmin_ini.or.y2.gt.ymax_ini) then
          write(*,*) "New model should be an extraction of old model"
          write(*,*) "On doit avoir new_x > ini_x (coin bas gauche)/// new_y < ini_y (coin haut droit) : ",x1,xmin_ini,y2,ymax_ini
          stop
        endif

        if (y1.lt.ymin_ini) then
          write(*,*) "New model should be an extraction of old model"
          write(*,*) "On doit avoir new_y > ini_y (coin bas gauche) : ",y1,ymin_ini
          stop
        endif


!       Recentrage de la bbox sur des points de grille existants
!       coin haut droit et coin bas gauche
!       bbbox (x1,y1 et x2,y2) 
!       on connait x0,y0 et ds - modele initial
        

        write(*,*) "Initial domain",xmin_ini,xmax_ini,ymin_ini,ymax_ini
        write(*,*) "New domain before shift",x1,x2,y1,y2
        i=0
        do while (x1.gt.(x0+ds*i))
                i=i+1
        enddo
        x1 = x0 +ds*i

        j=0
        do while (y1.gt.(y0+ds*j))
                j=j+1
        enddo
        y1 = y0 +ds*j

        i=0
        do while (x2.gt.(x0+ds*i))
                i=i+1
        enddo
        x2 = x0 +ds*i

        j=0
        do while (y2.gt.(y0+ds*j))
                j=j+1
        enddo
        y2 = y0 +ds*j

        write(*,*) "New domain after shift due to the grid",x1,x2,y1,y2



                




!       Pas d espace par defaut fixe en dur au niveau des parametres
!       Le pas de temps est adapte a partir de cette valeur


!       Position en abscisse

	imin = 1
	imax = int((x2-x1)/ds)
        xmin = ds*(imin)
        xmax = ds*(imax)



!       Position en ordonnee

	jmin=1
	jmax=int((y2-y1)/ds)
        ymin = ds*(jmin)
        ymax = ds*(jmax)

        kmin=1
        kmax=int((abs(z0)+1000)/ds)





!       Calcul du pas de temps MAX en fonction du pas d espace        
!        dt <= 0.5 ds/Vmax (say 8 km/s)
        dtmax = 0.5*ds/8000

!       Exple
!               Si dt=0.1 < dtmax on garde dt=0.1
!               Si dt=0.1 > dtmax on diminue dt avec dt=dt-0.00.5 en essayant de
!               s'approcher de dtmax progressivement 
        dt = 0.1d0
        do while ( dt.gt.0.0 ) 
          dt = dt - 0.005d0
          if (dt.lt.dtmax) goto 11
        enddo
 11     continue



!       Calcul de la duree de la simulation
!       Valeur arbitraire de 10  secondes
        itmax = int(tmax/dt) 
        write(*,*) "simulation duration = ", tmax, "number of time step = ", itmax,"time step = ",dt




!       Lecture des coordonnees de la localisation de la source
	OPEN(22, FILE='epicentre.in') 
        read(22,*) xhypo, yhypo,zhypo
	close(22)



!       Controle source dans la bbox.
        if (xhypo.le.x1.or.yhypo.le.y1) then
          write(*,*) "Mauvaise position de la source - hors de la bbox"
          stop
        endif
        if (xhypo.gt.x2.or.yhypo.gt.y2) then
          write(*,*) "Mauvaise position de la source - hors de la bbox"
          stop
        endif
        if (zhypo.lt.z1) then
          write(*,*) "Mauvaise position de la source - hors de la bbox"
          stop
        endif





!       Ajustement de la position de la source en fonction de la bbox
!       (x1,y1 //x2,y2)
!        if (xhypo.lt.(x1+wd)) xhypo = x1 + wd
!        if (xhypo.gt.(x2-wd)) xhypo = x2 - wd
!        if (yhypo.lt.(y1+wd)) yhypo = y1 + wd
!        if (yhypo.gt.(y2-wd)) yhypo = y2 - wd
!        if (zhypo.lt.(zmin+wd)) zhypo = zmin + wd


        if (xhypo.lt.(x1+wd)) x1 = xhypo  - wd
        if (xhypo.gt.(x2-wd)) x2 = xhypo  + wd
        if (yhypo.lt.(y1+wd)) y1 = yhypo  - wd
        if (yhypo.gt.(y2-wd)) y2 = yhypo  + wd
        if (zhypo.lt.(zmin+wd)) zhypo = zmin + wd

        write(*,*) "New domain (extended) after shift due to the source",x1,x2,y1,y2


        if (zhypo.gt.0.0) zhypo = -wd


!       Controle source dans la bbox.
        if (xhypo.le.x1.or.yhypo.le.y1) then
          write(*,*) "Mauvaise position de la source - hors de la bbox"
          stop
        endif
        if (xhypo.gt.x2.or.yhypo.gt.y2) then
          write(*,*) "Mauvaise position de la source - hors de la bbox"
          stop
        endif
        if (zhypo.lt.z1) then
          write(*,*) "Mauvaise position de la source - hors de la bbox"
          stop
        endif




        write(*,*) "Source location (adjusted)",xhypo, yhypo, zhypo
        ihypo = int(xhypo-x1)/ds
        jhypo = int(yhypo-y1)/ds
        khypo = int(zhypo-z1)/ds
!        write(*,*) "grid location ", ihypo, jhypo, khypo




!       Lecture de la magnitude
  	OPEN(23, FILE='magnitude.in') 
        read (23,*)  smw
	close(23)


!       Mechanisme (strike, dip, rake) defini par defaut
!       Possibilite de saisie manuelle par la suite strike (0-360), dip (0-90), rake (-90-90).
        strike = 237.
        dip = 46.
        rake = 79.


!c fault length simply following Wells and Coppersmith (1994)       
        fl = 10.**(-3.22 + 0.69 * smw)

!c source duration (roundtrip of rupture of 2.2 km/s)
        duration = 2.*fl/2.2 
        write(*,*) "Default values of 2001/02/25 M4.5 Nice earthquake"
        write(*,*) "strike = ",strike,"dip = ",dip,"rake = ",rake
        write(*,*) "magnitude = ",smw
        write(*,*) "Fault length = ", fl,"Given duration = ", duration




!        conversion seismic moment - magnitude
        sm0 = 10.**(1.5 * smw + 9.1)

        dts = 1./20.
        nmax = int(duration/dts) + 1
        tpeak = duration/2.
        dtsp = duration/4.
    
        sum = 0.0d0 

        do n=1, nmax
          t = dts* (n-1)
          stf(n) = sm0/dtsp*spline3(dble(t),dble(tpeak),dble(dtsp))
          sum = sum + stf(n)*dts
        enddo 
        write(*,'(a25, 2e12.5)') "Seismic moment (N.m) = ", sm0, sum
        write(*,*) "Duree de la source (etape de temps)",nmax

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
!c
!c OUTPUT for Ondes3D
!c

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
!c
!c	PRM File
!c
	open (31,FILE='nice2.prm')
	write(31,*) ":::::::::::: main parameters ::::::::::::::::"
	write(31,*) "<dt>",dt,"</dt>"
	write(31,*) "<ds>",ds,"</ds>"
	write(31,*) "<tmax>",itmax,"</tmax>"
	write(31,*) "<xMin>",1,"</xMin> "," <xMax>",imax,"</xMax>"
	write(31,*) "<yMin>",1,"</yMin> "," <yMax>",jmax,"</yMax>"
	write(31,*) "<zMin>",1,"</zMin> "," <zMax>",kmax,"</zMax>"

	write(31,*) "<fd0>",5.," </fd0> "
	write(31,*) ":::::::::::: output description :::::::::::::"

	write(31,*) "<dir> ",trim(dir2)," </dir>", " Directory of the Results"
	write(31,*) "# snapshots "
	write(31,*) "<i0> ",1," </i0>", " index of the output planes"
	write(31,*) "<j0> ",1," </j0>" 
	write(31,*) "<k0> ",1," </k0>"
	
	write(31,*) " # seismograms # "
	write(31,*) "<fstatMap> ",trim(dir1),"nice2_station.map",&
       &" </fstatMap>" 

	write(31,*) " ::::::::::: Source Description ::::::::::"
	write(31,*) "<fsrcMap> ",trim(dir1),"nice2.map"," </fsrcMap>"
	write(31,*) "<fsrcHist> ",trim(dir1),"nice2.hist"," </fsrcHist>"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	On prend x1,y1,z1 comme nouveau zero de la simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	write(31,*) " ::::::::::: Geology Description ::::::::::"
	write(31,*) "<fgeo> ",trim(dir1),"nice_geol_2.prm"," </fgeo>"
	write(31,*) "<x0> ",x1," </x0>"
	write(31,*) "<y0> ",y1," </y0>"
	write(31,*) "<z0> ",z1," </z0>"
	write(31,*) "<nlayer> ","6"," </nlayer>"

	write(31,*) "<layer1>"
	write(31,*) "       <name>"," out "," </name>"
	write(31,*) "       <rho>"," 1730 "," </rho>"
	write(31,*) "       <vp>"," 0.0 "," </vp>"
	write(31,*) "       <vs>"," 0.0 "," </vs>"
	write(31,*) "</layer1>"

	write(31,*) "<layer2>"
	write(31,*) "       <name>"," sea "," </name>"
	write(31,*) "       <rho>"," 1030 "," </rho>"
	write(31,*) "       <vp>"," 1530 "," </vp>"
	write(31,*) "       <vs>"," 0.0 "," </vs>"
	write(31,*) "</layer2>"

	write(31,*) "<layer3>"
	write(31,*) "       <name>"," M "," </name>"
	write(31,*) "       <rho>"," 2990 "," </rho>"
	write(31,*) "       <vp>"," 6930 "," </vp>"
	write(31,*) "       <vs>"," 400 "," </vs>"
	write(31,*) "</layer3>"

	write(31,*) "<layer4>"
	write(31,*) "       <name>"," CFV "," </name>"
	write(31,*) "       <rho>"," 2990 "," </rho>"
	write(31,*) "       <vp>"," 6930 "," </vp>"
	write(31,*) "       <vs>"," 4000 "," </vs>"
	write(31,*) "</layer4>"

	write(31,*) "<layer5>"
	write(31,*) "       <name>"," SA "," </name>"
	write(31,*) "       <rho>"," 2600 "," </rho>"
	write(31,*) "       <vp>"," 5720 "," </vp>"
	write(31,*) "       <vs>"," 3300 "," </vs>"
	write(31,*) "</layer5>"

	write(31,*) "<layer6>"
	write(31,*) "       <name>"," Cover "," </name>"
	write(31,*) "       <rho>"," 1730 "," </rho>"
	write(31,*) "       <vp>"," 3000 "," </vp>"
	write(31,*) "       <vs>"," 1400 "," </vs>"
	write(31,*) "</layer6>"


	close(31)




!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c
!c	nice.map File
!c

	open (32,FILE='nice2.map')
        write(32,'(i5)') 1
        write(32,'(3f12.1)') xhypo, yhypo, zhypo
        write(32,'(i5, 3f12.1)') 1, xhypo, yhypo, zhypo
	close (32)



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c
!c	nice.hist File
!c


	open (33,FILE='nice2.hist')
        write(33,'(a1, 2e12.5)') ">", sm0, sum
        write(33,'(a1, 2f12.3)') ">", smw, duration
        write(33,'(f12.3, f12.5)') 1.0, dts
        write(33,'(i5)') nmax
        write(33,'(i5, 3f12.3)') 1, strike, dip, rake
        write(33,*) (stf(n), n=1, nmax)
	close(33)







!C	Fichier Stations
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	open (34,FILE='nice2_station.map')

	cpt = 0
        nxgridsurf = int((x2-x1)/dsmax)
        nygridsurf = int((y2-y1)/dsmax)

        cpt = nxgridsurf*nygridsurf
        write(*,*) "D",cpt,nxgridsurf,nygridsurf
        write(34,*) cpt
        
        icpt = 0
        do i=1,nxgridsurf
                do j=1,nygridsurf
                        icpt = icpt + 1
                        x_surf = i*dsmax + x1
                        y_surf = j*dsmax + y1
                        write(34,*) icpt,x_surf,y_surf
                enddo
        enddo
        
!	i1 = dsmax*1
!	i2 = dsmax*imax
!	j1 = dsmax*1
!	j2 = dsmax*jmax
!	
!	cpt = 0
!	Do i=i1,i2,dsmax
!		Do j=j1,j2,dsmax
!		cpt = cpt + 1
!	Enddo
!	Enddo
	write(*,*) "Default number of surface receivers",cpt
!	write(34,*) cpt
!
!	cpt = 0
!	k=0
!	Do i=i1,i2,dsmax
!		Do j=j1,j2,dsmax
!		cpt = cpt + 1
!		write(34,*) cpt,i,j,0
!	Enddo
!	Enddo

	close(34)



        ratio = dt / dt_ori
        write(*,*) "////  RATIO temps ////",ratio






        stop
        end

	DOUBLE PRECISION FUNCTION spline3(x, x0, dx)
	DOUBLE PRECISION x, x0, dx, absx

	absx = abs(x-x0)
	if( absx.ge.(2.*dx) ) then
	  spline3 = 0.0d0
	else if ( absx.lt.(2.*dx).and.absx.ge.dx) then
	  spline3 = (2.*dx - absx)**3/(6.*dx**3)
	else
	  spline3 = 3.*absx**3 - 6.*dx*absx**2 + 4.*dx**3
	  spline3 = spline3/(6.*dx**3)
	endif

	return
	END

