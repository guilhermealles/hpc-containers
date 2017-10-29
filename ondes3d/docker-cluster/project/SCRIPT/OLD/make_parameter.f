c program make parameters for Ondes3d
c 2016/06/15 H. Aochi
        DOUBLE PRECISION stf(10000)
        DOUBLE PRECISION spline3, sum
        CHARACTER resp*3
        EXTERNAL spline3


c NOTE: (0,0,0) of simulaiton fixed at (x0, y0, z0) in nice_geol.prm       
        x0 = 976000.
        y0 = 1836000.
        z0 = -25000.
c maximum domain 500 x 350 x 180

c
c INPUT
c

        write(*,*) "lower, left corner ...  in m"
c        read(*,*) x1, y1
        x1 = x0
c978000.
        y1 = y0
c1840000.

        write(*,*) "Upper, right corner ... in m"
c        read(*,*) x2, y2
        x2 = 1006000
c1020000.
        y2 = 1866000.

c depth limit fixed at 20 km      
        zmin = z0
c20000.

c check     
        if (x2.le.x1.or.y2.le.y1) then
          write(*,*) "wrong position of model frame"
          stop
        endif

        
c
c ADJUSTEMENT OF MODEL PARAMETER
c
c apply default value (ds = 500 m)
        ds = 500.
        imin = int((x1-x0)/ds)
        imax = int((x2-x0)/ds)
        xmin = ds*(imin)
        xmax = ds*(imax)
        write(*,*) "along X", x2-x1, " (m) ->", imin, "to", imax

        jmin = int((y1-y0)/ds)
        jmax = int((y2-y0)/ds)
        ymin = ds*(jmin)
        ymax = ds*(jmax)
        write(*,*) "along Y", y2-y1, "(m) ->", jmin, "to", jmax

        kmin = int((zmin-z0)/ds)
        kmax = int((1000.-z0)/ds)
c	kmax = int((-z0)/ds)
        write(*,*) "along Z", zmin, "(m) ->", kmin, "to", kmax

        
c dt <= 0.5 ds/Vmax (say 8 km/s)
        dtmax = 0.5*ds/8000.
        write(*,'(a31, f12.5)') "possible time step = ", dtmax
        dt = 0.1
        
        do while ( dt.gt.0.0 ) 
          dt = dt - 0.005
          if (dt.lt.dtmax) goto 11
        enddo
 11     continue
        write(*,'(a31, f12.5)') "proposed time step = ", dt

c run a simulation for 10 sec
        tmax = 10.
        itmax = int(tmax/dt) 
        write(*,*) "simulation duration = ", tmax, "time step = ", itmax

c
c source
c
        write(*,*)
       write(*,*) "Source position ... x, y, z in m"
c       read(*,*) xhypo, yhypo, zhypo
        xhypo = 1014949.
        yhypo = 145168.
        zhypo = -5000.
        wd = 3000.
        if (xhypo.lt.(x1+wd)) xhypo = x1 + wd
        if (xhypo.gt.(x2-wd)) xhypo = x2 - wd
        if (yhypo.lt.(y1+wd)) yhypo = y1 + wd
        if (yhypo.lt.(y2-wd)) yhypo = y2 - wd
        if (zhypo.lt.(zmin+wd)) zhypo = zmin + wd
        if (zhypo.gt.0.0) zhypo = -wd

        write(*,*) "ajusted if needed (avoid model boundaries)..."
        write(*,*) xhypo, yhypo, zhypo
        ihypo = int(xhypo-x0)/ds
        jhypo = int(yhypo-y0)/ds
        khypo = int(zhypo-z0)/ds
        write(*,*) "grid number = ", ihypo, jhypo, khypo

        write(*,*)
        write(*,*) "Magnitude (0-7, for example)..."
c        read (*,*)  smw
        smw = 4.5        
c        smw = 6.0

        write(*,*) "Mecnaism ...."
        strike = 237.
        dip = 46.
        rake = 79.
        write(*,*) "default values of 2001/02/25 M4.5 Nice earthquake"
        write(*,'(5x, 3f12.2)') strike, dip, rake
        write(*,*) "Is this OK? (Y/N)" 
        read (*,*) resp
        if (resp.ne."Y".and.resp.ne."y") then
          write(*,*) 
     &  "input strike (0-360), dip (0-90), rake (-90-90)..."
          read(*,*) strike, dip, rake
        endif
        

c fault length simply following Wells and Coppersmith (1994)       
        fl = 10.**(-3.22 + 0.69 * smw)

c source duration (roundtrip of rupture of 2.2 km/s)
        duration = 2.*fl/2.2 

        write(*,*) "Magnitude = ", smw
        write(*,*) "Fault length = ", fl, 
     &          "Given duration = ", duration

c conversion seismic moment - magnitude
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
c
c OUTPUT for Ondes3D
c
        write(*,*)
        write(31,'(f12.3)') dt
        write(31,'(f12.3)') ds
        write(31,'(i5)') itmax
        write(31,'(2i5)') 0, imax
        write(31,'(2i5)') 0, jmax
        write(31,'(2i5)') kmax, 0
        write(*,*) "prm file was written in 31." 

        write(32,'(i5)') 1
        write(32,'(3f12.1)') xhypo, yhypo, zhypo
        write(32,'(i5, 3f12.1)') 1, xhypo, yhypo, zhypo
        write(*,*) "src file was written in 32."

        write(33,'(a1, 2e12.5)') ">", sm0, sum
        write(33,'(a1, 2f12.3)') ">", smw, duration
        write(33,'(f12.3, f12.5)') 1.0, dts
        write(33,'(i5)') nmax
        write(33,'(i5, 3f12.3)') 1, strike, dip, rake
        write(33,*) (stf(n), n=1, nmax)
        write(*,*) "hist file was written in 33."

        stop
        end
c
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

