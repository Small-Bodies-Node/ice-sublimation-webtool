       Program cgifastrot
c
c   This program calculates the average sublimation per unit area for
c   a rapidly rotating cometary nucleus.  For a
C   sufficiently rapid rotation, or equivalently for sufficiently high
C   thermal inertia, a parallel of latitude is an isotherm and this
c   is assumed by the program.
C
C      PROGRAM ITERATES ENERGY BALANCE EQUATION BY NEWTON-RAPHSON METHOD
C      TO GET EQUILIBRIUM TEMPERATURE.  SIMPSON'S RULE IS USED TO
C      INTEGRATE OVER LATITUDE.
C
C   The properties of the ice are handled in the separate subroutine
C   sublime which provides the vapor pressure and latent heat, as well
C   as their derivatives, for a given temperature.
C
c
c   sb is the set of values of sin(b) at which the sublimation is calculated.
c   b is latitude, frac is the effective value of cos(theta) for each latitude
c  In this version, there are 41 steps in latitude (variable nb in the
c    data statement).
c
c  Modification History
c    B. Prager 06/24
c      - Original fastrot.f renamed to cgifastrot.f to make a distinction.
c      - Removed file i/o from cgifastrot.f to avoid permissions issues with cgi scripts.
c      - Added command line input of parameters. Parameter 1: Species. Parameter 2: Visual
c        Albedo. Parameter 3: Infared Albedo. Parameter 4: Heliocentric Distance. Parameter 5:
c        inclination.
c      - Initialized the parameters such that they can accept nine characters of data.

       dimension sb(41), b(41), frac(41), Z(41) 
       data pi/3.1415926/, nb/41/, q/1./
       real incl
       character*1024 prefix
c   Variables to store the read in parameters, and the type converted parameters.
       character(len=1) :: index_in
       character(len=9) :: a0_in,a1_in
       character(len=9) :: rh_in
       character(len=9) :: inc_in
       real*4 a0,a1,rh,inc
       integer index
c   various constants
       data sigma/5.67e-5/, F0/1.39e6/, dynmm/1.333e3/, 
     :   boltz/1.38e-16/, ergcal/6.953e-17/, proton/1.67e-24/
c
c   Common block subl transfers data between the main program and the
c   subroutine sublime.  Index specifies the ice under consideration and
c   should be specified in the main program together with the temperature.
c   All other parameters are set in the subroutine.
c
       character*8, spec
       double precision press, pprim, prec
       real mass
       common /subl/ press, pprim, index, temp, spec, mass,xlt,xltprim
c
c   Calculate equally spaced values of sin(b)
c   Initialize the values of latitude from the values of sin(b)
c
       prec = 0
       delsb = 2.0/(nb-1)
       sb(1) = -1.
       b(1) = asin(sb(1))
       do 5 i=2,nb
         sb(i) = sb(1) + (i-1)*delsb
         b(i)=asin(sb(i))
   5   continue
c
c  Read in command line arguments, convert their type to real or integer, and display values to screen.
c
      call getenv('HOME',prefix)

      CALL GetArg(1,index_in)
      CALL GetArg(2,a0_in)
      CALL GetArg(3,a1_in)
      CALL GetArg(4,rh_in)
      CALL GetArg(5,inc_in)

      Read(index_in, '(i10)' )  index
      Read(a0_in,*)  a0
      Read(a1_in,*)  a1
      Read(rh_in,*)  rh
      Read(inc_in,*)  inc

c   Print the initial conditions.

      write (6,*) ' Ind = ',index,'  Avis = ',a0,'  Air = ',a1,
     &  '  r_H = ',rh,'  Incl = ',inc
  10  continue
      if (a0.lt.0.)  go to 800
      incl = (90.-inc)*pi/180.
c
c   Initialize things for the chosen species
c     subroutine picks a starting temperature dependent on the species
c     if temp <= 0.
c

      temp = -1.
      call sublime
      root = 1./sqrt(mass*2.*pi*boltz)

c
c   Loop on latitude
c

      nflag = 1
      do 600 n=1,nb
  100   continue
        tp = 0.
        rootT = sqrt(Temp)
        if (b(n).LE.-incl) then
          frac(n) = 0.
          z(n) = 0.
          go to 600
        else if (b(n).GT.incl) then
          frac(n) = sb(n)*cos(incl)
        else 
          x1 = cos(incl)*sb(n)*(acos(-tan(b(n))*(1./tan(incl))))/pi
          x2 = sin(incl)*cos(b(n))*sin(acos(-tan(b(n))/tan(incl)))/pi
          frac(n) = x1 + x2
        endif
        call sublime
        SUN=F0*FRAC(N)*(1.-A0)/RH**2
        RADIAT=(1.-A1)*SIGMA*Temp**4
        evap = q*root/rootT*press*xlt
        phi = radiat+evap-sun
        Z(n) = evap/xlt
        if (Z(n).lt.1.e-30)  Z(n) = 1.e-30
        drad = 4.*radiat/temp
        x1 = pprim*xlt
        x2 = p*xltprim
        devap = q*root/roott*(x1+x2)
        PHIPRI = drad + devap
        dt = sign(min(10.,abs(phi/phipri/2.)),phi/phipri)
        tp = temp - dt
        temp = tp
        if (abs(phi/sun).lt.1.e-4.or.abs(phi).lt.1.e-4) go to 300       
c  Added a flag to the script that prevents the server hanging if the code gets stuck. If the script
c  tries to loop more than 10 million times, the script kills itself.
        if (nflag.EQ.10000000) then
          subdis = Z(gd)*4*pi*rh*149.6*1E11
          if (perc.NE.0) then
            if ((subdis/perc).LT.1E-02) then
              Z(n)=0
              go to 300
            else
              print *,"Error calculating sublimation."
              stop
            endif
          endif
        endif
        nflag=nflag+1
        go to 100
  300   continue
9110    format(' Lat = ', f5.1, '  effcos = ', f6.3,
     &    '  Temp = ', f7.2, '  Z = ',1PE10.3)
        write (6,*) 'Incl: ',inc,'  Lat: ',b(n)
        gd=n
        perc=perc+Z(n)*4*pi*rh*149.6*1E11
  600 continue
      Zbar = 0.
      do 700 nn = 1, nb-1
        Zbar = Zbar + .5*(Z(nn)+Z(nn+1))*delsb
  700 continue
      Zbar=Zbar/2.
      Zlog = alog10(Zbar)
      rlog = alog10(rh)
      write (6,9112) spec, inc, rh, rlog, A0, A1, Zbar, Zlog
9112  format (1x,a8,2x,f4.1,2x,f5.2,1x,f6.3,1x,2(1x,f4.2),1x,
     &  1pe10.3,2x,0pf7.3)
c  Manually set a0 to -1. Originally fastrot went through a file with parameters until it reached a0=-1,
c  this ensures that it only runs the one time.
      a0 = -1
      go to 10
  800 continue
      END
