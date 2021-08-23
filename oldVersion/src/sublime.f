      subroutine sublime
C
c     Modification History
c       Brian Prager 06/25/2010
c         - Commented out the C02 pressures derived by G.N. Brown Jr. and W.T. Ziegler
c           in "Vapor Pressure and the Heats of Vaporization and Sublimation of Liquids
c           and Solids of interest in cryogenics below 1-atm pressure" (Advances in
c           Cryogenic Engineering, 25, 1980) [Also used to derive C0]. The script could
c           end up looking at temperatures below the minimum accepted value for their
c           derivation, and would crash the script. The pressure was then replaced with
c           the Antoine Equation, which is justified in another paper. (Given in comments
c           above the equation).
c
c
C     calculates the latent heat of sublimation and the vapor pressure of
C     the solid for various ices and the derivatives thereof
C     All arguments are transferred in the comon block /subl/.
C     Index specifies the ice: 1=h2o, 2=h2o clathrate, 3=co2, 4=co
C
      character*8 spec
      double precision press, pprim
      real mass
      common /subl/ press,pprim,index,temp,spec,mass,xlt,xltprim
c
      data sigma/5.67e-5/, f0/1.39e6/, dynmm/1.333e3/, boltz/1.38e-16/,
     :   ergcal/6.953e-17/, proton/1.67e-24/, aln10/2.30258/
      double precision p, pp
      character*8 species(4)
      data species/'H2O','H2O_CH4','CO2','CO'/
      real tstart(4)
      data tstart/190., 190., 100., 60./
C
c   initialize things
c   assume that, if subroutine is called with temp = 0, it is being done
c   to initialize parameters in the main program.  To speed convergence
c   in equilibrium vaporization calculations, a starting temperature 
c   must be picked that depends on the species under consideration.
c
c
      if (temp .LE. 0) then
        temp = tstart(index)
      endif
      aln10 = alog(10.)
      spec = species(index)
      t = temp
      t2 = t*t 
      t3 = t2*t
      t4 = t2*t2
      t5 = t4*t
      t6 = t3*t3

c
c  Now go to a specialized section that is different for each species.
c  The section first calculates the latent heat of vaporization in
c    calories per mole and its first derivative with respect to temperature.
c  It also calculates p, the log base 10 of the vapor pressure in mm 
c    of Hg, and the coefficient by which you must multiply the pressure 
c    to get its first derivative.  This is done because the vapor pressure 
c    is almost always an exponential, sometimes multiplied by a power of
c    temperature, and the conversion back to real units is done later.
c
      go to (100, 200, 300, 400) index
C
C     This is for ordinary water ice.  
C
  100 mass = 18.
      xlt = 12420. - 4.8*T
      xltprim = -4.8
c  This pressure formulation is from Marti & Mauersberger (1993 GRL 20, 363)
c  Note the factor 10 to convert from Pascals to dyn/cm^2.
      p = -2663.5/t + 12.537
      p = 10.*10.**p
      pp = (+2663.5/t2)*p
c     p = 10.**(1.20514d-5*T2-.01677006d0*T-2445.5646d0/T-6.757169d0)
c     p = dlog(p*T**8.2312)
c     pp = (8.2312/T+2.41028d-5*T-1.677006d-2+2445.5646d0/T**2)
c     p = 1.20514d-5*T2-.01677006d0*T-2445.5646d0/T-6.757169d0
c     p = p + 8.2312*alog10(t)
c     pp = 8.2312/T+aln10*(2.41028d-5*T-.01677006d0+2445.5646d0/T2)
      go to  5100
C
C   This is for clathrate hydrate of methane
c    We use the same formula for vapor pressure as for ordinary water
c    but a modified latent heat.
C
  200 mass = 18.
      xlt = 12160. + .5*T - .033*T2
      xltprim = 0.5 - 0.066*T
c  This pressure formulation is from Marti & Mauersberger (1993 GRL 20, 363)
c  Note the factor 10 to convert from Pascals to dyn/cm^2.
      p = -2663.5/t + 12.537
      p = 10.*10.**p
      pp = (+2663.5/t2)*p
c     p = 1.20514d-5*T2-.01677006d0*T-2445.5646d0/T-6.757169d0
c     p = p + 8.2312*alog10(t)
c     pp = 8.2312/T+aln10*(2.41028d-5*T-.01677006d0+2445.5646d0/T2)
      go to 5100
C
C  This is for dry ice, CO2,
C    from Brown & Ziegler for 194.694>T>40.0K
c  Sections commented out are for an older version of the
c    of the relevant quantities (used in our earlier work).
C
  300 mass = 44.
c     xlt = 6065.
c     xltprim = 0.
      xlt = 6269. + 9.877*t - .130997*t2 + 6.2735e-4*t3 - 1.2699e-6*t4
      xltprim = 9.877 -.261994*t + 1.88205e-3*t2 -5.0796e-6*t3
c     IF (T.GT.138.) then                                               CO2*****
c        P = (9.9082 - 1367.3/T)
c        PP = 1367.3/T2
c     else IF (T.LE.138.)  then                                         CO2*****
c        P = (8.307 + .00683*T - 1275.6/T)
c        PP = (1275.6/T2 + .00683)
c     endif

cc As of 06/23/2010 this portion was commented out due to errors in calculating pressures when
cc temperature reaches ~17 kelvin.
      p = 21.3807649d0 - 2570.647d0/t - 7.78129489d4/t2 + 4.32506256d6/
     :   t3 - 1.20671368d8/t4 + 1.34966306d9/t5
      pp = 2570.647d0/t2 + 1.556258978d5/t3 - 12.97518768d6/t4 +
     :  4.82685472d8/t5 - 6.7483153d9/t6
      p = dynmm*10.**p
      pp = pp*p
      go to 5100

cc As of 06/23/2010 Vapor Pressure is found using the Antoine Equation, with NIST Parameters
cc    A = 6.80657, B = 1301.679, C = 3.494. Justification for use with low temperatures is
cc    taken from 'Low Temperature Data for Carbon Dioxide' by Mustapha Azreg-Ainou (2005).
cc    This equation is supposedly valid even down to 0 kelvins, though the script will crash
cc    at around 4k.

C      if (t.lt.4.0) then
C        print *, 'error in CO2 temp, T less than 4k'
C        stop
C      else
C        p = 6.80657-(1301.679/(t-3.494))
C        pp = 1301.679/((t-3.494)**2)
C        p = dynmm*10.**p
C        pp = pp*p
C        go to 5100
C      endif
C
C   This is for solid CO, from Brown & Ziegler
C   beta form: for 68.127>T>61.544
C   alpha form: for 61.544>T>14.0K
C
  400 mass = 28.
      if (t.gt.68.127) then
        print *, 'error in CO temp, T =:', t
        stop
      elseif (t.gt.61.544) then
        xlt = 1855 + 3.253*t - .06833*t2 
        xltprim = 3.253 - .13666*t      
        p = 16.8655152d0 -748.151471d0/t - 5.84330795d0/t2 +
     :    3.93853859d0/t3
        pp = 748.15147d0/t2 + 11.6866159d0/t3 - 11.81561577d0/t4
        go to 5100
      elseif (t.ge.14.0) then
        xlt = 1893 + 7.331*t + .01096*t2 - .0060658*t3 +1.166e-4*t4
     :     - 7.8957e-7*t5
        xltprim = 7.331 + .02192*t - .0181974*t2 + 4.664e-4*t3
     :     - 3.94785e-6*t4
        p = 18.0741183d0 - 769.842078d0/t - 12148.7759d0/t2 +
     :    2.7350095d5/t3 - 2.9087467d6/t4 + 1.20319418d7/t5
        pp = 769.842078d0/t2 + 24297.5518/t3 - 820502.85d0/t4 +
     :    11634986.8d0/t5 - 60159709.d0/t6
      else
        print *, 'error in CO temp, T =:', t 
        stop
      endif
      p = dynmm*10.**p
      pp = pp*p
      go to 5100
C
 5100 mass = mass * proton
      xlt = xlt * ergcal
      xltprim = xltprim * ergcal
      press = p
      pprim = pp
C
      return
      end
