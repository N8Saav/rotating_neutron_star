! first line of file
!
! ----------------------------------------------------------------------------------
!  Modifications: 
!    07/15/93:  Parameter `ivonr' ><  computation of I(r)
!    08/02/93:  Spline, iv4
!    08/03/93:  761 engd(i) in sub _org_
!    08/02/93:  Changes in _spline_
!    08/18/93:  Changes in _swdf_
!    02/01/94:  Subroutine _resov_  ><  N(r), z(r) 
!    09/14/95:  Implemented do loop to compute stellar sequence with constant
!               baryon number, sub  _ovnrot_, _ht1000_ 
!    10/19/95:  Case dy=0 being treated separetely, sub _org_
!    05/22/96:  Density profiles of rotating stars, sub _IDEFST_
!    07/27/96:  Changed Do loop (see /07/27/1996)
!    08/27/96:  Tabulate stellar properties as done by WFF
!    10/17/96:  Tabulate mass quadrupole moment in km**3, sub _outrht_
!    12/13/96:  Mass of quark phase and mass of mixed phase as a function
!               of Omega, sub _ht1000_
!    06/05/05:  Write frame dragging frequency as a function of r to
!               output file omega_vs_r.out
!    06/13/05:  Add threshold densitites for CFL quark matter eos
!    06/29/06:  Introduce idefopt=1 or idefopt=2; amount of crust matter
!               is computed if idefopt=2
!    06/30/05:  Introduce common/ceinter/einter to compute stellar (crust)
!               mass below a density of einter=2.4e14 g/cm^3 
!    07/09/05:  Initialize iremov= 0, iremne= 0, iremht= 0, iov1=0
!    05/31/07:  Tabulate m(A), p(A), e(A), z(A), Phi(A) etc. as a function
!               of baryon number, A
!    09/20/12:  Made code Absoft fortran compiler compatible
!    10/16/14:  Fixed incorrect er_profile output. To compute e(r) 
!               set ioutse \= 0 in prns.dat. (The moment of inertial of 
!               the crust will not be computed in this case.)
!    04/09/15:  Changed all .out extensions to .dat
!    04/25/15:  Removed yes/no option from main program
!    04/25/15:  Values in rlomg no longer have to be monotonicaly increasing
!    04/30/15:  Fixed SIG 11 error (invalid memory reference) in HT1000 
!               (encountered during rotating white dwarf calculations)
!    05/08/15:  Removed unnecessary Do loop in _SEOS1_
!                                       DO 1233  I=1,L200
!                                 1233  AUX41(I,1)=AXNR1(I)  
!                   (wrong dimension -> AXNR1(52), L200=200)
!    05/08/15:  Introduced  iwhiteDwarf  as formal parameter
!    05/11/15:  Introduced subroutine SEOSWhiteDwarf (-> reads white dwarf
!               eos from external file (prns_eos.dat). Must set KST1=200 in prns.dat
!    05/13/15:  Compute internal energy if white dwarf eos is read from
!               external data file (prns_eos.dat)
!    05/03/15:  Change in _ht1000_ (look for 05/03/2015)
!    05/09/15:  Change in _INTER_ (look for 05/09/2015)
!    06/04/15:  Added additional output table (-> Compute Kepler frequency using 
!               Newton's expression, evaluated at the equator of rotating star)
!    06/05/15:  Relative fractions of quark matter and quark-hadron matter is
!               only calculated for neutron stars, but skipped for white dwarf   
!               calculations  
!    06/12/15:  Changed running index NAB 
!               Fixed padding of 4 bytes required before ee14 in COMMON/PowersTen/
!
! Version v10.4
!    06/18/15:  iwhiteDwarf and latent_heat are now read from prns.dat 
!               Issue with ioutse (particle thresholds) has been fixed       
!    06/19/15:  Several common blocks have been renamed (given more reasonable
!               names)
!    06/21/15:  Added one line which fixes rst=0 issue in _CRUST_
!
! Version v10.5
!    06/21/15:  Changed value of  npomg  from 40 to 50 -> all prns.dat files are
!               to be updated (iog40=50 now)
!    06/23/15:  Changed ALN (to speed up convergency for rotating white dwarf 
!               calculations
!    02/25/15:  Added common block /tabquad/tquad(. , .), which contains the mass
!               quadrupole moment of stars in km^3
!    03/05/16:  Added J/G M_sun^2 to output
!    06/04/16:  Introduced igravBary (=1,2,3)
!               1 = compute stellar model with a given gravitational mass
!               2 = compute stellar model with a given baryon mass          
!               3 = compute stellar model with a given baryon number
!      
!    09/19/17:  Program renamed to prns
!               'rotqs_ov.dat' changed to 'outputTOV.dat'      
!               'rotqs_ht.dat' changed to 'outputHT.dat'      
! Version v11.1
!    Introduced parameter for white dwarf calculations
!     iwhiteDwarf =   0  ><  Perform Neutron Star or Strange Star calculation       
!     iwhiteDwarf = -10  ><  Perform Strange Dwarf calculation
!     iwhiteDwarf =  10  ><  Perform White Dwarf calculation      
!
! Version v11.2
!    07/06/18: Did some clean-up
!    07/17/18: Show elapsed cpu_time; determine and print Date and Time 


! Version v12.1
!   08/09/2018 -- introduce threshold values for DD2
!   08/15/2018 -- output values of etasc6, etasc7, igivmx, iscmx7
!   08/16/2018 -- introduce common/BlobsRodsSlabs/thresh_brs(15,6), contains
!   threshold values of q-blobs, q-rods, q-slabs/h-slabs, h-rods, h-blobs, pure quark
!                                             
! Version v12.2
!   08/21/2018 -- value of ipzero (=20) changed to 40
!   04/11/2018 -- calculation for given e_c and freqency: use igravBary=0 in prns.dat

! Version v12.3
!   07/15/2020 -- typeos: changed character * 40 to character * 50
!   11/07/2020 -- density profiles will be written to mass_energy_density_profile.dat
!                 in subroutine profiles
!
! Version v13.0
!     02/16/2021 -- compiled code as
!     gfortran -fcheck=all, -ffpe-trap=invalid,zero,overflow -g -Wall  -pedantic 2> out.txt 
!     and removed all Warnings generated by the compiler; changes made in the code can be
!     found by searching for !*** symbol combination.
!
!     Note: 1> output.txt redirects stdout to a file 
!           2> error.txt redirects stderr to a file       
!----------------------------------------------------------------------------------
!
!  read input parameters from   prns.dat
!  read equation of state from  prns_eos.dat
!
!     gfortran -Wall -Wextra -Wconversion -o prns12.1.o prns12.1.f

!  compiler options (https://gcc.gnu.org/onlinedocs/gfortran/Invoking-GNU-Fortran.html):
!     -Wall: warn about all
!     -Wextra: warns about subroutine arguments that are never used
!     -Wconversion: warns about implicit conversion
!     -g: generate extra on debugging information usable by GDB
!     -g3: includes even more debugging information
!
!     gfortran -O -fno-align-commons -o prnsxx.f
!     gfortran -O -fbounds-check prnsxx.f
!     gfortran -std=legacy  prnsxx.f    ! turns of all warnings
!
! ------------------------------------------------------------------------
!
!
      program prns 

! ----------------------------------------------------------------------
!
! purpose:
!
!     formal main program for calculating the properties of (rotating)
!     neutron, hybrid, and quark matter stars with and without nuclear
!     crusts. The latter cover compact strange stars and strange white
!     dwarfs. 
!
!     parameters, arrays, constants:
!
!     2-dimensional grids:
!     pressure        p=p(j,i);  1 <= j <= nj,  1 <= i <= ni,
!     energy density  e=e(j,i);
!     radial distance r=r(j,i);
!
!     input for solving the stellar structure equations: equation of
!     state, P(e;rho), in units of  [p]=[e]=MeV/fm^3; [rho]=1/fm^3
!
!
!  P(e;rho)
!     |
!     |
!     |
!    \ /
!     |
!     |                /---- OV equations
!     |               /
!     ---------------/------ HT equations
!                    \
!                     \
!                      \---- classical Newtonian equations
!
!
!    optional:  1) stars rotating at their Kepler frequency
!               2) stars with given mass and rotational frequency 
!                  (smaller than Kepler frequency)
!               3) star with given baryon number and frequency
!               4) star with given central mass density and frequency
!
! ----------------------------------------------------------------------



      parameter (npr=600,npr1=npr+1,npt=80,npt1=npt+1,np200=600)
      parameter (np100=100,npvar=3*np100)

      real r(0:npr,0:npt),t(0:npt),p(0:npr,0:npt),e(0:npr,0:npt)
      real rth(0:npt,4),R_90(0:npr),R_00(0:npr),e_90(0:npr),e_00(0:npr)
      real ax1(0:npr),ax2(0:npr1),ax3(npr1),axx(np200),ayy(np200)

      real start, finish

      character (8) date
      character (10) time 

      call cpu_time(start)
      
      
 124  continue
      ni = npt
      nj = npr
      njp1 = npr1
      nip1 = npt1
      n200 = np200

!__ call organizing routine
      call org(r,t,p,e,ni,nj,rth,ax1,ax2,ax3,nip1,njp1,n200,axx,ayy,
     +         R_90,e_90,R_00,e_00)


      close ( 17 )

! determine and print CPU time      
      call cpu_time(finish)
      print '(" Elapsed CPU Time: ", F10.3, " seconds")', finish-start

! determine and print Date and Time
      call date_and_time(date, time)
      print '(" Date: ", A, ",", 1X, "Time: ", A, /)', date, time
      
      stop ' _prns_ --> regular stop encountered'

      end program prns
!org
      subroutine org(r,t,p,e,ni,nj,rth,ax1,ax2,ax3,nip1,njp1,n200,axx,ay
     +               y,R_90,e_90,R_00,e_00)


! ----------------------------------------------------------------------
!
! purpose:
!          o r g   calls the subroutines needed for the calculation
!          of general relativistic, non-rotating as well as rotating
!          neutron star models.
!
!  _ input
!  _ calculation of the grid points
!  _ preliminary calculations (starting values)
!
!  _ calculation of the different equations of state
!
!  _ computation of the parameters of  n o n -  rotating stars
!     (i.e., oppenheimer-volkoff treatment)
!
!  _ computation of the properties of  r o t a t i n g  stars in the
!      - classical newtonian neutron star models
!      - hartle-thorne treatment (GR models)
!
! ----------------------------------------------------------------------


      real * 8 rmsun,grav,gravk,gravkm,ee18,gcm2km,eem18
      real * 8 ee14,ee34,ee44
      real * 8 ee03,ee13,ee15
      
      real r(0:nj,0:ni),t(0:ni),p(0:nj,0:ni),e(0:nj,0:ni)
      real rth(0:ni,4),ax1(0:nj),ax2(0:njp1),ax3(njp1)
      real axx(n200),ayy(n200)
      real R_90(0:nj),e_90(0:nj),R_00(0:nj),e_00(0:nj)

      parameter (ipzero=40)
      real eczero(ipzero),yzero(ipzero)

      parameter (np100=100,npvar=3*np100,n27=50,n17=50,np4=26,npapp=10)
      real sc(4,npvar),zehw(n27)
      real * 8  c8(4,npvar)


      common/zzzz/zedy(n27),zmns(n27),zrns(n27),ztm(n27),zkepf(n27),
     +            zfrp(n27),fixmov(n27)
      common/zzz1/zedyov(n27),zedynw(n27),zedyht(n27)
      common/tabht1/tht1(n27,n17),tht2(n27,np4),tappx(n27,npapp)
      parameter (nptb3=10)
      common/tabbbb/tabb(n27,nptb3)
      common/ep1/engd(npvar),press(npvar),eintnl(npvar),einth(npvar),
     +           pden(npvar)
      common/ax1/aux1(npvar),aux3(npvar),aux4(npvar),aux31(1,npvar),
     +           aux41(1,npvar)
      common/ax2/aux100(npvar)

      common/MiscInt1/i100,i7,i27,ivar,inpog1,ikp1,icrust
      common/MiscInt2/j651,j650,igivmx,ivonr,kc20
      common/PowersTen/ee03,ee15,ee18,eem18,ee13,ee14,ee34,ee44
      common/MultiPi/pi,piq,zpi,vpi,zpi3,pid2
      common/SolarProperties/rmsun,rsun,rsunkm
      common/Gravitational/gravk,grav,gravkm,gcm2km
      common/ConversionFactors/uf1,uf12,uf13,uf4,uf6,uf10
      common/ConversionFactor30/ufe30
      common/MiscParameters/lim007,icro22,iine9,dzer9l,dzer9p
      common/RadialScalingFactor/r0r0
      common/ctype/typeos,teoscb
      common/NuclearMatterDensity/engnm0
      common/NucleonMasses/amnuc,amneut
      common/mishra/im05
      common/ConvergencyLimits/etasc6,etasc7,iscmx7
      common/iplot1/iradot,iout20,iotomg,iout33,ioutes,iotmp0,iop30,ioth
     +              v2,ioutpe,ioutse

      parameter (m651=160501)
      common/cnsta1/f251(m651),fout(m651),fst(m651),fst251(m651),feng1
     +              (m651),x251(m651),y251(m651),f251pm(m651),fdrip(m6
     +              51)
      common/cnsta2/aux51(1,m651),aux5(m651)
      common/cnsta3/aux61(1,m651),aux6(m651),dpdrax(m651)
      common/cnsta4/aux71(1,m651),aux7(m651),aux8(m651),aux9(m651)
      common/bn01/axbn1(m651),axbn2(m651),pdenex(m651),einthx(m651),
     +            deindp(m651),dehdrh(m651),pdent(m651)

      parameter (npomg=50,np10og=13)
      common/tabog1/rlomg(npomg),tfreq(npomg,np10og)

      parameter (nttby2=15)
      common/tbynb1/ttby(n27,nttby2)

      parameter (nthov1=12)
      common/tabsph/thov1(n27,nthov1)

      parameter (kpc20=20,kpc5=6)
      common/tc/tcore(kpc20,kpc5),r_core(kpc20)

      character * 30 argint
      character * 30 eass
      character * 19 selec0
      character * 15 selec1
      character * 25 xipov,xipht
      character * 50 typeos
      character * 45 teoscb
      character * 38 prkov,prkht
      character * 58 version
      character * 30 cthr

!__ values for eps(c) for hw test calculations  (in  g/cm^3)
      data (zehw(i),i=1,27)/ 1.0e03,1.e04,1.e05,1.e06,1.e07,1.e08,
     +                       2.5e08,1.e09,1.e10,1.e11,1.e12,1.e13,
     +                       1.5e13,1.7e13,2.e13,2.5e13,5.e13,1.e14,
     +                       3.e14,1.e15,3.e15,6.e15,1.e16,3.e16,
     +                       1.e17,3.e17,1.e18/


      prkov='Integration of OV eqs: not initialized'
      prkht='Integration of HRL eqs: not initialized'
!             '                                                        '
      version='program last modified on 30 June 2021 (prns13.0.f)      '

      it     = 0
      jsurf  = nj
      nxy    = 120  ! Note: nxy <= n200
      im05   = -10
      iremov = 0
      iremne = 0
      iremht = 0
      iov1   = 0



      call input(rmax,itmax,epsrad,epsmas,bagct,iqm,kst1,ldecov,
     +           ldecht,i4006a,i4006b,icmoi,selec0,selec1,edrip,
     +           pdrip,ecqmax,rcqmax,eqdrip,bord0,bord1,ejgcm3,x
     +           ipov,xipht,r5km,iq2020,iccaus,pcut,pcuts,ahw,bh
     +           w,anv,bnv,irukov,irukht,iccorr,isdwf,ico21a,ico
     +           21b,iarm,iwhiteDwarf,latent_heat,igravBary)
!     note: [bagct]=MeV,
!     return: axnr1=n_e,   [ .. ]=1/fm^3    \  use the HW+NV eos
!             axed1=e,     [ .. ]=MeV/fm^3  /  (low-density range)


      if (irukov == 0) then
      prkov='Integration of OV eqs via Euler method'
      else
      prkov='Integration of OV eqs via Runge-Kutta '
      end if
      if (irukht == 0) then
      prkht='Integration of HT eqs via Euler method'
      else
      prkht='Integration of HT eqs via Runge-Kutta '
      end if

      do 711 i = 1, i27
        zedyov(i) = zedy(i)
        zedynw(i) = zedy(i)
        zedyht(i) = zedy(i)
  711 continue

      i5651=1
      iipzo=ipzero
      if (igivmx > iipzo)   igivmx=iipzo
      if (zedy(1) == zedy(2)) i5651=igivmx



!__ compute equation of state

      if (kst1 == (-100)) then
!                *** *****
      call mishq(iep,icomb)

      else if (kst1 == (-90)) then
!                *** *****
      call can75(iep,icomb)

      else if (kst1 == (-80)) then
!                *** *****
      call ptrop(iep,icomb)

      else if ( (kst1 >= (-70)) .and. (kst1 <= (-50)) ) then

!                *** *****            /ep1/engd,press,pden,e_intnl
      call wff88(iep,icomb,sc,ivar,i100,kst1,iccorr)

      else if (kst1 == -40) then

!               *** *****            /ep1/engd,press,pden,e_intnl
      call fp81(iep,icomb,sc,ivar,i100,iccorr)

      else if (kst1 == -30) then

!                *** *****            /ep1/engd,press,pden,e_intnl
      call bps71(iep,icomb)

      else if (kst1 == -20) then

!               *** *****             /ep1/engd,press,pden,e_intnl
      call rd84(iep,icomb)

      else if (kst1 == -10) then
!                 *** *****           /ep1/engd,press,pden,e_intnl
      call hwht68(iep,icomb)

      else if (kst1 == 100) then

!                ***       *****      /ep1/engd,press,pden,e_intnl
      call seos1(iep,bagct,icomb,iqm,sc,ivar,rcqmax,ahw,bhw,anv,bnv)
!     return:  pden=n_e,       [ .. ]=1/fm^3    \
!                                                  over whole density range
!              eintnl=e_int,   [ .. ]=MeV/fm^3  /


      im05 = -10
      i200 = 200
      if (typeos=='Mishra et al.; T=0 MeV                            ')
     + im05=0
      if (typeos=='Mishra et al.; T=5 MeV                            ')
     + im05=5
      if (typeos=='Mishra et al.; T=0 hadronic + quark               ')
     + im05=10
      if (typeos=='E_SH93I                                           ')
     + im05=15

      if (im05 < 0) goto 5039

! access only if some pressure values are negative, i.e., P(e) < 0
      if (im05 == 0)  imi10 = 11
      if (im05 == 5)  imi10 = 10
      if (im05 == 10) imi10 = 11
      if (im05 == 15) imi10 = 33

      imi11 = imi10 + 1
      i210  = i200 + imi10
      i211  = i200 + imi11

! energy density
      a = alog10(engd(i210))
      b = alog10(engd(i211))
      argint='org: energy grid points Mishra'
!                        ****
      call intv(a,b,i210,aux7,0,argint)
      do 761 i=1,i210
  761 aux7(i) = 10**aux7(i)

! pressure
      aux6(1)    = press(1)
      aux6(i210) = press(i211)

! linear interpolation of P(e)
      xxm      = (aux6(i210)-aux6(1)) / (aux7(i210)-aux7(1))
      do 762  i=1,i210
  762 press(i) = xxm * (aux7(i)-aux7(1)) + press(1)

! baryon density
      a = pden(i210)
      b = pden(i211)*(1.-0.0001)
      argint='org: baryon density     Mishra'
!                        ****
      call intv(a,b,i210,aux5,0,argint)
      do 763 i=1,i210
  763 pden(i) = aux5(i)

! internal energy
      a =engd(i210) - amneut*pden(i210)
      b =engd(i211) - amneut*pden(i211)
      b = b*(1.-0.0001)
!     argint='org: internal energy    Mishra'
!                        ****
      call intv(a,b,i210,aux4,0,argint)
      do 764 i=1,i210
  764 eintnl(i) = aux4(i)

!      do 765 k=1,i200
!      engd(k)  = 0.
!      press(k) = 0.
!      pden(k)  = 0.
!  765 eintnl(k)= 0.

c  copy energy density
      do 766 i=1,i210
  766 engd(i) = aux7(i)


c    end of mishra // sharma
 5039 continue



      else if (kst1 == 110) then

!     eos of a quark star having a nuclear crust (BPS or HW)

!                ***       *****         ***** ***** ***** *****
      call seos2(iep,bagct,icomb,sc,ivar,edrip,pdrip,idria,idrib,
!                ******        ******
     +           itotal,ecqmax,eqdrip,ejgcm3,ahw,bhw,anv,bnv)
!     return:
!              /ep1/engd,press,pden,eintnl
!              pden=n_e,       [ .. ]=1/fm^3    \  
!                                                  whole density range
!              eintnl=e_int,   [ .. ]=MeV/fm^3  /

!     check if any of the e_c  values (input) lies in the range
!     e_drip  < e_c  < e_(d;quark)
        ichinp=0
      do i=1,i27
        if ( (zedy(i) < eqdrip) .and. (zedy(i) > edrip) ) then
        ichinp=-10
        zedy(i)  =edrip
        zedyov(i)=zedy(i)
        zedynw(i)=zedy(i)
        zedyht(i)=zedy(i)
        end if
      end do
      
        
      else if (kst1 == 200) then
!__ input eos of white dwarf matter from an external data file

      call SEOSWhiteDwarf(iep) 

         
      else

      write (6,'(///,'' :::: fatal error: kst1 not correctly initial'',
     +          ''ized;  kst1='',i4)') kst1
      stop

!__ determination of eos (determined by paramter kst1) ends here
      end if
!...........................................................................
      

      if (kst1 == -100) then
         teoscb='fmt+bps+bbp                                 '
         typeos='Misras quark eos: T=0 1_conserved_charge          '
      else if (kst1 == -90) then
         teoscb='fmt+bps+bbp                                 '
         typeos='CANUTOS favored eos: BJ_pure_neutron_gas          '
      else if (kst1 == -80) then
         teoscb='n=3/2 relativistic polytrope; [Tro65],[Har73]'
         typeos='relativistic polytrope                            '
      else if (kst1 == -70) then
         teoscb='fmt+bps+nv+wff88(av14+uvii)                  '
         typeos='WFF(AV14+UVII)                                    '
      else if (kst1 == -60) then
         teoscb='fmt+bps+nv+wff88(uv14+uvii)                  '
         typeos='WFF(UV14+UVII)                                    '
      else if (kst1 == -50) then
         teoscb='fmt+bps+nv+wff88(uv14+tni)                   '
         typeos='WFF(UV14+TNI)                                     '
      else if (kst1 == -40) then
         teoscb='fmt+bps+nv+fp81; cf. [FP81]                  '
         typeos='FP(81)                                            '
      else if (kst1 == -30) then
         teoscb='fmt+bps+bbp+pan; cf. [BPS71], [AB77]         '
         typeos='Pan(C)                                            '
      else if (kst1 == -20) then
         teoscb='fmt+bps+pps+bj(i); cf. [RD84], [AB77]        '
         typeos='BJ(I)                                             '
      else if (kst1 == -10) then
         teoscb='hw [HT68]                                    '
         typeos='HARRISON-WHEELER                                  '
      else if (kst1 == 100 .and. iqm <= 0) then
         teoscb='harrison-wheeler & negele-vautherin          '
      else if (kst1 == 100 .and. iqm == 10) then
         teoscb='pure quark matter eos (nothing to be joined) '
      else if (kst1 == 110 .and. icrust == 10) then
         teoscb='harrison-wheeler & negele-vautherin + quark  '
         typeos='QUARK MATTER (BAG)                                '
      else if (kst1 == 110 .and. icrust == 20) then
         teoscb='baym-pethick-sutherland + quark              '
         typeos='QUARK MATTER (BAG)                                '
      else if (kst1 == 200) then
         teoscb='White Dwarf eos (<-- prns_eos.dat)           '
! 07/15/20         typeos='None                                        '         
      end if

      if (im05 >= 0) 
     +   teoscb='star is self-bound!                          '
 

!     calculation of total baryon number of the non-rotating and
!     rotating neutron star:
      iatoty= 10
      if (pden(iep) == 0.) iatoty= - 10

!     check values of the internal energy (e_int > 0!)
      if377=10
      do 377 i=1,iep
      af=eintnl(i)
      if (af >= 0.) go to 377
      if377=-10
      write (6,'(/,'' :::: warning: internal energy < 0'',
     +           ''  e_int='',e14.6,''  i='',i3)') eintnl(i),i
  377 continue
! no access to line below
! if (if377 < 0) stop


      write (6,21)
   21 format (//' calculation of the eos successfully completed:'/)
      do 61 kk=1,iep,1

      if (iqm <= 0) then

      if (kst1 < 0) then
        if (kk == 1 .and. kst1 == -10) write (6,22)
        if (kk == 1 .and. kst1 == -20) write (6,26)
        if (kk == 1 .and. kst1 == -30) write (6,27)
        if (kk == 1 .and. kst1 == -40) write (6,28)
        if (kk == 1 .and. kst1 == -50) write (6,29)
        if (kk == 1 .and. kst1 == -60) write (6,30)
        if (kk == 1 .and. kst1 == -70) write (6,31)
        if (kk == 1 .and. kst1 == -80) write (6,32)
        if (kk == 1 .and. kst1 == -90) write (6,33)
        if (kk == 1 .and. kst1 == -100) write (6,34)
        else
        if (kk == 1 .and. iwhiteDwarf == 0) write (6,22)
        if (kk == 1 .and. iwhiteDwarf == 10) write (6,2211)
        if (kk == 1 .and. iwhiteDwarf == -10) write (6,2311)
        if (kk == (i100+1)) write (6,23)
        if (kk == (2*i100+1)) then
           write (6,24)
           write (6,244) typeos
           end if

      if (kk == icomb) write (6,25)

      end if
      
      else if (iqm == 10) then

      if (kk == 1) write (6,'(/,'' :::: bag constant  b**(1/4)='',f8.2,
     +   '' MeV'')') bagct
      if (kk == icomb) write (6,25)
      else if (iqm == 20) then
      if (kk == 1) write (6,'(/,'' :::: bag constant  b**(1/4)='',f8.2,
     +   '' MeV'',/,'' kst1='',i4,'' icrust='',i4,'' iqm='',i4,/,
     +   '' investigate crust thickness of a rotating neutron'',
     +   '' star'',/)') bagct,kst1,icrust,iqm
      if (kk == icomb) write (6,25)
      else
      write (6,'(///,'' :::: fatal input error; iqm out of range'',/,
     +           '' iqm='',i4,//)') iqm
      stop

      end if

 61   write (6,3) kk,engd(kk),press(kk),pden(kk),eintnl(kk)
   3  format (' k=',i3,2x,'e=',e12.5,2x,'p=',e12.5,2x,'n_e=',e12.5,
     +        2x,'e_i=',e12.5)
   22 format (/,23x,'< begin of  H W  eos >')
 2211 format (/,23x,'< begin of  White Dwarf  eos >')
 2311 format (/,23x,'< begin of  Strange White Dwarf  eos >')
 23   format (/,23x,'< begin of  N V  eos >')
   24 format (/,23x,'< begin of  H A D R O N I C  matter eos >')
  244 format (/,25x,a50,/)
   25 format (/,23x,'< begin of  Q U A R K  matter eos >')
   26 format (/,23x,'< begin of  BJ (I)  eos [rd84] >')
   27 format (/,23x,'< begin of  BPS (Pan (C))  eos [bps71] >')
   28 format (/,23x,'< begin of  FP81  eos: v14+TNI  [fp81] >')
   29 format (/,23x,'< begin of  WFF  eos: uv14+TNI  [wff88] >')
   30 format (/,23x,'< begin of  WFF  eos: uv14+uvii  [wff88] >')
   31 format (/,23x,'< begin of  Wff  eos: av14+uvii  [wff88] >')
   32 format (/,23x,'< n=3/2 RELATIV. POLYTROP. eos [tro65] >')
   33 format (/,23x,'< begin of BJ, neutron gas, Canuto 75  >')
   34 format (/,23x,'< BPS combined with Misras Quark eos   >')

      if (ioutse /= 0) then

      open(unit=91, file='seos.dat', status='unknown')

      do 6729 i=1,iep
      if (i == 1) write (91,'(/,''  e  and  P   (MeV/fm^3) '',/)')
 6729 write (91,6730) engd(i),press(i)
 6730 format (3x,e12.6,3x,e12.6)
      do 6792 i=1,iep
      if (i == 1)
     +write (91,'(/,''  rho (1/fm^3) ; e  and  P  (MeV/fm^3) '',/)')
 6792 write (91,6793) pden(i),engd(i),press(i)
 6793 format (3x,e12.6,3x,e12.6,3x,e12.6)
      end if

      close ( 91 )


!>>>>>  for hw test calculations only!
      ihwt6=10
      if (kst1 == -10 .and. ihwt6 /= 0) then
      write (6,'(/////,''>>>>>>>>>>>>>>>>>>>> hw test calculation '',
     +                 ''<<<<<<<<<<<<<<<<<<<<'',/)')
      do 667 i=1,i27
      zedy(i)=zehw(i)/uf6
      print*,' i=',i,' eps(c)=',zedy(i),' MeV/fm**3 ',zehw(i),' g/cm**3'
  667 continue
      end if
!>>>>>


      if (ichinp < 0)
     +write (6,'(///,''  NOTE: (some of the) input values of  e_c   '',
     +             ''have been set equal to'',/,8X,''e_c=e_join '')')


!__ check   c a u s a l i t y   of eos
      if (iccaus /= 0) call causal(engd,press,iep,sc,ivar)


      m11=0


!__ central energy density or gravitational mass loop
      do 4006 m = i4006a, i4006b


!__ write information about present status of the calculation to output file
!   prns_showStatus.dat
      if (isdwf == 0) call show(m,i4006a,i4006b,icab,ico21a,ico21b)


!__ initialize central energy density, [e_c]=MeV/fm^3
      ecov = zedyov(m)
      ecnw = zedynw(m)
      echt = zedyht(m)


!__ test segment for  e_c
      if ((ecov < engd(1)) .or. (ecnw < engd(1)) .or. (echt < engd(1)) )
     +   then
      write (6,'(///,'' ::::  e_center < eps(1)  '', / , ''Check in'',
     +          ''put file  prns_eos.dat <-- data may be wrong'', /,  
     +          '' m='',i3,/,'' e_c='',3e14.6,/,'' eps_1='',e14.6,/,
     +          '' continue calculation with next value of   e_c'',/
     +          /)') m, ecov, ecnw, echt, engd(1)
      end if



! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!             begin of spherical star calculations
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


      if (selec0 == 'oppenheimer-volkoff') then

!     iteration loop (calculation of e_c for a   g i v e n   NS mass)
      do 5651 ims=1,i5651


      if ( (i5651 /= 1) .and. (m == i4006a) .and. (ims == 1) ) then
      write (6,'(///)')
      write (6,712)
  712 format(1x,'----------------------------------------------',
     +                       '---------------------------------')
      write (6,'('' NOTE:'',/)')
      write (6,'('' Perform TOV calculation for a given mass or'',
     +           '' baryon number => e_center'')')
      write (6,712)
      write (6,'(///)')
      end if

      if (ims == 1) then
         if (i5651 == 1) then
            ecov=zedyov(m)
            eczero(ims)=zedyov(m)
            go to 7432
         end if
         ecov=zedyov(3)
         eczero(ims)=zedyov(3)
         else if (ims == 2) then
         ecov=zedyov(4)
         eczero(ims)=zedyov(4)
         else if (ims > 2) then
         ecov=ecnext
         eczero(ims)=ecnext
      end if

 7432 continue


! compute starting values for the star's radius, r(theta) by solving
! the ov equations for a non-rotating, spherically symmetric star
      if (isdwf == 0) then
      call ovnrot(engd,press,iep,ecov,aux3,aux4,ivar,ifail,rrns,rmqoz,
     +            wi,sc,aux31,aux41,rth,ni,m,iremov,xipov,edrip,eqdrip
     +            ,pdrip,idria,idrib,bord0,bord1,pden,ejgcm3,i4006a,i4
     +            006b,r5km,iq2020,sslt,ssgt,c8,pcut,pcuts,irukov,popm
     +            fm)

! write  | r | m(r) | epsilon(r) | n(r) | p(r)  to output file      
      call profiles(rrns,rmqoz,iremov,ecov,igravBary)
! mas been moved  

      else

      do 2503 icab=ico21a,ico21b
      bord0=r_core(icab)
         
! write information about present status of the calculation to output file
      call show(m,i4006a,i4006b,icab,ico21a,ico21b)
      call sdwf  (engd,press,iep,ecov,aux3,aux4,ivar,ifail,rrns,rmqoz,
     +            wi,sc,aux31,aux41,rth,ni,m,iremov,xipov,edrip,eqdrip
     +            ,pdrip,idria,idrib,bord0,bord1,pden,ejgcm3,i4006a,i4
     +            006b,r5km,iq2020,sslt,ssgt,c8,pcut,pcuts,irukov,popm
     +            fm,bagct,icab)
 2503 continue
      end if      

!     return: rrns  = radius (meter),
!             rmqoz = m_NS / m_sun,
!             wi    = log ( I (g cm^2) )


      if (ifail < 0) then
      write (6,'('' :::: _org_ failure to compute non-rotating star'',
     +           '' model (change input value of r0r0)'')')
      write (6,'('' e_c='',e13.6,3x,''rrns='',e13.6//)') ecov,rrns
!>>>>>
      if (m11 == 0) then
      write (6,'(//,'' In order to find a solution for the chosen'' ,
     +         '' e_c value, select'',/,'' a larger value for the '',
     +         ''parameter  r0r0 in input file  "prns.dat".'',/,
     +         '' The present value is r0r0='',e14.6,//)') r0r0
      m11=10
      end if
!>>>>>
      go to 5050
      end if


      if (i5651 == 1) go to 7433

      write(6,'(//,'' ==========================================='',
     +             ''======================'',/,''     Results of'',
     +              '' spherical mass determination:'')')
      write(6,'(   '' ==========================================='',
     +              ''======================'')')



      if (igravBary == 1 ) then 
         varMA=rmqoz      ! <-- compute stellar model with a given gravitational mass
         else if (igravBary == 2) then 
         varMA=ttby(m,8)  ! <-- compute stellar model with a given baryon mass          
         else if (igravBary == 3) then
         varMA=thov1(m,7) ! <-- compute stellar model with a given baryon number
         else
         print*, '::: Fatal input error:', 'igravBary=', igravBary,
     +           'must be 1, 2, or 3'
         call exit
      end if
      yzero(ims)=varMA-fixmov(m)

      if (ims <= 1) then
          go to 5651

          else if (ims == 2) then
          y1=yzero(ims-1)
          y2=yzero(ims)
          dy=y1-y2
          e1=eczero(ims-1)
          e2=eczero(ims)
          de=e1-e2
          if (dy /= 0.) then
             ecnext=e1 - de*y1/dy
             else 
             ecnext=(e1+e2)/2.
          end if
          else if (ims >= 3) then
               if (yzero(ims) < 0.) then
                   y1=yzero(ims)
                   y2=yzero(ims-1)
                   e1=eczero(ims)
                   e2=eczero(ims-1)
                   else
                   y1=yzero(ims-2)
                   y2=yzero(ims)
                   e1=eczero(ims-2)
                   e2=eczero(ims)
               end if
          de=e1-e2
          dy=y1-y2

      print*,' e1=',e1,' e2=',e2,' y1=',y1,' y2=',y2

          if (dy /= 0.) then
             ecnext=e1 - de*y1/dy
             else
             ecnext=(e1+e2)/2.
          end if
      end if

      diffec=(eczero(ims-1)-ecnext)/eczero(ims-1)
      diffec=abs(diffec)

         zedyov(m)=ecnext
         ecov     =ecnext

      print*,' ims=',ims,' M_(sph;ims)=',rmqoz,' e_c=',eczero(ims),
     +       ' diff=',diffec
      print*,' '

         if (diffec <= etasc6) go to 5652


 5651 continue
 5652          continue
 7433                   continue


!__ compute results                  return: /zzzz/
      call resov(m,rmqoz,rrns,wi,iremov,ecov,iarm)

      if ( (ldecov < 0) .and. (m >= (i4006a+1)) ) then
         if (zmns(m) < zmns(m-1)) goto 4477
      end if


!   save radius of the sperhical ov star
      do 10 i=0,ni
   10 rth(i,1)=rrns




!__ grid points (radial and theta direction)
      do 100 i=0,ni
      a = 0.
      b = rth(i,1)
      nab = njp1
      argint='org: radial grid points       '
!                       ***
      call intv(a,b,nab,ax2,0,argint)
      do 150 l=1,njp1
  150 ax1(l-1) = ax2(l)

! copy  ax1  to  r(j,.)
      do 200 j=0,nj
  200 r(j,i) = ax1(j)
  100 continue


! compute arc at the star's surface at  0<=theta<=#/2, phi=const

!               #/2
!   arc:  arc =  I d[theta] r(theta)
!                0

      a = 0.
      b = pid2
      nab = nip1
      argint='org: theta grid  points       '
!                       ***
      call intv(a,b,nab,ax2,0,argint)
      do 300 l=1,nip1
      ax1(l-1)=ax2(l)
  300 ax3(l) = r(jsurf,l-1)
      eass='org: call igrn 300            '
!                                        **
      call igrn(ax2,ax3,nip1,axx,ayy,nxy,u4,sc,ivar,10,2,aux51,eass)
!
! total circumference of the star for  phi=const
      uphic = 4.*u4


! compute arc at the star's surface at theta=#/2, 0<=phi<=2#

!               #/2
!   arc:  arc =  I d[phi] r(phi,theta=#/2)
!                0

      do 400 i=1,nip1
  400 ax3(i) = r(jsurf,ni)
      eass='org: call igrn 400            '
!                                      ****
      call igrn(ax2,ax3,nip1,axx,ayy,nxy,u4t0,sc,ivar,10,2,aux51,eass)

! total circumference at  theta=#/2
      utt0 = 4.*u4t0

! total circumference of a spherical star with radius r=r(theta=0)
      urth0 = zpi*r(jsurf,0)

! the circumference at theta=#/2 is trivially given by 2#r(theta=#/2)
      utt0a = zpi*r(jsurf,ni)



      ipt1=0
      if (ipt1 /= 0) then
      print*,'circumference at theta=#/2:'
      write (6,500) utt0,utt0a
  500 format (/' utt0=',e14.7,'(num) !=! utt0a=',e14.7,'(analyt)'//)

      print*,'circumference of a spherical star with radius r=r(theta=0)
     +:'
      write (6,510) urth0
  510 format (/' urth0=',e14.7//)

      yq = urth0/uphic
      print*,'circumference at phi=const:'
      write (6,520) uphic,yq
  520 format (/' uphic=',e14.7,3x,'u_spherical/u_deformed(theta=#/2)=',
     +        f8.5,//)
      end if

      end if


      

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!             begin of non-spherical star calculations
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 5050 continue


      write (6,'(//,18x,''begin of rotating NS calculation'')')
      write (6,'(18x,''treatment: '',a15)') selec1


! compute polar angle  theta lying in the range  [0 , #/2]
      a = 0.
      b = pid2
      nab = nip1
      argint='org: theta grid points        '
!                       ***
      call intv(a,b,nab,ax2,0,argint)
      do 307 l=1,nip1
  307 t(l-1)=ax2(l)


      if (selec1 == 'newtonian limit') then
! rotating neutron star treated in the newtonian limit

      omega=rlomg(1)*uf10
!                                                          ***** **
      call rnsnew(engd,press,iep,ecnw,aux3,aux4,ivar,ifail,rmqoz,wi,sc,
!                                       ***
     +            aux31,aux41,omega,p,r,rth,t,e,ni,nip1,nj,njp1,iov1,
     +            m,iremne,bord0,bord1)

      if (ifail < 0) then
      write (6,'('' _org_ failure to compute rotating star'',
     +           '' model (change input value of r0r0)'')')
      go to 4006
      end if

!__ compute results                              return: /tabn1/
      call resrn(m,rmqoz,rth,ni,wi,omega,ecnw)



      else if (selec1 == 'hartle-thorne  ') then

!   iteration loop (calculation of e_c for a   g i v e n   NS mass)
      do 5751 ims=1,i5651

      if ( (i5651 == 1) .and. (m == i4006a) .and. (ims == 1) ) then
      write (6,'(///)')
      write (6,712)
      if      (iscmx7 /= 1) then
      write (6,'(/,'' NOTE: perform self-consistent rotating'',
     +             '' NS/QS calculation for a given '',/,
     +             ''       central energy density value: =>'',
     +             '' O_K_GR(e_c;M_rot)'',/)')
      else if (iscmx7 == 1) then
      write (6,'(/,'' compute the properties of rotating neutron/'',
     +             ''quark stars for given values'',/,'' of '',
     +             ''e_c(i)  and  Omega(i) => M(e_c,Omega)'',/
     +           /,'' NOTE: not a self-consistent calculation '',/)')
      end if
      write (6,712)
      write (6,'(///)')
      end if

      if ( (i5651 /= 1) .and. (m == i4006a) .and. (ims == 1) ) then
      write (6,'(///)')
      write (6,712)
      write (6,'(''  NOTE:'',/)')
      write (6,'('' perform rotating NS/QS calculation for '',
     +           '' given values of mass or baryon number'',/,
     +           '' and rotational frequency => e_center'',/)')
      write (6,712)
      write (6,'(///)')
      end if

      if (ims == 1) then
         if (i5651 == 1) then
            echt=zedyht(m)
            eczero(ims)=zedyht(m)
            go to 7532
         end if
         echt=zedyht(3)
         eczero(ims)=zedyht(3)
         else if (ims == 2) then
         echt=zedyht(4)
         eczero(ims)=zedyht(4)
         else if (ims > 2) then
         echt=ecnext
         eczero(ims)=ecnext
      end if


 7532 continue
! solve the hartle-thorn equations for a rotating relativistic star
!                                                          **** ****
      call ht1000(engd,press,iep,echt,aux3,aux4,ivar,ifail,rsph,rqoz,
!                 **                ***    ******
     +            wi,sc,aux31,aux41,rth,ni,rqozsp,t,nip1,iqm,m,irem
     +            ht,sum124,ogc,xipht,eintnl,einth,pden,iatoty,icmo
!                         ******       ******
     +            i,sumw1,rmtoth,i5651,ogcoms,edrip,eqdrip,pdrip,id
     +            ria,idrib,r,nj,njp1,ax1,ax2,bord0,bord1,ejgcm3,i4
     +            006a,i4006b,r5km,iq2020,sslt,ssgt,pcut,pcuts,c8,i
     +            rukht,popmfm,R_90,e_90,R_00,e_00,cthr,anuefip,k_qmc,
     +            k_mix,iwhiteDwarf)

      tabb(m,1)=ogcoms
      if (ifail < 0) then
      write (6,'('' _org_ failure to compute rotating star'',
     +           '' model (change input value of r0r0)'')')
      go to 4006
      end if
      if (i5651 == 1) go to 7533

      write(6,'(//,'' ==========================================='',
     +             ''======================'',/,''     Results of'',
     +              '' rotational mass determination:'')')
      write(6,'(   '' ==========================================='',
     +              ''======================'')')

      
      if (igravBary == 1 ) then 
         varMA=rmtoth     ! <-- compute stellar model with a given gravitational mass
         else if (igravBary == 2) then 
         varMA=ttby(m,7)  ! <-- compute stellar model with a given baryon mass
      bymrlg=alog10(amneut) + ttby(m,4) - dlog10(rmsun)-alog10(ufe30)
      varMA=10.**bymrlg
      else if (igravBary == 3) then
         varMA=ttby(m,4)  ! <-- compute stellar model with a given baryon number   
         else
         print*,'::: Fatal input error:', 'igravBary=', igravBary,
     +           'must be 1, 2, or 3'
         call exit
      end if
      
      yzero(ims)=varMA-fixmov(m)

      
      if (ims <= 1) then
          go to 5751

          else if (ims == 2) then
               y1 = yzero(ims-1)
               y2 = yzero(ims)
               dy = y1-y2
               e1 = eczero(ims-1)
               e2 = eczero(ims)
               de = e1-e2
                  if (dy /= 0.) then
                     ecnext=e1 - de*y1/dy
                     else
                     ecnext=(e1+e2)/2.
                  end if
          else if (ims >= 3) then
                  if (yzero(ims) < 0.) then
                      y1 = yzero(ims)
                      y2 = yzero(ims-1)
                      e1 = eczero(ims)
                      e2 = eczero(ims-1)
                      else
                      y1 = yzero(ims-2)
                      y2 = yzero(ims)
                      e1 = eczero(ims-2)
                      e2 = eczero(ims)
                    end if
          de = e1-e2
          dy = y1-y2

          print*,' e1=',e1,' e2=',e2,' y1=',y1,' y2=',y2, 'varMA=',varMA

          if (dy /= 0.) then
             ecnext=e1 - de*y1/dy
             else
             ecnext=(e1+e2)/2.
          end if
      end if

      diffec     = (eczero(ims-1)-ecnext) / eczero(ims-1)
      diffec     = abs(diffec)
      tht2(m,13) = diffec

      zedyht(m) = ecnext
      echt      = ecnext

      print*,' ims=',ims,' M(rot;ims)=',rmtoth,' e_c=',eczero(ims),
     +       ' diff_ec=',diffec
      print*,' '
      if (diffec <= etasc6) go to 5752


 5751 continue
 5752          continue
 7533                   continue



!__ compute results                       return: /tabht1/
      call resrht(m,rth,ni,echt)

      if ( (ldecht < 0) .and. (m >= (i4006a+1)) ) then
         if (tht1(m,11) < tht1(m-1,11)) goto 4477
      end if


      else
      write (6,'(/,''  _org_ parameter  selec1='',a15)') selec1


      end if

      thov1(m,1)=tht1(m,1)
      thov1(m,2)=zedyov(m)/engnm0
      thov1(m,3)=zmns(m)
      thov1(m,4)=tht2(m,4)
      thov1(m,5)=ttby(m,2)
      thov1(m,6)=tht1(m,19)


 4006 continue



!__ write output to disk file(s)
 4477 continue

      call outov(nj,ni,iqm,iov1,kst1,selec0,iremov,bagct,xipov,i5651,
     +           bord0,bord1,ejgcm3,edrip,sslt,ssgt,r5km,iq2020,prkov
     +           ,popmfm,iccorr,isdwf,ico21a,ico21b,version,iwhiteDwarf,
     +           latent_heat,igravBary,ahw,bhw,anv,bnv)

      if (selec1 == 'newtonian limit') then
      call outrn(iov1,omega,iremne)

      else if (selec1 == 'hartle-thorne  ') then
      call outrht(iremht,sum124,bagct,ogc,xipht,icmoi,sumw1,i5651,iqm,
     +            bord0,bord1,ejgcm3,edrip,sslt,ssgt,r5km,iq2020,prkht
     +            ,popmfm,iccorr,kst1,echt,version,cthr,iwhiteDwarf,
     +            latent_heat,igravBary,i4006a,i4006b,ahw,bhw,anv,bnv)

      end if


!__ write output data needed for latent heat studies to disk file(s)
      if (latent_heat /= 0) then
      if (iwhiteDwarf == 0) then ! Not for white dwarfs
      RMoI = tappx(i4006b,9)
      call heating_project(iremht,echt,anuefip,k_qmc,k_mix,i4006a,RMoI)
      end if
              end if

! **********************************************************************
! write informations about:
!                         - accuracy of the calculated functions
!                         - smallest step size encountered by
!                           interpolating functions
!                         - number of necessary  e x t r a polations
!                           in subroutine  inter  and failing searches
!                           r zero points in subroutine  zero22
      write (6,'(///,'' infos:'',/)')
      write (6,777) iine9,icro22,dzer9l,dzer9p,sum124,iremov,iremne,irem
     +              ht,sumw1
      write (6,778) xipov,xipht
  777 format (15x,'iine9=',i6,/,15x,'icro22=',i6,/,15x,'dzer9l=',e14.6,/
     +      ,15x,'dzer9p=',e14.6,/,15x,'accuracy of function  omega_bar'
     +       ,'(r): ',e12.6,/,15x,'iremov=',i8,2x,'iremne=',i8,2x,'irem'
     +       ,'ht=',i8,/,15x,'accuracy of function  w(1;r): ',e12.6)
  778 format (15x,'interpolation in ovnrot: ',a25,/,
     +        15x,'interpolation in ht1000: ',a25,/)

! write input parameters:
      write (6,'(//,'' input parameters:'')')
      write (6,788) j651,nj,ni,iqm,iov1,kst1,selec1,iccorr,iwhiteDwarf,
     +              latent_heat
 788  format (19x,'j651=',i6,/,19x,'nj=',i4,';  ni=',i3,/,19x,'iqm=',
     +        i3,/,19x,'iov1=',i3,/,19x,'kst1=',i4,/,19x,'treatment: ',
     +        a13,/,19x,'iccorr=',i3,/,19x, 'iwhiteDwarf=',i3,/,19x, 
     +        'latent_heat=',i3, /)
! **********************************************************************



      if (ioutse /= 0) then

      open (unit=27, file='er_profile.dat', status='unknown')
      write (27,'(/,'' jobname = '',a50)') typeos
!      write (27,'('' Omega='',f8.2,'' 1/sec'')') tabb(i4006a,1)
!      write (27,'('' log_10 A_rot='',f9.5)') ttby(i4006a,4)

! equatorial radial distance and mass density profile: theta=pi/2
      write (27,'(/, ''  radial distance  | mass density in'',
     +               '' equatorial direction'')')
      write (27,'('' --------------------------------------------'')')      
      do 2324 j=0,nj,1
      R_90(j) = r(j,ni) * r0r0/ee18
      e_90(j) = aux7(j+1) * echt
      write (27,334) R_90(j),e_90(j)
 2324 continue

! polar radial distance and mass density profile: theta=0
      write (27,'(//, ''  radial distance  | mass density in polar'',
     +                '' direction'')')
      write (27,'('' --------------------------------------------'')')      
      do 2334 j=0,nj,1
      R_00(j) = r(j,0) * r0r0/ee18
      e_00(j) = aux7(j+1) * echt
      write (27,334) R_00(j),e_00(j)
 2334 continue
 334  format (3x,e14.8,4x,e14.8)
      stop '_IDEFST_ regular stop encountered'
      end if

             close ( 27 )

      return

      end subroutine org


!CAUSAL
      SUBROUTINE CAUSAL(ENGD,PRESS,IEP,SC,IVARDM)


! ----------------------------------------------------------------------
!     PURPOSE:
!              COMPUTE DERIVATIVE OF EOS,  d P / d e, AND CHECK IF IT
!              BECOMES ACAUSAL BEYOND A CERTAIL ENERGY DENSITY VALUE
! ----------------------------------------------------------------------


      REAL SC(4,IVARDM),ENGD(IVARDM),PRESS(IVARDM)


      PARAMETER (NP100=100,NPVAR=3*NP100)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)

      PARAMETER (m651=160501)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA3/AUX61(1,M651),AUX6(M651),DPDRAX(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)

      common/iplot1/iradot,iout20,iotomg,iout33,ioutes,iotmp0,iop30,ioth
     +              v2,ioutpe,ioutse

      COMMON/NUCLEARMATTERDENSITY/ENGNM0
      COMMON/MISCINT2/J651,J650,IGIVMX,IVONR,KC20

      CHARACTER * 27 KAUS
      CHARACTER * 30 ARGINT
      CHARACTER * 30 EASS
      CHARACTER * 32 RMKAUS

      DATA ELOWL/400./

      DO 1 I=1,IEP
      E=ENGD(I)
      ISAVE=I-1
      IF (E >= ELOWL) GO TO 2
    1 CONTINUE
      WRITE (6,'(///,'' CAUSALITY OF THE EOS CANNOT BE CHECKED'',/,
     +          '' CONTINUE CALCULATION !!!!!!!!!'',///)')
      RETURN
    2 CONTINUE

      DO 5 I=ISAVE+1,IEP
      IL=I-ISAVE
      AUX3(IL) =  ENGD(I)
      AUX4(IL) = PRESS(I)
    5 CONTINUE
      XEA=AUX3(1)
      XEB=AUX3(IL)

      IIT=J651
      ARGINT='CAUSAL: ENERGY DENSITY GRID   '
!                           ****      ENERGY DENSITY GRID  E(I)
      CALL INTV(XEA,XEB,IIT,AUX5,0,ARGINT)


      IPRT=0
      KAUS  =' CAUSALITY  N O T  VIOLATED'
      RMKAUS=' MATTER MICROSCOPICALLY STABLE  '


      JR35=-10
!  COMPUTE  DERIVATIVE   dP/de
      DO 35 J=1,IIT
      XNODE=AUX5(J)
      EASS='CAUSAL: CALL SELSPL 35        '
!                                    *****                RETURN: dP/de
      CALL SELSPL(AUX3,AUX4,IL,SC,10,DERIV,XNODE,20,EASS)
      AUX6(J)=DERIV
      IF (DERIV > 1.) THEN
      IF (JR35 < 0) JR35=J
      KAUS=' CAUSALITY  V I O L A T E D'
      END IF
      IF (DERIV < 0.) RMKAUS=' MATTER MICROSCOPICALLY UNSTABLE'
      IF (IPRT /= 0) WRITE (6,'('' I='',I6,''  e='',F7.2,''  dP/de='',
     +                      E10.4)') J,XNODE,DERIV
   35 CONTINUE

      WRITE (6,'(////,''  CAUSALITY CHECK ->'',A27,/)') KAUS
      WRITE (6,'(     ''  STABILITY CHECK ->'',A32,//)') RMKAUS
      IF ( KAUS == (' CAUSALITY  V I O L A T E D') ) THEN
      JC=0
      JCMX=50
      DO 45 J=JR35,IIT
      DV=AUX6(J)
      JC=JC+1
      IF (JC > JCMX) GO TO 55
      QZ=AUX5(J)/ENGNM0
      WRITE (6,'('' J='',I6,'' d P / d e='',E10.4,'' >1  FOR  e='',F7.1,
     +           '' MeV/fm^3;  e/e_0='',F5.2)') J,AUX6(J),AUX5(J),QZ
   45 CONTINUE
   55           CONTINUE
      END IF


      RETURN

      END


!CRUST
      SUBROUTINE CRUST(R,NI,NJ,NIP1,NJP1,THT,IREM,IDPREM,NHU,C,AX1,AX2,
!                          ** ****
     +                 RTH,WT,WT28,MLF,EC)


! ----------------------------------------------------------------------
!
!  PURPOSE: CALCULATION OF THE (HADRONIC) MASS CONTAINED IN THE CRUST
!           OF A ROTATING NEUTRON STAR.
!           NOTE: THE CRUST EXTENDS FROM THE NEUTRON DRIP DENSITY OF THE
!           ROTATING STAR OUTWARDS TO THE STAR'S SURFACE (P=0).
!
!  RETURN:  WT = MASS IN UNITS OF THE SOLAR MASS
!           WT28 = MASS IN UNITS OF 10**28 g
!
! ----------------------------------------------------------------------


      REAL * 8  DGRAV,A,AR0,RMSUN,GRAV,GRAVK,GRAVKM,EE18,EE55,GCM2KM,
     +          EEM18,EE14,EE34,EE44
      real * 8  ee03,ee13,ee15
      
      DIMENSION C(4,NHU),THT(0:NI),R(0:NJ,0:NI)
      DIMENSION AX1(0:NJ),AX2(0:NJP1),RTH(0:NI,4)


      PARAMETER (m651=160501)
      COMMON/CNSTA1/F251(M651),FOUT(M651),FST(M651),FST251(M651),
     +              FENG1(M651),X251(M651),Y251(M651),F251PM(M651),
     +              FDRIP(M651)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA3/AUX61(1,M651),AUX6(M651),DPDRAX(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)
      COMMON/CNSTA5/AUX10(M651),AUX11(M651),AUX12(M651),AUX13(M651)
      COMMON/CNSTA6/AUX11A(M651),AUX12A(M651),AUX13A(M651)
      COMMON/CNSTA7/AUX14(M651),AUX14A(M651)
      COMMON/BN01/AXBN1(M651),AXBN2(M651),PDENEX(M651),EINTHX(M651),
     +            DEINDP(M651),DEHDRH(M651),PDENT(M651)
      COMMON/HT68A/DNUEDR(M651),RLJ(M651),DOGDR(M651),DJ2DR(M651)
      COMMON/HT68B/DEDP(M651),RNUED(M651),OGDR(M651),DLNJ(M651)
      COMMON/HT68C/RM0D(M651),P0SD(M651),DM0DR(M651),DP0SDR(M651)
      COMMON/HT68D/DMSDR(M651),DPHDRH(M651),P2SD(M651),RM2H(M651)
      COMMON/HT68E/H0OUT(M651),H0IN(M651),DH0IN(M651)
      COMMON/HT68F/DV2DR(M651),DH2DR(M651),V2RIN(M651),H2RIN(M651)


      PARAMETER (N27=50,N17=50,NP4=26,NPAPP=10)
      COMMON/TABHT1/THT1(N27,N17),THT2(N27,NP4),TAPPX(N27,NPAPP)


      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      COMMON/Exponent55/EE55
      common/ConversionFactor30/ufe30
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/SOLARPROPERTIES/RMSUN,RSUN,RSUNKM
      COMMON/RADIALSCALINGFACTOR/R0R0
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/GRAVITATIONAL/GRAVK,GRAV,GRAVKM,GCM2KM


      CHARACTER * 30 ARGINT
      CHARACTER * 30 EASS


!>> CHECK DIMENSION
      IF (NHU < IVAR) THEN
      WRITE (6,'(///,'' FATAL ERROR IN _CRUST_; ERROR STOP'',//,
     +           '' NHU='',I4,3X,'' IVAR='',I4,///)') NHU,IVAR
      STOP 'ERROR STOP IN _CRUST_'
      END IF



      AR0   =R0R0
      A     =EC*AR0**3/(RMSUN*ufe30)
      DGRAV =(GRAV/ufe30)*(RMSUN*ufe30)/(UF1*AR0)

      IDPP1=IDPREM



!__ GRID POINTS (RADIAL AND THETA DIRECTIONS, r* and O*)
      NAB = NJP1
      DO 100 I=0,NI
!                                r*(drip)  and  r*(surface)
      AA = RTH(I,4)*EE18/R0R0
      BB = RTH(I,2)*EE18/R0R0


!>>> CHECK NUMERICAL ACCURACY OF THE INTEGRATION SCHEME
      ICHECK=0
      IF ( ICHECK /= 0 ) THEN
      WRITE (6,'(////,'' _CRUST_ CHECK NUMERICAL ACCURACY'',/,
     +                '' NOTE: LIMITS OF INTEGRATION ARE CHANGED'',/,
     +                '' M28_sph;d-s  = M28_rot;d-s  !!!! <--- see'',
     +                '' output file  outputHT.dat'',////)')
      AA = FST251(IDPP1)
      BB = FST251(IREM)
      END IF
!>>> END OF CHECKING


      ARGINT='CRUST:  [r*(drip),r*(surface)]'
!                         ***    GRID  [r*(drip),r*(surface)]
      CALL INTV(AA,BB,NAB,AX2,0,ARGINT)
      DO 150 L=1,NAB
      AX1(L-1) = AX2(L)
  150 CONTINUE


! COPY  AX1  TO  R(J,.)
      DO 200 J=0,NJ
!                           2-dim  r* - theta*  grid
  200 R(J,I) = AX1(J)
  100 CONTINUE



!__ r(s,d) - r(s,surface) GRID
!>>>      AA = THT2(MLF,21)*EE18/R0R0
!>>>      BB = THT1(MLF,5) *EE18/R0R0
      AA = FST251(IDPP1)
      BB = FST251(IREM)
      NAB = NJP1
      ARGINT='CRUST: [r*(s,d),r*(s,surface)]'
!                         ***    GRID  [r*(s,drip),r*(s,surface)]
      CALL INTV(AA,BB,NAB,AX2,0,ARGINT)
      DO 250 L=1,NAB
  250 AX1(L-1) = AX2(L)



! COPY ENERGY DENSITY TO AUXILIARY ARRAY (concerns region beyond quark matter core)
!  ->  x251(i) and y251(i) then run from i=1 to i's value at surface      
        IDPJJ = IDPREM
      DO 300 J=IDPP1,IREM
      X251(J- (IDPJJ-1)) = FST251(J)
      Y251(J- (IDPJJ-1)) = FENG1(J)
 300  CONTINUE


      IDEPS=IREM - (IDPJJ - 1)


      EASS='CRUST: CALL SELSPL 40012      '
 ! 06/12/2015: Subtraction of -1 to avoid error exit for rotating qs calculations
      DO 400 K=1,NAB - 1 
      XXX=AX2(K)
!                                      ****
      CALL SELSPL(X251,Y251,IDEPS,C,10,ERET,XXX,10,EASS)
      AUX7(K)=ERET
  400 CONTINUE




      WT=0.
      DO 1000 I=1,NI

         SUM   = 0.

         DELTH = THT(I)-THT(I-1)
         ASC   = THT(I)

         SINTH = SIN(ASC)
         DS    = DELTH*SINTH


         DO 2000 J=1,NJ

            DELR=R(J,I)-R(J-1,I)
            RST =R(J,I)
            RST2=RST*RST


!        RADIAL INTEGRATION OF CRUST MASS
            SUM  = SUM  + DELR*RST2*AUX7(J)

 2000       CONTINUE

         WT = WT + DS * SUM

 1000    CONTINUE

      WT           = VPI * A     * WT
      WT28         = WT  * RMSUN * EE55




      RETURN

      END


!IDEFST
      SUBROUTINE IDEFST(R,NI,NJ,NIP1,NJP1,THT,IREM,IDPREM,NHU,C,AX1,AX2,
     +                  RTH,MLF,EC,OGLTD,M5FP,R_90,e_90,R_00,e_00)


! ----------------------------------------------------------------------
!
!  PURPOSE: CALCULATION OF THE MOMENT OF INERTIA OF A ROTATIONALLY
!           DEFORMED STAR, WHICH CONSISTES OF A QUARK CORE AND A HADRONIC
!           MATTER CRUST (BELOW NEUTRON DRIP DENSITY).
!
!  ON RETURN: WIT  = TOTAL MOMENT OF INERTIA OF THE STAR, [WI]=g cm^2
!             WIC  = MOMENT OF INERTIA OF THE CRUST,  [WIC]=g cm^2
!
! ----------------------------------------------------------------------


      REAL * 8  DGRAV,A,AR0,RMSUN,GRAV,GRAVK,GRAVKM,EE18,EE55,GCM2KM,
     +          EEM18,EE14,EE34,EE44
      real * 8  ee03,ee13,ee15
      
      DIMENSION C(4,NHU),THT(0:NI),R(0:NJ,0:NI),R_90(0:NJ)
      DIMENSION R_00(0:NJ),e_00(0:NJ)
      DIMENSION AX1(0:NJ),AX2(0:NJP1),RTH(0:NI,4),e_90(0:NJ)


      PARAMETER (m651=160501)
      COMMON/CNSTA1/F251(M651),FOUT(M651),FST(M651),FST251(M651),
     +              FENG1(M651),X251(M651),Y251(M651),F251PM(M651),
     +              FDRIP(M651)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA3/AUX61(1,M651),AUX6(M651),DPDRAX(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)
      COMMON/CNSTA5/AUX10(M651),AUX11(M651),AUX12(M651),AUX13(M651)
      COMMON/CNSTA6/AUX11A(M651),AUX12A(M651),AUX13A(M651)
      COMMON/CNSTA7/AUX14(M651),AUX14A(M651)
      COMMON/BN01/AXBN1(M651),AXBN2(M651),PDENEX(M651),EINTHX(M651),
     +            DEINDP(M651),DEHDRH(M651),PDENT(M651)
      COMMON/HT68A/DNUEDR(M651),RLJ(M651),DOGDR(M651),DJ2DR(M651)
      COMMON/HT68B/DEDP(M651),RNUED(M651),OGDR(M651),DLNJ(M651)
      COMMON/HT68C/RM0D(M651),P0SD(M651),DM0DR(M651),DP0SDR(M651)
      COMMON/HT68D/DMSDR(M651),DPHDRH(M651),P2SD(M651),RM2H(M651)
      COMMON/HT68E/H0OUT(M651),H0IN(M651),DH0IN(M651)
      COMMON/HT68F/DV2DR(M651),DH2DR(M651),V2RIN(M651),H2RIN(M651)


      PARAMETER (N27=50,N17=50,NP4=26,NPAPP=10)
      COMMON/TABHT1/THT1(N27,N17),THT2(N27,NP4),TAPPX(N27,NPAPP)


      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      COMMON/Exponent55/EE55
      common/ConversionFactor30/ufe30
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/SOLARPROPERTIES/RMSUN,RSUN,RSUNKM
      COMMON/RADIALSCALINGFACTOR/R0R0
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/GRAVITATIONAL/GRAVK,GRAV,GRAVKM,GCM2KM

      CHARACTER * 30 ARGINT
      CHARACTER * 30 EASS


!>> CHECK DIMENSION
      IF (NHU < IVAR) THEN
      WRITE (6,'(///,'' FATAL ERROR IN _IDEFST_; ERROR STOP'',//,
     +           '' NHU='',I4,3X,'' IVAR='',I4,///)') NHU,IVAR
      STOP 'ERROR STOP IN _IDEFST_'
      END IF



      IFAIL1=10
      AR0   =R0R0
      A     =EC*AR0**3/(RMSUN*ufe30)
      DGRAV =(GRAV/ufe30)*(RMSUN*ufe30)/(UF1*AR0)


      DO 555 M5=1,M5FP

!  M5=1  ><  calculation of  I_total
!  M5=2  ><  calculation of  I_crust


      IDPP1=1
      IF (M5 == 2) IDPP1=IDPREM

      do 578 i=0,nj
             ax1(i)=0.
             ax2(i)=0.
 578         continue
             ax2(njp1)=0.


!__ GRID POINTS (RADIAL AND THETA DIRECTIONS, r* and O*)
      NAB = NJP1
      DO 100 I=0,NI
!                                r*(drip)  and  r*(surface)
                   AA = FST251(1)
      IF (M5 == 2) AA = RTH(I,4)*EE18/R0R0
                   BB = RTH(I,2)*EE18/R0R0




      ARGINT='IDEFST: [r*(drip),r*(surface)]'
!                         ***    GRID  [r*(drip),r*(surface)]
      CALL INTV(AA,BB,NAB,AX2,0,ARGINT)
      DO 150 L=1,NAB
      AX1(L-1) = AX2(L)
  150 CONTINUE


! COPY  AX1  TO  R(J,.)
      DO 200 J=0,NJ
!                           2-dim  r* - theta*  grid
  200 R(J,I) = AX1(J)
  100 CONTINUE



!__ r(s,d) - r(s,surface) GRID
      AA = FST251(IDPP1)
      BB = FST251(IREM)
      NAB = NJP1
      ARGINT='IDEFST:[r*(s,d),r*(s,surface)]'
!                         ***    GRID  [r*(s,drip),r*(s,surface)]
      CALL INTV(AA,BB,NAB,AX2,0,ARGINT)
      DO 250 L=1,NAB
  250 AX1(L-1) = AX2(L)



!  COPY ENERGY DENSITY, PRESSURE, FRAME DRAGGING FUNCTION  omega(r), ALL
!  ARE FUNCTIONS OF RADIAL DISTANCE, TO AUXILIARY ARRAYS
                   IDPJJ = 1
      IF (M5 == 2) IDPJJ = IDPREM
      DO 300 J=IDPP1,IREM
      X251(J- (IDPJJ-1)) = FST251(J)
      Y251(J- (IDPJJ-1)) = FENG1(J)
      AUX5(J- (IDPJJ-1)) = FOUT(J)
      AUX6(J- (IDPJJ-1)) = OGDR(J)
      AUX8(J- (IDPJJ-1)) = F251(J)
      AUX9(J- (IDPJJ-1)) = RNUED(J)
      AUX10(J-(IDPJJ-1)) = H0IN(J)
      AUX11(J-(IDPJJ-1)) = H2RIN(J)
      AUX12(J-(IDPJJ-1)) = V2RIN(J)
      AUX13(J-(IDPJJ-1)) = RM0D(J)
      AUX14(J-(IDPJJ-1)) = RM2H(J)
  300 CONTINUE

      IDEPS=IREM - (IDPJJ - 1)

! 06/12/2015: Subtraction of -1 to avoid error exit for rotating qs calculations      
      NABM1=NAB -1
      
      EASS='IDEFST: CALL SELSPL 400       '
      DO 400 K=1,NABM1
      XXX=AX2(K)
!                                      ****
      CALL SELSPL(X251,Y251,IDEPS,C,10,ERET,XXX,10,EASS)
      AUX7(K)=ERET
  400 CONTINUE



      EASS='IDEFST: CALL SELSPL 401       '
      DO 401 K=1,NABM1
      XXX=AX2(K)
!                                      ****
      CALL SELSPL(X251,AUX5,IDEPS,C,10,PRET,XXX,10,EASS)
      AXBN1(K)=PRET
  401 CONTINUE

      EASS='IDEFST: CALL SELSPL 402       '
      DO 402 K=1,NABM1
      XXX=AX2(K)
!                                      ****
      CALL SELSPL(X251,AUX6,IDEPS,C,10,DRAG,XXX,10,EASS)
      AXBN2(K)=DRAG
  402 CONTINUE

      EASS='IDEFST: CALL SELSPL 403       '
      DO 403 K=1,NABM1
      XXX=AX2(K)
!                                      ****
      CALL SELSPL(X251,AUX8,IDEPS,C,10,RMAS,XXX,10,EASS)
      AUX51(1,K)=RMAS
  403 CONTINUE

      EASS='IDEFST: CALL SELSPL 404       '
      DO 404 K=1,NABM1
      XXX=AX2(K)
!                                      ***
      CALL SELSPL(X251,AUX9,IDEPS,C,10,RNU,XXX,10,EASS)
      AUX61(1,K)=RNU
  404 CONTINUE

      EASS='IDEFST: CALL SELSPL 405       '
      DO 405 K=1,NABM1
      XXX=AX2(K)
!                                       *****
      CALL SELSPL(X251,AUX10,IDEPS,C,10,RH0IN,XXX,10,EASS)
      AUX71(1,K)=RH0IN
  405 CONTINUE

      EASS='IDEFST: CALL SELSPL 406       '
      DO 406 K=1,NABM1
      XXX=AX2(K)
!                                       ******
      CALL SELSPL(X251,AUX11,IDEPS,C,10,RH2RIN,XXX,10,EASS)
      AUX11A(K)=RH2RIN
  406 CONTINUE

      EASS='IDEFST: CALL SELSPL 407       '
      DO 407 K=1,NABM1
      XXX=AX2(K)
!                                       ******
      CALL SELSPL(X251,AUX12,IDEPS,C,10,RV2RIN,XXX,10,EASS)
      AUX12A(K)=RV2RIN
  407 CONTINUE

      EASS='IDEFST: CALL SELSPL 408       '
      DO 408 K=1,NABM1
      XXX=AX2(K)
!                                       *****
      CALL SELSPL(X251,AUX13,IDEPS,C,10,RRM0D,XXX,10,EASS)
      AUX13A(K)=RRM0D
  408 CONTINUE

      EASS='IDEFST: CALL SELSPL 409       '
      DO 409 K=1,NABM1
      XXX=AX2(K)
!                                       *****
      CALL SELSPL(X251,AUX14,IDEPS,C,10,RRM2H,XXX,10,EASS)
      AUX14A(K)=RRM2H
  409 CONTINUE




      IF (M5 == 1) WIT = 0.
      IF (M5 == 2) WIC = 0.

! polar angle (theta) loop
      DO 1000 I=1,NI

         SUMIT = 0.
         SUMIC = 0.

         RTHDP = RTH(I,4)*EE18/R0R0


         DELTH = THT(I)-THT(I-1)
         ASC   = THT(I)

         SINTH = SIN(ASC)
         SINTH2= SINTH*SINTH
         COSTH = COS(ASC)
         COSTH2= COSTH*COSTH
         DS    = DELTH*SINTH

         PL2   = (3.*COSTH2-1.) / 2.


! radial loop (theta=const)
         DO 2000 J=1,NJ

            DELR = R(J,I)-R(J-1,I)
            RST  = R(J,I)
       if(rst == 0.) rst=r(j-1,i)  ! 06/21/2015:  addes this line as rst=0 for j=NJ
            RST2 = RST*RST


!        METRIC FUNCTIONS:  exp(lambda), exp(nue), exp(mue), exp(psi)
!        PRELIMINARY CALCULATIONS:
            UPSIL=2.*AUX51(1,J)*DGRAV/RST
            ONEM = 1. - UPSIL
            IF (ONEM <= 0.) THEN
            IF (IFAIL1 > 0)
     +      WRITE (6,'(///,'' 1-UPSILON < 0; SET MOMENT OF INERTIA OF'',
     +              '' CRUST = 0 AND CONTINUE CALCULATION'',///)
     +             ')
            UPSIL = 0.
            ONEM  = 1. - UPSIL
            IF (IFAIL1 > 0) IFAIL1 = -10
            END IF


!        exp(lambda)
            ZM0M2= AUX13A(J)*DGRAV + AUX14A(J)*PL2
            EXPL = (1. + 2.*ZM0M2/(RST*ONEM)) / ONEM

!        exp(nue)
            EXPN = EXP(AUX61(1,J)) * (1.+2.*(AUX71(1,J)+AUX11A(J)*PL2))

!        exp(mue)/r(0)**2
            EXPMR = RST2 * (1. + 2.*(AUX12A(J) - AUX11A(J)) * PL2)


!        r(0)**2 * exp(nue-psi)
            R2EXP = EXPN / (SINTH2 * EXPMR)


            EXP3 = EXPL * EXPMR * EXPN
            SEXP3= SQRT(EXP3)
            RATOG= AXBN2(J) / OGLTD
            RATIO= ( (AUX7(J) + AXBN1(J)) * RATOG ) /
     +             ( R2EXP - AXBN2(J)*AXBN2(J) )



      IF (M5 == 1) THEN
!        RADIAL INTEGRATION
!           TOTAL MOMENT OF INERTIA
            SUMIT = SUMIT + DELR * RST * SEXP3 * RATIO

      ELSE
!           MOMENT OF INERTIA OF HADRONIC CRUST
            SUMIC = SUMIC + DELR * RST * SEXP3 * RATIO

      END IF


 2000       CONTINUE


         IF (M5 == 1) WIT = WIT + DS * SUMIT
         IF (M5 == 2) WIC = WIC + DS * SUMIC


 1000    CONTINUE
  555               CONTINUE

      IF (M5FP == 1) THEN
      WIT = VPI * WIT
      ELSE
      RIIC= WIC / WIT
      WIT = VPI * WIT
      WIC = VPI * WIC
            END IF

!  log I_crust(rotationally deformed star) / g cm**2
      IF (IFAIL1 < 0) THEN
      THT2(MLF,22) = -1.
      THT2(MLF,23) = -1.
      THT2(MLF,24) = -1.
      GO TO 12345
      ELSE
      IF (M5FP == 2) THEN
      RLOGIC = ALOG10(WIC) + ALOG10(EC*UF6) + 5.*DLOG10(R0R0/EE13)
      ELSE
      RLOGIC = 0.
               END IF
                       END IF

      THT2(MLF,22) = RLOGIC

!///////////////////////////// no access \\\\\\\\\\\\\\\\\\\\\\\\\\\\
!      zjqq  = fst251(idprem+1)**4 * dogdr(idprem+1)*rlj(idprem+1)
!     +        / (6.*dgrav*a)
!      rlzjqq= ALOG10(zjqq) + ALOG10(EC*UF6) + 5.*ALOG10(R0R0/EE13)
!     +       -alog10(ogltd)
!      print*,'  log I_quark=',rlzjqqc

!      zjtt  = fst251(irem)**4 * dogdr(irem)*rlj(irem)
!     +        / (6.*dgrav*a)
!      rlzjtt= ALOG10(zjtt) + ALOG10(EC*UF6) + 5.*ALOG10(R0R0/EE13)
!     +       -alog10(ogltd)
!      print*,'  log I_total=',rlzjtt,'  compare with table'
!//////////////////////// end of no access \\\\\\\\\\\\\\\\\\\\\\\\\\


!  log I_total(rotationally deformed star) / g cm**2
      RLOGIT = ALOG10(WIT) + ALOG10(EC*UF6) + 5.*DLOG10(R0R0/EE13)
      THT2(MLF,24) = RLOGIT


!  I_crust / I_total
      IF (M5FP == 2) THEN
      DIFFLG = RLOGIC - RLOGIT
      THT2(MLF,23) = 10.**DIFFLG
      ELSE
      THT2(MLF,23) = 0.
                        END IF



12345 RETURN

      END



!input
      subroutine input(rmax,itmax,epsrad,epsmas,bagct,iqm,kst1,ldecov,
     +                 ldecht,i4006a,i4006b,icmoi,selec0,selec1,edrip,
     +                 pdrip,ecqmax,rcqmax,eqdrip,bord0,bord1,ejgcm3,x
     +                 ipov,xipht,r5km,iq2020,iccaus,pcut,pcuts,ahw,bh
     +                 w,anv,bnv,irukov,irukht,iccorr,isdwf,ico21a,ico
     +                 21b,iarm,iwhiteDwarf,latent_heat,igravBary)

! ----------------------------------------------------------------------
!
!  purpose: input routine, calculation of constants
!
!  [gravk] = fm/MeV, [grav] = fm^2, [gravkm] = km/MeV, [uf1] = MeV fm,
!
!  [rmsun] = MeV, [rsun] = km, [rmax] = km
!
!  assign preliminary values to variables  pdrip, edrip, eqdrip
!
! ----------------------------------------------------------------------


      real * 8 rmsun,grav,gravk,gravkm,ee55,ee18,gcm2km,eem18,
     +         ee14,ee34,ee44
      real * 8 ee03,ee13,ee15

      parameter (np100=100,npvar=3*np100,n7=7,n27=50,kp1=1)
!__ npvar = 3*np100!
      common/zzzz/zedy(n27),zmns(n27),zrns(n27),ztm(n27),zkepf(n27),
     +            zfrp(n27),fixmov(n27)
      common/zzz1/zedyov(n27),zedynw(n27),zedyht(n27)
      common/ax1/aux1(npvar),aux3(npvar),aux4(npvar),aux31(1,npvar),
     +           aux41(1,npvar)


      parameter (m651=160501)
      common/SolarProperties/rmsun,rsun,rsunkm
      common/Gravitational/gravk,grav,gravkm,gcm2km
      common/PowersTen/ee03,ee15,ee18,eem18,ee13,ee14,ee34,ee44
      common/Exponent55/ee55
      common/ConversionFactor30/ufe30
      common/MultiPi/pi,piq,zpi,vpi,zpi3,pid2
      common/MiscInt1/i100,i7,i27,ivar,inpog1,ikp1,icrust
      common/MiscInt2/j651,j650,igivmx,ivonr,kc20
      common/RadialScalingFactor/r0r0
      common/ExpoMaxMin/exmax,exmin
      common/ConversionFactors/uf1,uf12,uf13,uf4,uf6,uf10
      common/NucleonMasses/amnuc,amneut
      common/bag123/bag1,bag2,bag3
      common/MiscParameters/lim007,icro22,iine9,dzer9l,dzer9p
      common/ConvergencyLimits/etasc6,etasc7,iscmx7
      common/NuclearMatterDensity/engnm0
      common/Thresholds/thresh(15,2),rmtotm,rmtotq,nthrcb
      common/BlobsRodsSlabs/thresh_brs(15,6)
      common/iplot1/iradot,iout20,iotomg,iout33,ioutes,iotmp0,iop30,ioth
     +              v2,ioutpe,ioutse
      common/iplot2/ioutd4,ipolxx,icmi99
      common/GivenInteriorDensity/einter
      
      parameter (npomg=50,np10og=13)
      common/tabog1/rlomg(npomg),tfreq(npomg,np10og)

      parameter (kp52=52)
      common/caxnr1/axnr1(kp52),axed1(kp52)

      parameter (kpc20=20,kpc5=6)
      common/tc/tcore(kpc20,kpc5),r_core(kpc20)


      
      character * 16 eodf
      character * 19 selec0
      character * 15 selec1
      character * 18 bodf
      character * 25 xipov,xipht
      character * 30 dd


      namelist/nainp1/i7,i100,ivar,j650,j651,ikp1,iq202l,r5kml


!     Harrison-Wheeler particle density for the low-density region
!     n_e  (in 1/fm^3) for 0.44e-11 < e/(MeV fm^-3) < 56
      data (axnr1(i),i=1,kp52)/0.4681805E-14,0.4697894E-14,
     +0.4713982E-14,0.4802470E-14,0.4858781E-14,0.5084022E-14,
     +0.6853776E-14,0.9894537E-14,0.1882375E-13,0.5260998E-13,
     +0.1254917E-12,0.2373080E-12,0.4199145E-12,0.6853775E-12,
     +0.1230783E-11,0.2984448E-11,0.5960847E-11,0.1496242E-10,
     +0.3764726E-10,0.5960801E-10,0.1190544E-09,0.2984360E-09,
     +0.5960566E-09,0.1190472E-08,0.2984060E-08,0.5959713E-08,
     +0.1190238E-07,0.2983219E-07,0.5957581E-07,0.1189701E-06,
     +0.2981417E-06,0.5953157E-06,0.1188624E-05,0.2977877E-05,
     +0.5944445E-05,0.9381951E-05,0.1874966E-04,0.3746651E-04,
     +0.5928276E-04,0.9353877E-04,0.1485950E-03,0.2354980E-03,
     +0.3733539E-03,0.5907762E-03,0.9322584E-03,0.1481124E-02,
     +0.2347579E-02,0.3721392E-02,0.5886203E-02,0.1173025E-01,
     +0.2925967E-01,0.5808013E-01/

!     energy density corresponding to n_e, ranging from
!     0.44e-11 < eps/(MeV fm**-3) < 56
      data (axed1(i),i=1,kp52)/0.4399024E-11,0.4414141E-11,
     +0.4429258E-11,0.4512401E-11,0.4565310E-11,0.4776947E-11,
     +0.6439808E-11,0.9296906E-11,0.1768680E-10,0.4943233E-10,
     +0.1179120E-09,0.2229746E-09,0.3945516E-09,0.6439808E-09,
     +0.1156444E-08,0.2804189E-08,0.5600819E-08,0.1405874E-07,
     +0.3537360E-07,0.5600819E-07,0.1118652E-06,0.2804189E-06,
     +0.5600819E-06,0.1118652E-05,0.2804189E-05,0.5600819E-05,
     +0.1118652E-04,0.2804189E-04,0.5600819E-04,0.1118652E-03,
     +0.2804189E-03,0.5600819E-03,0.1118652E-02,0.2804189E-02,
     +0.5600819E-02,0.8843399E-02,0.1768680E-01,0.3537359E-01,
     +0.5600819E-01,0.8843399E-01,0.1405874E+00,0.2229746E+00,
     +0.3537360E+00,0.5600819E+00,0.8843399E+00,0.1405874E+01,
     +0.2229746E+01,0.3537359E+01,0.5600819E+01,0.1118652E+02,
     +0.2804189E+02,0.5600819E+02/



      data drmax/15./,itmaxd/1/,depsra/1.e-05/,depsma/1.e-05/
!       data rmsun/1.115829e60/,rsun/6.96e05/,rsunkm/1.476136/
!       data grav/2.611984e-40/,gravk/1.32367e-42/,gravkm/1.32367e-60/
!       data gcm2km/7.421943e-44/
      data rmsun/1.115829e30/,rsun/6.96e05/,rsunkm/1.476136/
      data grav/2.611984e-10/,gravk/1.32367e-12/,gravkm/1.32367e-30/
      data gcm2km/7.421943e-14/
      data uf1/197.329/,uf4/1.6022e33/,uf6/1.782687e12/
      data exmax/200./,exmin/-200./
      data amnuc/939./,amneut/939.6/
      data engnm0/140./,ufe30/1.e30/
      data ee03/1.e03/,ee15/1.e15/,ee18/1.e18/,eem18/1.e-18/
!       data ee55/1.782531e-55/,ee13/1.e13/
!       data ee14/1.e14/,ee34/1.e34/,ee44/1.e44/
      data ee55/1.782531e-25/,ee13/1.e13/
      data ee14/1.e14/,ee34/1.e34/,ee44/1.e14/

! 06/13/2005 -- threshold values for ioutse <= 70
!                              start  end of mixed phase      
!                              |      |
      data (thresh(1,j),j=1,2)/252.41,1391.25/  ! G_240^B180       ioutse=10
      data (thresh(2,j),j=1,2)/235.20,951.88/   ! G_300^B180       ioutse=20
      data (thresh(3,j),j=1,2)/242.06,1008.87/  ! G_300^B180_m7_a0 ioutse=30
      data (thresh(4,j),j=1,2)/239.20,953.83/   ! G_350^B180_m73   ioutse=40
      data (thresh(5,j),j=1,2)/240.24,957.40/   ! G_290^B180_m7    ioutse=50
      data (thresh(6,j),j=1,2)/262.00,930.00/   ! Alford CFL1      ioutse=60
      data (thresh(7,j),j=1,2)/901.83,1750.00/  ! epNL_GM1_Gv009   ioutse=70
      data (thresh(8,j),j=1,2)/421.159,1263/    ! DD2_npemu_ms119_gv00_hybrid  ioutse=80

      
! 08/16/2018 -- q-blob, q-rod, q-slab threshold values

!                                  q-blobs q-rods q-slabs h-rods       
!                                  |       |      |       |
      data (thresh_brs(8,j),j=1,6)/421.159,633.34,765.94,1005.34,
     +                                                   1107.56,1263.0/  ! ioutse=80
!                                                        |       |
!                                                        h-blobs pure quark   
      
      data einter/2.4e14/       ! nuclear matter density
!      data einter/4.3e14/   ! neutron drip density
      
      xipov = 'has not been invoked'
      xipht = 'has not been invoked'
  
      iccorr = 0
      lim007 = -5
      icro22 = 0
      iine9  = 0
      dzer9l = 1.0e10
      dzer9p = 1.0e10

      i7   = n7
      i100 = np100
      ivar = npvar
      ikp1 = kp1
      kc20 = kpc20


      rmax   = drmax
      itmax  = itmaxd
      epsrad = depsra
      epsmas = depsma

      pi   = acos(-1.)
      pid2 = pi/2.
      zpi  = 2.*pi
      vpi  = 4.*pi
      zpi3 = zpi**3
      piq  = pi**2

      uf12 = uf1**2
      uf13 = uf1**3
      uf10 = 1./2.9979e23   ! c=3x10^8 m/s, 1 m=10^15 fm => c=3x10^23 fm/s 

      pdrip  = -1.
      edrip  = -1.
      eqdrip = -2.



      open (unit=71, file='prns.dat', status='unknown')


!__ begin of data file
      read ( 71 ,*) dd, bodf


!__ dimension of work arrays
      read ( 71 ,*) dd, j651
      j650   = j651-1
      j651la = m651
      if (j651 > j651la) then
      write (6,'(////,'' fatal input error in _input_'',/,'' j651 >'',
     +                '' m651'',/)')
      write (6,nainp1)
      stop
      end if

!__ compute I(r)  ><  ivonr = 10
      read ( 71 ,*) dd, ivonr


!__ models of strange dwarfs with R_core -> 0   ><   isdwf /= 0
      read ( 71 ,*) dd, isdwf   

      
!__ White dwarf parameter
!     iwhiteDwarf =   0  ><  Perform NS or SS calculations       
!     iwhiteDwarf = -10  ><  Perform Strange Dwarf calculations
!     iwhiteDwarf =  10  ><  Perform White Dwarf calculations       
      read (71, *) dd, iwhiteDwarf

      
!__ Latent heat parameter
!     latent_heat = 0   ><  Table for latent heat calculations will not be created
!     latent_heat /= 0  ><  Table for latent heat calculations will be created
      read (71, *) dd, latent_heat 
      
!__ tabulate   r, z(r), n(r), A(r), e_c   ><   iarm /= 0 
!   (data for Armen Sedrakian's cooling calculations)  
      read ( 71 ,*) dd, iarm 
  

      do 593 i=1,kc20
  593 read ( 71 ,*) dd, r_core(i)

      read ( 71 ,*) dd, ico21a,ico21b
      if (ico21b > kc20) then
      write (6,'(////,'' fatal input error in _input_'',/,
     +                '' ico21b > kc20'',/)')
      stop
      end if


!__ parameters determining the quantities to be written to output files
!   during rotating star calculation
      read ( 71 ,*) dd,iradot,iout20,iotomg,iout33,ioutes,iotmp0,iop30,
     +             iothv2,ioutpe,ioutse
! NOTE:  ioutse=10  ><  calculation of how much mass is sitting in the 
!                       various phases for G_240^B180 eos 
!        ioutse=20  ><  same, but for G_300^B180 eos 
!        ioutse=30  ><  same, but for G_300^B180_m7_a0
!        ...            see table in subroutine input
!        ...
!        If ioutse =/ from values listed in input file -> this part of 
!        the calculation is skipped.

!__ minimum pressure of harrison-wheeler eos
      read ( 71 ,*) dd, ahw, bhw


!__ join hw eos with hv eos at the values given in the following
      read ( 71 ,*) dd, anv, bnv


!__ print output during calculation on screen (if # 0)
      read ( 71 ,*) dd, ioutd4


!__ method to interpolate the eos in subroutines OVNROT, HT1000
!   ipolxx < 0  ><  log
!   ipolxx = 1  ><  linear
!   ipolxx = 2  ><  polynomial
!   ipolxx = 3  ><  spline
      read ( 71 ,*) dd, ipolxx


!__ check causality and microscopic stability if  iccaus # 0
      read ( 71 ,*) dd, iccaus


!__ non-rotating neutron star calculation or not
!   selec0 = 'oppenheimer-volkoff'
!   selec0 = 'empty              '
      read ( 71 ,*) dd, selec0


!__ rotating neutron star calculation or not, performed in the
!   framework of:
!   (1) general raltivity (HT)
!   selec1 = 'hartle-thorne  '
!   (2) classical newtonian physics
!   selec1 = 'empty          ,
      read ( 71 ,*) dd, selec1



!__ OV calculation:
!   ldecov < 0: stop calculation when star mass  M(i+1) < M(i)
!   ldecov > 0: do not stop
      read ( 71 ,*) dd, ldecov

!__ HT calculation:
!   ldecht < 0: stop calculation when star mass  M(i+1) < M(i)
!   ldecht > 0: do not stop
      read ( 71 ,*) dd, ldecht


!__ moment of inertia of a:
!   spherically symmetric star           ><  icmoi = 0
!   rotationally  d e f o r m e d  star  ><  icmoi # 0
!
!   calculation of  I  in subroutine OVNROT   ><   icmi99 > 0
!                                                  icmi99 < 0, otherwise
      read ( 71 ,*) dd, icmoi, icmi99


!__ select eos and eos of crust region
!   - iqm=20  ><   hadronic crust + quark matter core
!   - iqm=10  ><   pure quark matter eos
!   - iqm=0   ><   pure hadronic matter eos (combined or not combined)
!   - iqm<0   ><   quark + hadronic matter eos combined

!   - icrust=10   ><   HW + NV
!   - icrust=20   ><   BPS
      read ( 71 ,*) dd, iqm, icrust


!__ for hadronic crust calculation, select value of the energy density
!   of the crust at the inner surface:
!
!                                             /  = 0  <==>  e = e_drip
!                                     ejgcm3
!                                             \  # 0  <==>  e < e_drip
      read ( 71 ,*) dd, ejgcm3


!__ specify  kst1
!   - kst1=-90 ><   canuto's favored eos [ca74,75],
!   - kst1=-80 ><   n=3/2 relativistic polytrope [tro65],
!   - kst1=-70 ><   wiringa-fiks-fabrocini eos  [wff88], av14+uvii
!   - kst1=-60 ><   wiringa-fiks-fabrocini eos  [wff88], uv14+uvii
!   - kst1=-50 ><   wiringa-fiks-fabrocini eos  [wff88], uv14+tni
!   - kst1=-40 ><   friedman-pandharipande eos  [fp81], combined
!   - kst1=-30 ><   Pandharipande eos C  [bps71], combined
!   - kst1=-20 ><   bethe-johnson eos  [rd84], model  I, combined
!   - kst1=-10 ><   harrison-wheeler eos only (according to  [ht68])
!   - kst1=100 ><   different eos will be calculated and appropriately
!                   combined (i.e., hw-, nv-, hadronic-, and quark
!                   matter eos)
!   - kst1=110 ><   hadronic matter (crust) + quark matter core
!   - kst1=200 ><   input white dwarf eos from an external file

      read ( 71 ,*) dd, kst1


!__ bag constant (in MeV) and maximum central energy density of the
!   quark star: [ecqmax] = MeV/fm^3               (subroutine SEOS2)
!               [rcqmax] = rho/rho(nucl. matter)   (subroutine SEOS1)
      read ( 71 ,*) dd, bagct, ecqmax, rcqmax


!__ scaling factor  r0r0  (fm)
!   r0r0=3.0e19  if neutron stars having radii < 20 km are concerned;
!   otherwise a larger value must be chosen
      read ( 71 ,*) dd, r0r0


!                               /\
!__ lower and upper bounds  of   r = [bord0,bord1]
      read ( 71 ,*) dd, bord0, bord1


!__ parameters determining the step size in the case of quark matter
!   stars possessing a nuclear crust;
!   note: iqm2020 < j651, and  r5km  denotes that radial distance at
!         which the step size changes in favor of a larger grid spacing
!         (typically, r5km=5 km)
      read ( 71 ,*) dd, iq2020, r5km
      if (j651 < (2*iq2020)) then
      write (6,'(////,'' input error in _input_'',/,'' j651 <'',
     +                '' 2*iq2020'',/)')
      iq202l=iq2020
      r5kml =r5km
      write (6,nainp1)
      stop
      end if

!__ maximum difference of:
!                           / (e(c;i+1)-e(c;i))/e(c;i) / < etasc6
!                           / (omega(i+1)-omega(i))/omega(i) / < etasc7
      read ( 71 ,*) dd, etasc6, etasc7


!__ maximum number of iteration steps
      read ( 71 ,*) dd, igivmx, iscmx7


!__ central energy densities (in MeV/fm**3)
      read ( 71 ,*) dd, i27
      do 200 i=1,i27
      read ( 71 ,*) dd, idummy, zedy(i)
  200 continue

      
!__ how many stellar models?   
      read ( 71 ,*) dd, i4006a, i4006b

      
!__ compute stars with a given
!                                  gravitatinal mass >< igravBary = 1
!                                  baryon mass       >< igravBary = 2
!                                  baryon number     >< igravBary = 3 
      read ( 71 ,*) dd, igravBary 

      
!__ gravitational neutron star masses
      do 202 i=1,i27
      read ( 71 ,*) dd, idummy, fixmov(i)
  202 continue



!__ rotational frequencies (in  1/sec) as seen by an observer at infinity
      read ( 71 ,*) dd, iog40
      read ( 71 ,*) dd, inpog1
      do 300 i=1,iog40
  300 read ( 71 ,*) dd, idummy, rlomg(i)


!__ pcut:  change from spline to logarithmic or spline-douple precision
!          interpolation for  p < pcut (in MeV/fm^3)

      read ( 71 ,*) dd, pcut


!__ pcuts: stop integration of stellar structure equations for surface
!          pressure values  p < P(e=7.86 g/cm^3) = pcuts (in MeV/fm^3)

      read ( 71 ,*) dd, pcuts


!__ integration of OV and HT equations via Euler or fourth order
!   Runge-Kutta method
!     I) OV equations
!     irukov = 0  ><  Euler method
!     irukov # 0  ><  Runge-Kutta
!     I) HT equations
!     irukht = 0  ><  Euler method
!     irukht # 0  ><  Runge-Kutta

      read ( 71 ,*) dd, irukov,irukht



      read ( 71 ,*) dd, eodf

      if ( (eodf /= ('end of data file'))   .or.
     +     (bodf /= ('begin of data file'))       ) then
      write(6,'(//,1x,''fatal input error'',/)')
      do 401 i=1,i27
  401 write(6,'('' i='',i3,3x,''eps(c)='',e14.8)') i,zedy(i)
      do 402 i=1,iog40
  402 write(6,'('' i='',i3,6x,''omega='',e14.8)') i,rlomg(i)
      stop
      end if


      close ( 71 )


      if ( (ivar /= (3*i100)) .or. (inpog1 > iog40) .or.
     +     (iog40 /= npomg)  .or. (i27 /= n27) ) then
      print*,' _input_  -->  fatal error detection'
      print*,' i100=',i100,'  ivar=',ivar,'  inpog1=',inpog1
      print*,' npomg=',npomg,'  iog40=',iog40,'  i27=',i27
      print*,' n27=',n27
      stop 'error stop'
      end if


!__ check whether elements in array  rlomg  are monotoniously increasing
      call sortor(rlomg,npomg,-1,ierso)
      if (ierso /= 0) then
      write (6,'(/////,'' Warning:  array rlomg not ordered! '')')
      do 177 i=1,npomg
  177 write (6,178) i,rlomg(i)
  178 format (' i=',i3,3x,' omega=',e14.6,' 1/sec')
!      call exit
      end if


      return

      end

!resov
      subroutine resov(mmm,rmqoz,rrns,wi,irem,ec,iarm)

! ----------------------------------------------------------------------
!
! purpose:
!          computation of the results
!
! foregoing calculation: oppenheimer-volkoff
!
! ----------------------------------------------------------------------


      real * 8 rmsun,grav,gravk,gravkm,ee18,gcm2km,eem18,ar0,
     +         ee14,ee34,ee44
      real * 8 ee03,ee13,ee15
      

      parameter (n27=50)
      common/zzzz/zedy(n27),zmns(n27),zrns(n27),ztm(n27),zkepf(n27),
     +            zfrp(n27),fixmov(n27)
      common/zzz1/zedyov(n27),zedynw(n27),zedyht(n27)
      common/PowersTen/ee03,ee15,ee18,eem18,ee13,ee14,ee34,ee44
      common/MultiPi/pi,piq,zpi,vpi,zpi3,pid2
      common/RadialScalingFactor/r0r0

      common/SolarProperties/rmsun,rsun,rsunkm
      common/Gravitational/gravk,grav,gravkm,gcm2km
      common/ConversionFactors/uf1,uf12,uf13,uf4,uf6,uf10

      parameter (m651=160501)
      common/cnsta1/f251(m651),fout(m651),fst(m651),fst251(m651),feng1
     +              (m651),x251(m651),y251(m651),f251pm(m651),fdrip(m6
     +              51)
      common/bn01/axbn1(m651),axbn2(m651),pdenex(m651),einthx(m651),
     +            deindp(m651),dehdrh(m651),pdent(m651)
      common/cnsta2/aux51(1,m651),aux5(m651)
      common/cnsta4/aux71(1,m651),aux7(m651),aux8(m651),aux9(m651)


      open (unit=19, file='nzr.dat', status='unknown')

      ar0 = r0r0

      zmns(mmm) = rmqoz
      zrns(mmm) = rrns
      ztm(mmm)  = wi

      a2 = zmns(mmm)
      a3 = zrns(mmm)/ee03
      a4 = ztm(mmm)

!     fractional redshift (a5)
      h1 = 2.*zmns(mmm)*rmsun*gravkm/a3
      h1 = 1.-h1
      a5 = 1./sqrt(h1)-1.
      zfrp(mmm) = a5

!     maximum keplerian angular velocity (rad/sec) as determined by Friedman, Ipser & Parker
      a6 = 24.*sqrt(a2/a3**3)*1.0e04
      zkepf(mmm) = a6


!     write r, z(r), n(r), N(r), and e(r) on output file
      iarm = 0
      if (iarm /= 0) then 
      write (19,302)
  302 format (/,4x,'r (km)',9x,'z',10x,'n (1/fm^3)',10x,'A',11x,'e (MeV',
     +             '/fm^3)')
      write (19,888)
  888 format (' ---------------------------------------------------',
     +        '----------------------------')
      do 10 i=1,irem
      rr = fst251(i)*r0r0/ee18
      sr = sqrt(1.-2.*f251(i)*rmsun*gravkm/rr)
      z  = 1./sr - 1.
! red shift z(r)
      aux5(i)=z

      as = 3.*dlog10(ar0) + alog10(vpi*axbn1(i))

! total number of baryons within sphere of radius r, A(r)
      aux8(i)=as

! baryon density at radial distance r, n(r), in 1/fm3
      aux7(i)=pdenex(i)

! energy density at radial distance r, e(r), in MeV/fm^3
      aux9(i)=feng1(i)*ec
      write (19,20) rr,aux5(i),aux7(i),aux8(i),aux9(i)
 20   format(2x,f9.5,3x,e11.5,3x,e10.4,3x,f15.11,3x,e12.6) 
   10 continue

      end if     

      close ( 19 )



      return

      end

!resrn
      subroutine resrn(mmm,rmqoz,rth,ni,wi,omega,ec)

! ----------------------------------------------------------------------
!
! purpose:
!          computation of the results
!
! foregoing calculation: newtonian treatment
!
! ----------------------------------------------------------------------


      real * 8 rmsun,grav,gravk,gravkm,ee18,gcm2km,eem18,ee14,ee34,
     +         ee44
      real * 8 ee03,ee13,ee15
      
      real rth(0:ni,4)

      parameter (n27=50,n10=15,n17=50,np4=26,npapp=10)
      common/tabn1/tn1(n27,n10)
      common/tabht1/tht1(n27,n17),tht2(n27,np4),tappx(n27,npapp)

      common/PowersTen/ee03,ee15,ee18,eem18,ee13,ee14,ee34,ee44
      common/MultiPi/pi,piq,zpi,vpi,zpi3,pid2
      common/SolarProperties/rmsun,rsun,rsunkm
      common/Gravitational/gravk,grav,gravkm,gcm2km
      common/ConversionFactors/uf1,uf12,uf13,uf4,uf6,uf10

      a3r = rth(ni,2)/ee03
      a33r = rth(0,2)/ee03

      tn1(mmm,1) =ec
      tn1(mmm,2) =ec*uf6
      tn1(mmm,3) =alog10(tn1(mmm,2))
      tn1(mmm,13)=a3r
      tn1(mmm,4) =a33r
      tn1(mmm,7) =rmqoz
      tn1(mmm,10)=wi

!     fractional redshift (a5r)
      h1r = 2.*rmqoz*rmsun*gravkm/a33r
      h1r = 1.-h1r
      if (h1r <= 0.) then
      write (6,'(/,'' _resrn_ -info-  h1r<0'',/)')
      a5r = 0.
      else
      a5r = 1./sqrt(h1r)-1.
      end if
      tn1(mmm,11)=a5r

!     eccentricity
      if (rth(ni,2) /= 0.) then
      wz2=1.-(rth(0,2)/rth(ni,2))**2
      ecc=sqrt(wz2)
      else
      ecc=-1.
      end if
      tn1(mmm,12)=ecc

      rfreq=omega/uf10
      tn1(mmm,5) =rfreq

      return

      end


!resrht
      subroutine resrht(mmm,rth,ni,ec)

! ----------------------------------------------------------------------
!
! purpose:
!          determine results of HT calculation
!
! ----------------------------------------------------------------------


      real * 8 rmsun,grav,gravk,gravkm,ee18,gcm2km,eem18,ee14,ee34,ee44
      real * 8 ee03,ee13,ee15
      
      real rth(0:ni,4)

      parameter (n27=50,n17=50,np4=26,npapp=10)
      common/tabht1/tht1(n27,n17),tht2(n27,np4),tappx(n27,npapp)

      parameter (nptb3=10)
      common/tabbbb/tabb(n27,nptb3)

      parameter (nttby2=15)
      common/tbynb1/ttby(n27,nttby2)

      common/PowersTen/ee03,ee15,ee18,eem18,ee13,ee14,ee34,ee44
      common/ConversionFactor30/ufe30
      common/MultiPi/pi,piq,zpi,vpi,zpi3,pid2
      common/SolarProperties/rmsun,rsun,rsunkm
      common/NucleonMasses/amnuc,amneut
      common/Gravitational/gravk,grav,gravkm,gcm2km
      common/ConversionFactors/uf1,uf12,uf13,uf4,uf6,uf10
      common/redsep/h0rcbe,h0rcbp,phircp,phirce,rdcbep,rdcbee,odcbee,
     +              h2rcbe,h2rcbp,v2rcbe,v2rcbp
      
      tht1(mmm,1)=ec
      tht1(mmm,2)=ec*uf6
      tht1(mmm,3)=alog10(tht1(mmm,2))


      do 100 i=1,2
!       -  i=1  ><  fractional redshift at the star's pole
!       -  i=2  ><  fractional redshift at the star's equator

!     radius
      if (i == 1) a33r = tht1(mmm,15)
      if (i == 2) a33r = tht1(mmm,16)

      h1r = 2.*tht1(mmm,11)*rmsun*gravkm/a33r
      h1r = 1.-h1r

!RRR recalculate using the definition given by  kapoor and datta
!RRR  if (i == 1) h1r = exp(phircp) * (1.+2.*(h0rcbp+h2rcbp))
      if (i == 2) h1r =   exp(phirce) * (1.+2.*(h0rcbe-h2rcbe/2.))
     +                  - (rdcbee*odcbee)**2 * (1.-v2rcbe+h2rcbe)

      if (h1r <= 0.) then
      write (6,'(//,'' _resrht_ -info-  h1r<0'',//)')
      a5r = 0.
      else
      a5r = 1./sqrt(h1r)-1.
      end if

      if (i == 1) tht1(mmm,23)=a5r
      if (i == 2) tht1(mmm,24)=a5r
  100 continue

!   critical angular velocity Omega (in 1/sec) as used by Datta
!   in [Dat88] (i.e., onset of the secular instability in a rotating
!   Maclaurin spheroid; eccentricity = 0.81267)
      ff=5984.97
      r10km=tht1(mmm,5)/10.
      srt=sqrt(tht1(mmm,10) / r10km**3)
      tabb(mmm,4)=ff*srt
!   same as above but for an onset of the dynamical instability at e=0.95289
      ff=6616.63
      tabb(mmm,10)=ff*srt

!   baryon mass of 
!   a) spherically symmetric neutron star:
      bymslg=alog10(amneut) + ttby(mmm,2) - dlog10(rmsun)-alog10(ufe30)
      bymsph=10.**bymslg
      ttby(mmm,8)=bymsph
      
!   b) non-spherical (i.e., rotating) neutron star:
      bymrlg=alog10(amneut) + ttby(mmm,4) - dlog10(rmsun)-alog10(ufe30)
      bymrot=10.**bymrlg
      ttby(mmm,7)=bymrot

!   injection energy
      beta = 1. / (tht1(mmm,23)+1.)**2
      tht2(mmm,1)=beta

!   ratio of rotational energy to gravitational energy: t/w
      tdw= 1. / ( 1. + (tht2(mmm,14)-tht1(mmm,11))/tht2(mmm,2) )
      tht2(mmm,5)=tdw


!   binding energy of non-rotating star (in units of the solar mass),
!   sign definition according to hartle [1967]
      ttby(mmm,9)=ttby(mmm,8) - tht1(mmm,10)
      if (ttby(mmm,8) == 0.) ttby(mmm,9)=0.

!   binding energy of rotating star (in units of the solar mass),
!   sign definition according to hartle [1967]
      tht2(mmm,6)=tht2(mmm,14) + tht2(mmm,2) - tht1(mmm,11)

!cc  corrected expression:
      tht2(mmm,6)=ttby(mmm,7) - tht2(mmm,2) - tht1(mmm,11)
      if (ttby(mmm,8) == 0.) tht2(mmm,6)=0.

!     packing coefficients: I) non-rotating star
!     alpha_1
      ttby(mmm,10)=(tht2(mmm,4) - tht1(mmm,10)) / tht1(mmm,10)
      if (tht2(mmm,4) == 0.) ttby(mmm,10)=0.
!     alpha_2
      ttby(mmm,11)=(ttby(mmm,8) - tht1(mmm,10)) / tht1(mmm,10)
      if (ttby(mmm,8) == 0.) ttby(mmm,11)=0.

!     packing coefficients: II) rotating star
!     alpha_1
      ttby(mmm,12)=(tht2(mmm,14) - tht1(mmm,11) + tht2(mmm,2)) /
     +             tht1(mmm,11)
      if (tht2(mmm,14) == 0.) ttby(mmm,12)=0.
!     alpha_2
      ttby(mmm,13)=(ttby(mmm,7) - tht1(mmm,11) - tht2(mmm,2))  /
     +             tht1(mmm,11)
      if (ttby(mmm,7) == 0.) ttby(mmm,13)=0.



!     compute interglitch period (according to star quake model for glitches, see
!     Diaz and Ibanez, ApJ 291 (1985) 308)
      rig=3.553e+06
      
!   I) non-rotating star
      xz=alog10(rig * (ttby(mmm,10)*tht1(mmm,10))**2)
      xn=3.*alog10(tht1(mmm,5)) + tht1(mmm,19) - 45.
      xd=xz-xn
      tq=10**xd
      ttby(mmm,14)=tq
      
!   II) rotating star
      xz=alog10(rig * (ttby(mmm,12)*tht1(mmm,11))**2)
      xn=3.*alog10(tht1(mmm,16)) + tht1(mmm,19) - 45.
      xd=xz-xn
      tq=10**xd
      ttby(mmm,15)=tq


!   surface dipole magnetic field (in units of 10**8 Gauss)
      xmdot=1.
      rus  =tht1(mmm,16)/10.
      b8min=0.9 * tht1(mmm,11)**0.25 * sqrt(xmdot) / rus**1.25
      perod=zpi/tabb(mmm,1)
      per76=(perod*1.e04 / 6.)**(7./6.)
      b8max=per76 * tht1(mmm,11)**(5./6.) * sqrt(xmdot) / rus**3
      tabb(mmm,2)=b8min
      tabb(mmm,3)=b8max

!   compute   d P/d t  and  log [d P(s s^-1)/d t]
      yilg  =tht1(mmm,19)
      yilg45=yilg-45.
      delg45=10.**yilg45
      ymrons=tht1(mmm,11)
      ym53  =ymrons**(5./3.)
      per43 =(1.0e03*perod)**(4./3.)
      ylgpdo=-20. + alog10(3.2 * ym53 * xmdot * per43 ) - yilg45
      pdo18  =10.**ylgpdo
      tabb(mmm,5)=ylgpdo
      tabb(mmm,7)=pdo18

!   compute   log (P (s))  and  T(s) (lifetime of pulsar in years)
      tabb(mmm,8)=perod*1.0e03
      tabb(mmm,6)=alog10(perod)

      per13 =(1.0e03*perod)**(1./3.)
      ts=1.0e09 * delg45 / (xmdot * ym53 * per13)
      tabb(mmm,9)=ts

      return

      end

!outov
      subroutine outov(nj,ni,iqm,iov1,kst1,selec0,irem,bagct,xipov,i5
     +                 651,bord0,bord1,ejgcm3,edrip,sslt,ssgt,r5km,iq
     +                 2020,prkov,popmfm,iccorr,isdwf,ico21a,ico21b,v
     +                 ersion,iwhiteDwarf,latent_heat,igravBary,ahw,b
     +                 hw,anv,bnv)

! ----------------------------------------------------------------------
!
! purpose:
!     output results of  n o n - rotating compact star calculation
!     (oppenheimer-volkoff treatment) 
!
! ----------------------------------------------------------------------

      real * 8 ee18,eem18,ee14,ee34,ee44
      real * 8 ee03,ee13,ee15
      
      parameter (n27=50)
      common/zzzz/zedy(n27),zmns(n27),zrns(n27),ztm(n27),zkepf(n27),
     +            zfrp(n27),fixmov(n27)
      common/zzz1/zedyov(n27),zedynw(n27),zedyht(n27)

      parameter (npovdp=4)
      common/tabccc/tovdp(n27,npovdp)

      parameter (nthov1=12)
      common/tabsph/thov1(n27,nthov1)

      parameter (kpc20=20,kpc5=6)
      common/tc/tcore(kpc20,kpc5),r_core(kpc20)


      common/ConversionFactors/uf1,uf12,uf13,uf4,uf6,uf10
      common/MiscInt1/i100,i7,i27,ivar,inpog1,ikp1,icrust
      common/MiscInt2/j651,j650,igivmx,ivonr,kc20
      common/ConvergencyLimits/etasc6,etasc7,iscmx7
      common/ctype/typeos,teoscb
      common/NuclearMatterDensity/engnm0
      common/PowersTen/ee03,ee15,ee18,eem18,ee13,ee14,ee34,ee44
      common/RadialScalingFactor/r0r0

      character * 19 selec0
      character * 25 xipov
      character * 50 typeos
      character * 45 teoscb
      character * 38 prkov
      character * 58 version

      real start, finish

      character (8) date
      character (10) time 

      call cpu_time(start)

      open (unit=66, file='outputTOV.dat', status='unknown')

      r0r0km = r0r0/ee18

      write (66,333)
  333 format (' *************************************************',
     +     '**************************')

      if(iwhiteDwarf == 0)  write (66,'('' properties of non-rotating'',
     +                                  '' compact (TOV) stars'')')
      if(iwhiteDwarf == 10) write (66,'('' properties of non-rotating'',
     +                                  '' white dwarfs (using TOV)'')')
      if(iwhiteDwarf  < 0)  write (66,'('' properties of rotating str'',
     +                                  ''ange (white) dwarfs'')')
      write (66,'('' oppenheimer-volkoff treatment'')')
      write (66,'(1x,a58)') version
      if (iwhiteDwarf /= 0 .and. kst1 == 200) then
       write (66,'('' input eos from external file (prns_eos.dat)'')')
       write (66,'('' eos at supernuclear densities: '', a50)') typeos
      else
       write (66,'('' eos at supernuclear densities: '', a50)') typeos
       write (66,'('' eos at subnuclear densities: '', a45)') teoscb
      end if
       if (iccorr /= 0)
     +write (66,'('' eos has been corrected for causality'')')
      write (66,333)

      if (iqm  < 0)   write (66,'('' combine quark + hadronic eos'')')
      if (iqm  == 0)  write (66,'('' hadronic eos only'')')
      if (iqm  == 10) write (66,'('' quark matter eos only'')')
      if (iqm  == 20) then
      if (ejgcm3 == 0.) then
      write (66,'('' hadronic crust + quark matter core: e_join='',
     +            e12.6,'' g/cm^3    = e_drip'')') edrip*uf6
      write (66,'(''                                     P_join='',
     +            e12.6,'' MeV/fm^3  = P_drip'')') popmfm
      else
      write (66,'('' hadronic crust + quark matter core: e_join='',
     +            e12.6,'' g/cm^3    < e_drip'')') ejgcm3
      write (66,'(''                                     P_join='',
     +            e12.6,'' MeV/fm^3  < P_drip'')') popmfm
      end if
              end if

      if (iov1 < 0) write (66,'('' rotating star: pseudo oppenhei'',
     +                           ''mer volkoff (pov) calculation'')')
      if (iov1 > 0) write (66,'('' rotating star: newtonian calcu'',
     +                           ''lation'')')

      if (i5651 /= 1 .and. iscmx7 == 1)
     +write (66,'(/,'' perform compact star calculation (TOV) of a'',
     +              '' star with a given mass'',/,'' or baryon number'', 
     +              '' => e_c self-consistent (root-finding)'')')

      if (iscmx7 == 1 .and. igravBary == 1) then
         write(66,'('' -> gravitational mass kept constant'')')
      else if (iscmx7 == 1 .and. igravBary == 2) then
         write(66,'('' -> baryon mass kept constant'')')
      else if (iscmx7 == 1 .and. igravBary == 3) then
         write(66,'('' -> baryon number (log_10 A) kept constant'')')
      end if
      
      if (i5651 == 1 .and. iscmx7 == 1) then
      write (66,'(///,'' FATAL INPUT ERROR !!!'',//,'' calculation of'',
     +                '' the properties of a rotating neutron/quark'',/,
     +                '' star for given values of'',
     +                '' e_c(i)  and  Omega_star(i) not possible '')')
      go to 3412
      end if


! determine and print CPU time      
      call cpu_time(finish)
      write(66,'(/, " Elapsed CPU Time: ", F10.3," seconds")')
     +     finish-start

! determine and print Date and Time
      call date_and_time(date, time)
      write(66,'(" Date: ", A, ",", 1X, "Time: ", A, /)') date, time

      
! write input-parameters:
      write (66,'(//,'' input parameters:'')')
      write (66,788) j651,nj,ni,iqm,iov1,kst1,selec0,irem,bagct,icrust,
     +               xipov,r0r0,r0r0km,bord0,bord1,iq2020,r5km,prkov,
     +               latent_heat,igravBary,zedy(3),zedy(4),ahw,bhw,anv,
     +               bnv
 788  format (19x,'j651=',i6,/,19x,'nj=',i4,';  ni=',i3,/,19x,'iqm=',
     +     i3,/,19x,'iov1=',i3,/,19x,'kst1=',i4,/,19x,'treatment: ',
     +     a19,/,19x,'irem=',i6,/,19x,'bag constant =',e12.6,' MeV',
     +     /,19x,'icrust=',i3,/,19x,'interpolation technique: ',
     +     a25,/,19x,'r0r0=',e12.6,' (=> R_max=',f8.2,' km)',
     +     /,19x,'bord0=',e12.6 ,' bord1=',e12.6,/,19x,'iq2020=',i6,
     +     '  r5km=',e12.6,' km',/,19x,a38,/,19x,'latent_heat=',i3,2x,
     +     'igravBary=',i2,/,19x,'e_c(1)=',f7.2,'  e_c(2)=',f7.2,
     +     ' MeV/fm^3',/,' Harrison-Wheeler boundaries:',
     +     ' a=',e12.6,2x,' b=',e12.6,' g/cm^3',/,' Negele-Vautheri',
     +     'n boundaries:',' a=',f10.7,2x,' b=', f10.6,' MeV/fm^3')

      if (edrip == (-1.)) then
      write (66,'(1x,''step size for radial integration: '',e12.6,'' m''
     +            ,//)') sslt
      else if (edrip /= (-1.)) then
      write(66,'(1x,''step size for  r < '',e12.6,'' km: '',f7.4,'' m'')
     +          ') r5km,sslt
      write(66,'(1x,''step size for  r > '',e12.6,'' km: '',f12.6,'' m''
     +            ,//)') r5km,ssgt
      end if


      write (66,302)
  302 format (/,3x,'   e (MeV/fm^3)     e/e_0       M/M_s         R (k',
     +             'm)  log I/g cm^2    z')
      write (66,888)
  888 format (' ---------------------------------------------------',
     +        '----------------------------')
  899 format (' ---------------------------------------------------',
     +        '----------------------------------------------------',
     +        '---------')

      do 1000 m=1,i27
      if (zmns(m) == 0.) go to 1000
      a1=zedyov(m)
      a3=zrns(m)/ee03
      eunm=zedyov(m)/engnm0
      write (66,202) m,a1,eunm,zmns(m),a3,ztm(m),zfrp(m)
 1000 continue
  202 format (1x,i3,1x,e15.9,1x,e9.3,1x,e14.8,1x,f11.5,1x,f9.6,2x,e11.5)



      write (66,402)
  402 format (//,4x,'R_core/km   M28_d-s/g    M_d-s/M_sun   M_prop/',
     +           'M_sun Om_K_FIP (1/s) M_c/M_sun')
      write (66,888)

      do 2000 m=1,i27
      if (zmns(m) == 0.) go to 2000
      diffm=thov1(m,3)-tovdp(m,3)
      write (66,203) m,tovdp(m,1),tovdp(m,2),tovdp(m,3),tovdp(m,4),
     +                 zkepf(m),diffm
 2000 continue
  203 format (1x,i3,e10.3,1x,e13.7,1x,e13.7,2x,f10.7,6x,f8.2,2x,e10.4)



      write (66,422)
  422 format (//,7x,' lg_10 A_sph   lg_10 A_core   lg_10 A_crust',
     +              '     n_c     P_cut (MeV/fm^3)')
      write (66,888)

      do 2022 m=1,i27
      if (zmns(m) == 0.) go to 2022
      write (66,223) m,thov1(m,7),thov1(m,8),thov1(m,9),thov1(m,10),
     +                 thov1(m,11)
 2022 continue
  223 format (1x,i3,3x,f12.8,3x,f12.8,3x,f12.8,3x,e10.4,4x,e12.5)


      if (isdwf == 0) go to 3412

      write (66,1402)
 1402 format (//,4x,' R_core/m          M/M_sun        M_core/M_sun  ',
     +            'M_crust/M_sun',10x,'R (km)',5x,'  lg A_core    lg ',
     +            'A_crust')
      write (66,899)

      do 2028 m=ico21a,ico21b
      write (66,1203) m,r_core(m),tcore(m,1),tcore(m,2),tcore(m,3),
     +                  tcore(m,4),tcore(m,5),tcore(m,6)
 2028 continue
 1203 format (1x,i3,e10.3,2x,e20.13,2x,e12.6,2x,e12.6,2x,f19.13,2x,
     +           f11.8,3x,f11.8)



      close ( 66 )

 3412 return

      end

!outrn
      subroutine outrn(iov1,omega,irem)

! ----------------------------------------------------------------------
!
! purpose:
!          output result of   r o t a t i n g  compact star calculation
!          perfomed in the (pseudo) newtonian limit.
!
! ----------------------------------------------------------------------


      parameter (n27=50,n10=15)
      common/tabn1/tn1(n27,n10)

      common/MiscInt1/i100,i7,i27,ivar,inpog1,ikp1,icrust
      common/ctype/typeos,teoscb

      character * 50 typeos
      character * 45 teoscb


      open (unit=41, file='rnsnew.dat', status='unknown')

      write (41,333)
  333 format (/,' *************************************************',
     +        '**************************',/)

      write (41,'('' properties of rotating (non-spherical)'',
     +           '' neutron/hybrid/quark star models'',/)')

      if (iov1 < 0) then
      write (41,'('' (pseudo) oppenheimer-volkoff treatment'')')
      else if (iov1 > 0) then
      write (41,'('' classical newtonian calculation'')')
      end if
      write (41,'('' eos at supernuclear densities: '',a50)') typeos
      write (41,'('' eos at subnuclear densities: '',a45)') teoscb
      write (41,333)

      write (41,410) omega,irem
  410 format (/,' omega=',e14.6,'1/sec',/,' irem=',i6,/)

      write (41,802)
  802 format (/,8x,'e_c   MeV/fm^3  g/cm^3   log  r_sph  km  og (1/s)',
     +        '  og_shedding ')
      write (41,888)
  888 format (' ---------------------------------------------------',
     +        '----------------------------')

      do 1000 m=1,i27
      write (41,800) m,tn1(m,1),tn1(m,2),tn1(m,3),tn1(m,4),tn1(m,5),
     +                 tn1(m,6)
 1000 continue
  800 format (1x,i3,2x,e13.6,2x,e13.6,2x,f5.2,2x,f5.2,2x,f8.2,2x,f8.2)


      write (41,902)
  902 format (//,7x,'m_sp/m_sun    m_tot/m_sun    del m/m_sp  log (I)',
     +        '  z   ecc  r_eq  (km)')
      write (41,888)

      do 2000 m=1,i27
      write (41,900) m,tn1(m,7),tn1(m,8),tn1(m,9),tn1(m,10),tn1(m,11),
     +                 tn1(m,12),tn1(m,13)
 2000 continue
  900 format (1x,i3,2x,f6.2,2x,f6.3,2x,f6.4,2x,f5.2,2x,f5.4,2x,f4.2,2x,
     +     f5.2)

      close ( 41 )

      return

      end

!outrht
      subroutine outrht(irem,sum124,bagct,ogc,xipht,icmoi,sumw1,i5651,
     +                  iqm,bord0,bord1,ejgcm3,edrip,sslt,ssgt,r5km,iq
     +                  2020,prkht,popmfm,iccorr,kst1,ec,version,cthr,
     +                  iwhiteDwarf,latent_heat,igravBary,i4006a,i4006b,
     +                  ahw,bhw,anv,bnv)

! ----------------------------------------------------------------------
!
! purpose:
!          output result of a  r o t a t i n g  compact star calculation
!          perfomed for the hartle-thorne treatment.
!
! ----------------------------------------------------------------------

      real * 8 rmsun,ee18,eem18,ee14,ee34,ee44
      real * 8 ee03,ee13,ee15
      
      parameter (n27=50,n17=50,np4=26,npapp=10)
      common/tabht1/tht1(n27,n17),tht2(n27,np4),tappx(n27,npapp)
      common/tabquad/tquad(n27,npapp)
      common/tabstr/tstrob(8,n27)
      common/tabmm/tabmm(8,n27)

      parameter (nttby2=15)
      common/tbynb1/ttby(n27,nttby2)

      parameter (nthov1=12)
      common/tabsph/thov1(n27,nthov1)

      parameter (nptb3=10)
      common/tabbbb/tabb(n27,nptb3)

      parameter (m651=160501)
      common/spl1/st(m651),yout(m651),xstrob(m651),ystrob(m651)
      common/bn01/axbn1(m651),axbn2(m651),pdenex(m651),einthx(m651),
     +            deindp(m651),dehdrh(m651),pdent(m651)
      common/cnsta1/f251(m651),fout(m651),fst(m651),fst251(m651),feng1
     +              (m651),x251(m651),y251(m651),f251pm(m651),fdrip(m6
     +              51)
      COMMON/ZZZZ/ZEDY(N27),ZMNS(N27),ZRNS(N27),ZTM(N27),ZKEPF(N27),
     +            ZFRP(N27),FIXMOV(N27)

      common/iplot1/iradot,iout20,iotomg,iout33,ioutes,iotmp0,iop30,ioth
     +              v2,ioutpe,ioutse
      common/MiscInt1/i100,i7,i27,ivar,inpog1,ikp1,icrust
      common/MiscInt2/j651,j650,igivmx,ivonr,kc20
      common/ConversionFactors/uf1,uf12,uf13,uf4,uf6,uf10
      common/ConvergencyLimits/etasc6,etasc7,iscmx7
      common/ctype/typeos,teoscb
      common/NuclearMatterDensity/engnm0
      common/RadialScalingFactor/r0r0
      common/PowersTen/ee03,ee15,ee18,eem18,ee13,ee14,ee34,ee44
      common/SolarProperties/rmsun,rsun,rsunkm
      common/Thresholds/thresh(15,2),rmtotm,rmtotq,nthrcb
      common/BlobsRodsSlabs/thresh_brs(15,6)
      common/MultiPi/pi,piq,zpi,vpi,zpi3,pid2
      common/GivenInteriorDensity/einter
      
      character * 25 xipht
      character * 50 typeos
      character * 45 teoscb
      character * 38 prkht
      character * 58 version
      character * 30 cthr
      character * 3  convergent
      character * 10 blobs_rods_slabs
      
      real start, finish

      character (8) date
      character (10) time 
      
      call cpu_time(start)
      
      open (unit=43, file='outputHT.dat', status='unknown')


      r0r0km = r0r0/ee18

      write (43,333)
  333 format (' *************************************************',
     +        '**************************')

      if(iwhiteDwarf == 0)  write (43,'('' properties of rotating com'',
     +                                  ''pact stars'')')
      if(iwhiteDwarf == 10) write (43,'('' properties of rotating whi'',
     +                                  ''te dwarfs'')')
      if(iwhiteDwarf  < 0)  write (43,'('' properties of rotating str'',
     +                                  ''ange (white) dwarfs'')')
      write (43,'('' hartle-thorne  treatment'')')
      write (43,'(1x,a58)') version
      if (iwhiteDwarf /= 0 .and. kst1 == 200) then
       write (66,'('' input eos from external file (prns_eos.dat)'')')
       write (66,'('' eos at supernuclear densities: '', a50)') typeos
      else  
       write (43,'('' eos at supernuclear densities: '',a50)') typeos
       write (43,'('' eos at subnuclear densities: '',a45)') teoscb
      end if
      if (iccorr /= 0)
     +write (43,'('' eos has been corrected for causality'')')

      write (43,333)

      if (iqm  < 0)   write (43,'('' combine quark + hadronic eos'')')
      if (iqm  == 0)  write (43,'('' hadronic eos only'')')
      if (iqm  == 10) write (43,'('' quark matter eos only'')')
      if (iqm  == 20) then
      if (ejgcm3 == 0.) then
      write (43,'('' hadronic crust + quark matter core: e_join='',
     +            e12.6,'' g/cm^3    = e_drip'')') edrip*uf6
      write (43,'(''                                     P_join='',
     +            e12.6,'' MeV/fm^3  = P_drip'')') popmfm
      else
      write (43,'('' hadronic crust + quark matter core: e_join='',
     +            e12.6,'' g/cm^3    < e_drip'')') ejgcm3
      write (43,'(''                                     P_join='',
     +            e12.6,'' MeV/fm^3  < P_drip'')') popmfm
      end if
              end if

      if (i5651 /= 1 .and. iscmx7 == 1)
     +   write (43,'(/,'' rotating NS/QS calculation for given mass'',
     +                 '' or baryon number and rotational frequency'', 
     +               /,'' => e_c self-consistent (root finding)'')')
      if (i5651 == 1 .and. iscmx7 == 1)
     +   write (43,'(/,'' compute the properties of rotating neutron/'',
     +             ''quark stars for given'', /, '' values of '',
     +             '' e_c and Omega => M(e_c;Omega)'', //,
     +             '' NOTE: not a self-consistent calculation '',/)')      

      if (i5651 == 1 .and. iscmx7 /= 1)
     +   write (43,'(/,'' perform NS/QS calculation for a star rotat'',
     +                 ''ing at the general relativistic'',/, '' Kep'',
     +                 ''ler frequency => O_K self-consistent'')')

      if (iscmx7 == 1 .and. igravBary == 1) then
         write(43,'('' -> gravitational mass kept constant'')')
      else if (iscmx7 == 1 .and. igravBary == 2) then
         write(43,'('' -> baryon mass kept constant'')')
      else if (iscmx7 == 1 .and. igravBary == 3) then
         write(43,'('' -> baryon number (log_10 A) kept constant'')')
      end if
      
      
      if (ioutse == 0 ) 
     + write (43,'(/, '' NOTE: ioutse=0 -> particle threholds not'',
     +                '' provided for this nuclear EoS'', /)')      

! determine and print CPU time      
      call cpu_time(finish)
      write(43,'(/, " Elapsed CPU Time: ",F10.3, " seconds")')
     +     finish-start

! determine and print Date and Time
      call date_and_time(date, time)
      write(43,'(" Date: ", A, ",", 1X, "Time: ", A, /)') date, time


      
      write (43,410) irem,sum124,bagct,ogc,xipht,icmoi,sumw1,kst1,icru
     +               st,iqm,j651,r0r0,r0r0km,bord0,bord1,iq2020,r5km,p
     +               rkht,j651,iccorr,einter,latent_heat,igravBary
  410 format (///,' irem=',i6,/,' accuracy of omega_bar(r): ',e12.6,/,
     +        ' bag constant = ',e12.6,'  MeV',/,' og(bar;center)= ',
     +        e12.6,'  1/sec',/,' interpolation technique: ',
     +        a25,/,' moment of inertia: icmoi=',i3,/,' accuracy of',
     +        ' w(1,hat;r): ',e12.6,/,' kst1=',i4,' icrust=',i3,' i',
     +        'qm=',i3,' j651=',i6,/,' r0r0=',e12.6,' (=> R_max=',f8.2,
     +        ' km)',/,' bord0=',e12.6 ,' bord1=',e12.6,/,' iq2020=',i6,
     +        '  r5km=',e12.6,' km',/,' ',a38,/,' j651=',i6,'  iccorr=',
     +        i3,/,' chosen value for inner crust density: ',e10.4,2x,
     +        'latent_heat=',i3,2x,'igravBary=',i2)

      write (43,483) iradot,iout20,iotomg,iout33,ioutes,iotmp0,iop30,
     +               iothv2,ioutpe,ioutse,etasc6,etasc7,igivmx,iscmx7,
     +               i4006a,i4006b,zedy(3),zedy(4),ahw,bhw,anv,bnv 
 483  format (' iradot=',i3,' iout20=',i3,' iotomg=',i3,' iout33=',i3,
     +        ' ioutes=',i3,' iotmp0=',i3,/,' iop30=',i3,' iothv2=',i3,
     +        ' ioutpe=',i3,' ioutse=',i3,/,' etasc6=',e9.2,1x,' etasc',
     +        '7=',e9.2,1x,' igivmx=',i3,' iscmx7=',i3,/,' i4006a=',i2,
     +        2x, 'i4006b=',i2,2x,' e_c(3)=',f7.2,' MeV/fm^3',2x, ' e_',
     +        '(4)=',f7.2,' MeV/fm^3',/,' Harrison-Wheeler boundaries:',
     +        ' a=',e12.6,2x,' b=',e12.6,' g/cm^3',/,' Negele-Vautheri',
     *        'n boundaries:',' a=',f10.7,2x,' b=', f10.6,' MeV/fm^3')

      if (edrip == (-1.)) then
      write (43,'(1x,''step size for radial integration: '',e12.6,'' m''
     +            ,//)') sslt
      else if (edrip /= (-1.)) then
      write(43,'(1x,''step size for  r < '',e12.6,'' km: '',f7.4,'' m'')
     +          ') r5km,sslt
      write(43,'(1x,''step size for  r > '',e12.6,'' km: '',f12.6,'' m''
     +            ,//)') r5km,ssgt
      end if

      write (43,802)
  802 format (/,2x,'  e_c: MeV/fm^3     g/ccm      log  R_gy/R_s R_s ',
     +        '(km)   n_c (1/fm^3)  E_int (MeV/fm^3)')
      write (43,888)
  888 format (' ---------------------------------------------------',
     +        '-----------------------------------')

      do 1000 m=1,i27
      if (tht1(m,5) == 0.) go to 1000
      write (43,800) m,tht1(m,1),tht1(m,2),tht1(m,3),tht1(m,4),tht1(m,5)
     +                ,tht1(m,25),tht1(m,26)
 1000 continue
  800 format (1x,i3,1x,e12.6,1x,e12.5,1x,f6.2,1x,f5.3,1x,e12.6,
     +        2x,e12.6,5x,f8.2)



      write (43,842)
 842  format (/,2x,'     O_K_FIP   P_K_FIP    f_K_FIP',1x,
     +             '    O_K_Newton  P_K_Newton     f_K_Newton')
      write (43,'(2x,''     1/sec      msec        Hz'', 8x,
     +               ''  1/sec        msec            Hz''/ , 17x,
     +                '' (using R=R_spherial star)'') ')
      write (43,888)

      do 1550 m=1,i27 
      if (tht1(m,5) == 0.) go to 1550
         
!  compute rotational period P_K_FIP (msec) using radius of non-rotating star
      P_K_FIP = zpi * 1.0E03 / tht1(m,6)

!  compute rotational frequency f_K_FIP (Hz) using radius of non-rotating star
      f_K_FIP = 1.0E03 /  P_K_FIP 

!  use OK_FIP = 2.4E05 sqrt(M/R**3) to compute OK_Newton = 3.63731E05 * sqrt(M/R**3)
!  (rad/sec) using radius of non-rotating star
      O_K_Newton = tht1(m,6) * 1.5155

!  compute rotational period P_K_FIP (msec) using radius of non-rotating star
      P_K_Newton = zpi * 1.0E03 / O_K_Newton
      
!  compute rotational frequency f_K_FIP (Hz) using radius of non-rotating star
      f_K_Newton = 1.0E03 / P_K_Newton
      
      write (43,808) m, tht1(m,6), P_K_FIP, f_K_FIP, O_K_Newton, 
     +                  P_K_Newton, f_K_Newton
 1550 continue
 808  format (1x,i3,1x,f8.2,2x,f9.3,3x,f8.3,7x,f7.2,3x,f9.3,6x,f9.3)



!  06/04/2015: output Kepler frequencies/periods using R_eq of rotating star
      write (43,843)
 843  format (/,2x,'     O_K_FIP   P_K_FIP    f_K_FIP',1x,
     +        '    O_K_Newton  P_K_Newton     f_K_Newton     f_Input')
      write (43,'(2x,''     1/sec      msec        Hz'', 8x,
     +               ''  1/sec        msec            Hz           Hz'',
     +            / ,17x, '' (using R=R_equator of rotating star)'') ')
      write (43,777)
 777  format (' ---------------------------------------------------',
     +        '-------------------------------------')      
      
      do  1560 m=1,i27 
      if ( tht1(m,5) == 0. .or. (tht1(m,16) < tht1(m,15)) ) go to 1560
      O_K_FIP = 2.4E05*sqrt(tht1(m,11)/tht1(m,16)**3)
!  compute rotational period P_K_FIP (msec) using radius of rotating star
      P_K_FIP = zpi * 1.0E03 / O_K_FIP

!  compute rotational frequency f_K_FIP (Hz) using radius of rotating star
      f_K_FIP = 1.0E03 /  P_K_FIP 

!  use OK_FIP = 2.4E04 sqrt(M/R**3) to compute OK_Newton = 3.63731E05 * sqrt(M/R**3)
!  (rad/sec) using radius of rotating star
      O_K_Newton = O_K_FIP * 1.5155

!  compute rotational period P_K_FIP (msec) using radius of rotating star
      P_K_Newton = zpi * 1.0E03 / O_K_Newton
      
!  compute rotational frequency f_K_FIP (Hz) using radius of rotating star
      f_K_Newton = 1.0E03 / P_K_Newton
      
      write (43,818) m, O_K_FIP, P_K_FIP, f_K_FIP, O_K_Newton, 
     +                  P_K_Newton, f_K_Newton, tabb(m,1)/zpi
 1560 continue
 818  format (1x,i3,1x,f8.2,2x,f9.3,3x,f8.3,7x,f7.2,3x,f9.3,6x,f9.3,3x,
     +        f9.3)
      

      write (43,902)
  902 format (//,6x,'og_s/og    og_c/og  del R/R_sp  M_sp   ',
     +        '  M_tot    dM/M_s  P_cut  (MeV/fm^3)')
      write (43,888)

      do 2000 m=1,i27
      if (tht1(m,7) == 0.) go to 2000
      write (43,900) m,tht1(m,7),tht1(m,8),tht1(m,9),tht1(m,10),tht1(m,1
     +                 1),tht1(m,12),thov1(m,12)
 2000 continue
  900 format (1x,i3,2x,e9.3,2x,e9.3,2x,f5.3,2x,f8.4,1x,f9.7,2x,f7.4,2x,
     +        e12.5)


      write (43,910)
  910 format (//,3x,' e_HT68    e_FIP86   R_pole    R_eq     Q/',
     +          '(Ms Rs^2)   log(I_negl)  log(I_exact)')
      write (43,888)

      do 3000 m=1,i27
      if (tht1(m,15) == 0.) go to 3000
      write (43,920) m,tht1(m,13),tht1(m,14),tht1(m,15),tht1(m,16),tht1
     +                (m,17),tht1(m,18),tht1(m,19)
 3000 continue
  920 format (1x,i3,1x,f5.2,1x,f9.3,1x,f9.3,1x,f9.3,3x,e9.3,6x,f9.6,3x,
     +        f9.6)


      write (43,930)
  930 format (//,1x,'dep. par.: e_FIP    R_p  R_eq;  Z_pole    Z_eq ',
     +           '  v+/c     v_/c        Omeg_K_ ')
      write (43,888)

      do 4000 m=1,i27
      if (tht1(m,21) == 0.) go to 4000
      tstrob(8,m)=tht1(m,23)
      write (43,940) m,tht1(m,20),tht1(m,21),tht1(m,22),tht1(m,23),tht1
     +                (m,24),tht2(m,8),tht2(m,25),tht2(m,26)
 4000 continue
  940 format (1x,i3,1x,f5.3,1x,f9.3,1x,f9.3,3x,f7.5,1x,f7.5,1x,f6.3,2x,
     +        f6.3,3x,f9.2)

      write (43,950)
  950 format (//,6x,'beta   T/M_sun     J/G M_rot^2  J/G M_sun^2 ',
     +          '    T/W      BE_rot/M_sun')
      write (43,888)

      do 5000 m=1,i27
      if (tht2(m,1) == 0.) go to 5000
      write (43,960) m,tht2(m,1),tht2(m,2),tht2(m,3),tht2(m,3)*tht1(m,11
     +               )**2,tht2(m,5),tht2(m,6)
 5000 continue
  960 format (1x,i3,1x,f5.3,1x,e12.6,1x,e12.6,1x,e12.6,1x,e12.6,1x,
     +        e12.6)


      write (43,9501)
 9501 format (//,6x,'BE_ov/M_sun',2x,' alph_1_ov ',2x,' alph_2_ov ',
     +           2x,' alph_1_rot',2x,' alph_2_rot')
      write (43,888)

      do 5011 m=1,i27
      if (tht2(m,1) == 0.) go to 5011
      write (43,966) m,ttby(m,9),ttby(m,10),ttby(m,11),ttby(m,12),
     +                 ttby(m,13)
 5011 continue
  966 format (1x,i3,2x,e11.4,2x,e11.4,2x,e11.4,2x,e11.4,2x,e11.4)


      write (43,1950)
 1950 format (//,6x,'R_sp        R(d_sp)  R_p       R(d_p)    R_eq  ',
     +          '     R(d_eq)      Q/km^3  ')   
      write (43,888)

      do 5001 m=1,i27
      if (tht2(m,1) == 0.) go to 5001
      qmoment=tht1(m,17)*(tht1(m,5)**2)*tht1(m,10)*rsunkm
      write (43,1960) m,tht1(m,5),tht2(m,21),tht1(m,15),tht2(m,19),
     +                  tht1(m,16),tht2(m,20),qmoment
 5001 continue
 1960 format (1x,i3,1x,f8.3,2x,f8.3,1x,f8.3,2x,f8.3,3x,f8.3,3x,f8.3,4x,
     +           f8.3) 


      write (43,2950)
 2950 format (//,6x,'M28_sph;d-s   M_sph;d-s/M_sun    M28_rot;d-s ',
     +          '  M_rot;d-s/M_sun ')
      write (43,888)

      do 5002 m=1,i27
      if (tht2(m,1) == 0.) go to 5002
      write (43,2960) m,tht2(m,15),tht2(m,18),tht2(m,17),tht2(m,16)
 5002 continue
 2960 format (1x,i3,2x,e11.5,4x,e11.5,7x,e11.5,5x,e11.5)


      write (43,3950)
 3950 format (//,6x,'log I_cr/g cm2    I_cr/I_tot    log I_tot/g ',
     +          'cm2    log I_H/g cm2  convergent')
      write (43,888)

      do 5062 m=1,i27
      if (tht2(m,1) == 0.) go to 5062
      convergent='yes'   
      if (tht2(m,13) > etasc7) convergent='no'
      write (43,1963) m,tht2(m,22),tht2(m,23),tht2(m,24),tht1(m,27),
     +                convergent
 5062 continue
 1963 format (1x,i3,2x,f10.6,7x,e11.5,5x,f10.6,10x,f10.6,5x,a3)


      write (43,951)
  951 format (//,2x,'  Og(K_GR_I)  v+/c  Og(K_GR_II)   v+/c ',
     +        '   Z_bw     Z_fw      d_%      v-/c')
      write (43,888)

      do 5005 m=1,i27
      if (tht1(m,7) == 0.) go to 5005
      write (43,961) m,tht2(m,7),tht2(m,8),tht2(m,9),tht2(m,10),
     +                 tht2(m,11),tht2(m,12),tht2(m,13),tht2(m,25)
 5005 continue
  961 format (1x,i3,2x,f8.2,1x,f5.3,3x,f8.2,4x,f5.3,2x,f7.4,2x,f7.4,
     +        2x,e9.2,1x,e10.3)


      write (43,751)
  751 format (//,1x,'del E_B/M_sun  lg(A_s)  lg(del A) ',
     +        'lg(A_tot)  del A/A_s  A_tot/E57  M_by/M_sun')
      write (43,888)

      do 5751 m=1,i27
      if (tht1(m,7) == 0.) go to 5751
      write (43,753) m,ttby(m,1),ttby(m,2),ttby(m,3),ttby(m,4),
     +                 ttby(m,5),ttby(m,6),ttby(m,7)
 5751 continue
  753 format (1x,i3,2x,e10.3,2x,f8.4,2x,f8.4,2x,f8.4,2x,f7.4,
     +        2x,f7.3,2x,f7.3)


      write (43,'(///,3x,''Properties of  n o n  rotating stars:'')')
      write (43,22)
   22 format (/,6x,'e_c',12x, 'e_c/e_0     M_gr   M_by   M_pr     R (k',
     +'m)   lg A_s   lg I/(g cm2)')
      write (43,888)

      do 23 m=1,i27
      if (tht1(m,7) == 0.) go to 23
      eumg=thov1(m,1)/engnm0
      write (43,24) m,thov1(m,1),eumg,tht1(m,10),ttby(m,8),thov1
     +                (m,4),tht1(m,5),thov1(m,5),thov1(m,6)
   23 continue
   24 format (1x,i3,1x,e12.6,2x,e12.6,1x,f6.3,1x,f6.3,1x,f6.3,1x,
     +        f9.3,2x,f9.6,1x,f9.6)


      write (43,'(///,3x,''Properties of  r o t a t i n g  stars:'')')
      write (43,7151)
 7151 format (/,6x,'e_c',12x, 'e_c/e_0     M_gr   M_by   M_pr     R (k',
     +'m)  lg A_rot  lg I/(g cm2)')
      write (43,888)

      do 5175 m=1,i27
      if (tht1(m,7) == 0.) go to 5175
      eumg=tht1(m,1)/engnm0
      tstrob(1,m)=tht1(m,1)*uf6/ee14
      tstrob(3,m)=tht1(m,11)
      tstrob(4,m)=ttby(m,7)-tht1(m,11)
      tstrob(5,m)=tht1(m,16)
      tstrob(7,m)=10**(tht1(m,19)-44.)
      write (43,7153) m,tht1(m,1),eumg,tht1(m,11),ttby(m,7),tht2(m,14),
     +                  tht1(m,16),ttby(m,4),tht1(m,19)
 5175 continue
 7153 format (1x,i3,1x,e12.6,2x,e12.6,1x,f6.3,1x,f6.3,1x,f6.3,1x,
     +        f9.3,2x,f9.6,1x,f9.6)


      write (43,65)
 65   format (///,6x,'t_quake_sp',2x,'t_quake_rot',4x,'Og(K_GR_II)',4x,
     +        'Omega')
      write (43,'(6x,'' (years)'',6x,''(years)'',9x, ''(Hz)'', 9x, 
     +           ''(Hz)'')')
      write (43,888)

      do m=1,i27
      if (tht1(m,7) /= 0.0) 
     +write (43,175) m,ttby(m,14),ttby(m,15),tht2(m,9)/zpi,tabb(m,1)/zpi
      end do
  175 format (1x,i3,2x,e10.4,3x,e10.4,5x,f8.2,4x,f8.2)

      write (43,15)
   15 format (///,6x,'   e_c    e_c/e_0   Om (rad/s)  R_eq ',
     +        ' M_rot/M_sun    B_min      B_max  ')
      write (43,'(58x,''  (1E+08  Gauss)'')')
      write (43,888)

      do 16 m=1,i27
      if (tht1(m,7) == 0.) go to 16
      eresc=tht1(m,1)/engnm0
      write (43,17) m,tht1(m,1),eresc,tabb(m,1),tht1(m,16),
     +                tht1(m,11),tabb(m,2),tabb(m,3)
   16 continue
   17 format (1x,i3,2x,f8.2,2x,f8.3,2x,f8.2,1x,f9.2,1x,f7.4,
     +        4x,f8.5,4x,e12.6)


      write (43,18)
   18 format (//,3x,'  Og_Macl     lg P(dot) lg P   P(dot)      P',
     +          '      T_s     Og_Macl  1/s')
      write (43,'(3x,''  e=0.81267    s s**-1   s    s s**-1    '',
     +          '' ms       yr     e=0.95289'')')
      write (43,888)

      do 20 m=1,i27
      if (tht1(m,7) == 0.) go to 20
      write (43,21) m,tabb(m,4),tabb(m,5),tabb(m,6),tabb(m,7),tabb(m,8),
     +                tabb(m,9),tabb(m,10)
   20 continue
   21 format (1x,i3,2x,f8.2,2x,f6.2,2x,f6.2,2x,e8.2,2x,f8.2,2x,e8.2,
     +        2x,f8.2)


      write (43,1818)
 1818 format (//,8x,'I',8x,'nue_eq',8x,'psi_eq',3x,'dnue/dr',3x,
     +          'dpsi/dr',6x,'f',9x,'g')
      write (43,'(7x,''(km^3)'',24x,''(1/km)    (1/km)'')')
      write (43,'(17x,''Metric as defined by Friedman-Ipser-Parker'',
     +                '' [FIP86]'')')
      write (43,888)

      do 1820 m=1,i27
      if (tht1(m,7) == 0.) go to 1820
      write (43,1821) m,tappx(m,9),tappx(m,1),tappx(m,2),tappx(m,3),
     +                  tappx(m,4),tappx(m,5),tappx(m,6)
 1820 continue
 1821 format (1x,i3,1x,e12.6,2x,f8.5,2x,f8.5,2x,f8.5,2x,f8.5,2x,f8.5,
     +        2x,f8.5)


      write (43,1828)
 1828 format (//,8x, 'f*g', 8x, 'O_K', 12x, 'Q(hat)', 10x, 'Q (km^3)')
      write (43,'(16x,''(rad/sec)'')')
      write (43,888)

      do 1830 m=1,i27
      if (tht1(m,7) == 0.) go to 1830
      write (43,1831) m, tappx(m,7), tappx(m,8), tquad(m,1), tquad(m,2)
 1830 continue
 1831 format (1x, i3, 2x, f8.5, 2x, f9.2, 6x, e12.6, 5x, e12.6)


      if (ivonr /= 0) then
      write (43,1848)
 1848 format (//,6x,'radial distance from origin',2x,
     +              'log_10 (I(r)/g cm^2)',2x,'log_10 n(r)',2x,
     +              '   e_c')
      write (43,'(14x,''(meters)'',36x,''(1/fm^3)'',5x,''(g/'',
     +                ''c^m3)'')')
      write (43,888)
      iadd=1
      do 1849 m=1,irem
      if (m == iadd .or. m == irem) then
      iadd=iadd+10
      write (43,1851) m,st(m),yout(m),pdent(m),feng1(m)*ec*uf6
      end if
 1849 continue
 1851 format (1x,i4,8x,f8.2,17x,e14.8,3x,e12.6,2x,e10.4)
      end if


      write (43,'(///,3x,''Table as computed in WFF88:'')')
      write (43,822)
  822 format (/,6x,'  e_c      P_c        M_G       M_A-M_G     R',
     +'     Delta_c    I        z')
      write (43,'(6x,''10**14'',3x,''10**34'',44x,''10**44'')')
      write (43,'(4x,''(g/cm**3)'',1x,''(dyn/cm**2)'',2x,''(M_sun)'',4x,
     +           ''(M_sun)'',4x,''(km)'',5x,''(km)'',2x,''(g cm**2)'')')
      write (43,888)

      do 823 m=1,i27
      if (tht1(m,7) == 0.) go to 823
      write (43,1879) m,tstrob(1,m),tstrob(2,m),tstrob(3,m),tstrob(4,m)
     +                 ,tstrob(5,m),tstrob(6,m),tstrob(7,m),tstrob(8,m)
 823  continue
 1879 format (1x,i3,1x,f6.3,2x,f9.4,4x,f6.3,4x,f8.4,2x,f7.2,1x,f7.2,1x,
     +        e12.6,1x,f7.4)

      
! 06/05/2015 -- For neutron stars with quark-hybrid composition only      
      if (ioutse /= 0) then
      write (43,'(///,3x,''Output mass and geometrical structure of '',
     +                   ''quark matter core'',/)')
      write (43,'(''    thresholds for '', a20, '' eos:'', /, 
     +            ''    mixed phase begins at '', f7.2, '' MeV/fm^3'',/, 
     +            ''    mixed phase ends at '', f7.2, '' MeV/fm^3'',  /, 
     +            ''    q-blobs    at '', f7.2, '' MeV/fm^3'', /, 
     +            ''    q-rods     at '', f7.2, '' MeV/fm^3'', /, 
     +            ''    q/h-slabs  at '', f7.2, '' MeV/fm^3'', /, 
     +            ''    h-rods     at '', f7.2, '' MeV/fm^3'', /, 
     +            ''    h-blobs    at '', f7.2, '' MeV/fm^3'', /, 
     +            ''    pure quark at '', f7.2, '' MeV/fm^3'')') 
     +   cthr,thresh(nthrcb,1),thresh(nthrcb,2),thresh_brs(nthrcb,1),
     +   thresh_brs(nthrcb,2),thresh_brs(nthrcb,3),thresh_brs(nthrcb,4),
     +   thresh_brs(nthrcb,5), thresh_brs(nthrcb,6)
      if (cthr == ('ioutse not specified          ,')) go to 7042
      
      write (43,854)
  854 format (/,6x,' e_c       M_qcore    M_mixed    M_rest  deviati',
     +             'on   Omega    M_total   Lattice')
      write (43,888)

      do 855 m=1,i27
      if (tht1(m,7) == 0.) go to 855

! 09/08/2018 -- needed for output table
         if(tht1(m,1)<thresh_brs(8,1)) blobs_rods_slabs='----------'
         if(thresh_brs(8,1)<=tht1(m,1) .and. tht1(m,1)<=thresh_brs(8,2))
     +       blobs_rods_slabs='q-blobs   '     
         if(thresh_brs(8,2)<tht1(m,1) .and. tht1(m,1)<=thresh_brs(8,3))
     +       blobs_rods_slabs='q-rods    '     
         if(thresh_brs(8,3)<tht1(m,1) .and. tht1(m,1)<=thresh_brs(8,4))
     +       blobs_rods_slabs='q/h-slabs '
         if(thresh_brs(8,4)<=tht1(m,1) .and. tht1(m,1)<=thresh_brs(8,5))
     +       blobs_rods_slabs='h-rods    '     
         if(thresh_brs(8,5)<tht1(m,1) .and. tht1(m,1)<=thresh_brs(8,6))
     +       blobs_rods_slabs='h-blobs   '     
         if (tht1(m,1)>thresh_brs(8,6)) blobs_rods_slabs='----------'
     +       
         

      write (43,1839) m,tht1(m,1),tabmm(1,m),tabmm(2,m),tabmm(3,m),
     +                  tabmm(4,m),tabb(m,1),tht1(m,11),blobs_rods_slabs
 855  continue
 1839 format (1x,i3,2x,f7.2,4x,f7.5,4x,f7.5,3x,f7.5,3x,e7.1,3x,f7.2,3x,  
     +        f6.4,3x,A10)

      write (43,864)
  864 format (/,6x,' e_c      M_qcore/M  M_mixed/M  M_rest/M   Omega ',
     +             '   M_total')
      write (43,888)

      
      do 865 m=1,i27
      if (tht1(m,7) == 0.) go to 865
      col1=tabmm(1,m)/tht1(m,11)
      col2=tabmm(2,m)/tht1(m,11)
      col3=tabmm(3,m)/tht1(m,11)
      write (43,1869) m,tht1(m,1),col1,col2,col3,tabb(m,1),tht1(m,11)
 865  continue
 1869 format (1x,i3,2x,f7.2,4x,f7.5,4x,f7.5,3x,f7.5,3x,f7.2,3x,f6.4)

      end if

 7042 continue
      close ( 43 )

      return

      end



!SEOS1                 ***       *****
      SUBROUTINE SEOS1(IEP,BAGCT,ICOMB,IQM,SC,IVARDM,RCQMAX,AHW,BHW,ANV,
     +                 BNV)


! ----------------------------------------------------------------------
!
!  PURPOSE:
!
!     COMPUTATION OF THE EQUATIONS OF STATE OF
!
!     1) HARRISON AND WHEELER (1965), (7.86<MASS DENSITY/(g cm**-3)<1E11),
!                                  (4.4E-12<ENERGY D./MeV fm**-3<5.6E-02),
!     2) NEGELE AND VAUTHERIN (1973), (1E11<     .  .   .   .      <1E13),
!                                     (5.6E-02<  .  .   .   .      <5.6 ).
!     3) HADRONIC MATTER  (LEPTONS AND BARYONS),
!
!     4) QUARK MATTER
!
!  RETURN: ENGD = ENERGY DENSITY (IN MeV/fm**3),
!          PRESS = PRESSURE (IN MeV/fm**3)
!
! ----------------------------------------------------------------------

      REAL * 8 EE18,EEM18,EE14,EE34,EE44
      real * 8 ee03,ee13,ee15
      
      REAL ERG(1),SC(4,IVARDM)

      PARAMETER (NP100=100,N7=7)
      REAL AUX1HW(NP100),AUX2HW(NP100),AUX3HW(NP100),AUX4HW(NP100)
      REAL AUX1NV(NP100),AUX2NV(NP100),AUX3NV(NP100),AUX4NV(NP100),
     +     AUX5NV(NP100)



      COMMON/PARNV1/CNVEOS(N7),A1NV,C0NV

      PARAMETER (N27=50,NPVAR=3*NP100)
      COMMON/ZZZZ/ZEDY(N27),ZMNS(N27),ZRNS(N27),ZTM(N27),ZKEPF(N27),
     +            ZFRP(N27),FIXMOV(N27)
      COMMON/ZZZ1/ZEDYOV(N27),ZEDYNW(N27),ZEDYHT(N27)
      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)

      PARAMETER (KP52=52)
      COMMON/CAXNR1/AXNR1(KP52),AXED1(KP52)

      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/ExpoMaxMin/EXMAX,EXMIN
      COMMON/B12345/B1,B2,B3,B4,B5
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/CTYPE/TYPEOS,TEOSCB

      CHARACTER * 50 TYPEOS
      CHARACTER * 45 TEOSCB
      CHARACTER * 11 CHEOD
      CHARACTER * 30 ARGINT
      CHARACTER * 30 EASS

      NAMELIST/NAME1/I100,I7,B1,B2,B3,B4,B5,UF1,C,L2,IVDY,IVAR
      NAMELIST/NAME2/I100,BAGC4,AUX1HW,AUX2HW,AUX3HW


      OPEN (UNIT=50, FILE='prns_eos.dat', status='old')



!__ PARAMETERS OF THE NV-EOS
!   GROUND-STATE CNVEOS(I):
      DATA (CNVEOS(I),I=1,7)/2.8822899E-01,5.9150523E-01,
     +                       9.0185940E-02,-1.1025614E-01,
     +                       2.9377479E-02,-3.2618465E-03,
     +                                     1.3543555E-04/
!   UNIFORM GAS C(I):
!     DATA (CNVEOS(I),I=1,7)/1.4821424,-4.0373482E-02,
!    +                       6.0455728E-02,-1.5307639E-02,
!    +                       3.4774416E-03,-4.3627154E-04,
!    +                       2.3383473E-05/
      DATA A1NV/1.E04/,C0NV/-4.0/
      DATA B1/3.0271E12/,B2/1.4415/,B3/1.50319E11/,B4/2.4616E15/
      DATA B5/2.3888E03/


      CHEOD='CHECK -EOF-'

      IF (IQM <= 0) THEN
      READ ( 50 , 515 ) TYPEOS
 515  FORMAT (A50)
      READ ( 50 , 516 ) L2
 516  FORMAT (I3)
         IF ((L2 > I100) .OR. (NPVAR /= (3*I100)) .OR.
     +       (IVAR /= IVARDM)) THEN
         IVDY=IVARDM
         WRITE (6,'(4X,'' _SEOS1_ ERROR STOP   '')')
         WRITE (6,NAME1)
         STOP
         END IF
      ELSE IF (IQM > 0) THEN
      TYPEOS='QUARK MATTER (BAG)                                '
      END IF


      I777=0




      IF (IQM == 0) THEN
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                        HARRISON-WHEELER EOS
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      A=AHW
      B=BHW
      AL=ALOG10(A)
      BL=ALOG10(B)
      NAB=I100
      ARGINT='SEOS1: LOG MASS DENSITY OF HW '
      CALL INTV(AL,BL,NAB,AUX1HW,0,ARGINT)
!__ RETURN: AUX1HW, i.e. LOG(MASS DENSITY, IN  g/cm**3)

      DO 324   LO=1,NAB
      AUX1HW(LO)=10**AUX1HW(LO)
  324 CONTINUE


!        NAB=I100
!        ARGINT='SEOS1: MASS DENSITY OF  HW EOS'
!        CALL INTV(A,B,NAB,AUX1HW,0,ARGINT)
!__ RETURN: AUX1HW (MASS DENSITY IN  g/cm**3)
      DO 10 I=1,NAB
      AUX2HW(I)=AUX1HW(I)/UF6
      EMEVFM=AUX2HW(I)
      CALL HWEOS(EMEVFM,PMEVFM)
!__ RETURN: PRESSURE  (MeV/fm**3)
      AUX3HW(I)=PMEVFM
      AUX4HW(I)=PMEVFM*UF4
   10 CONTINUE



! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                          NEGELE-VAUTHERIN EOS
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      A=ANV
!>>>>>
!     B=0.1E-02
!     B=5.0E-02
      B=BNV
      NAB=I100
      ARGINT='SEOS1: BARYON DENS.  OF NV EOS'
      CALL INTV(A,B,NAB,AUX1NV,0,ARGINT)
!__ RETURN: AUX1NV  (BARYON DENSITY IN  1/fm**3)
      DO 20 I=1,NAB
      RFM=AUX1NV(I)
      CALL NVEOS(AMNEUT,RFM,EMEVFM,PMEVFM)
!__ RETURN: ENERGY DENSITY  EMEVFM  (MeV/fm**3)
!           PRESSURE        PMEVFM  (MeV/fm**3)

!__ COMPARISON: ENERGY PER NUCLEON
!     EPN=EMEVFM/RFM-AMNEUT
!     D=RFM*1.E39
!     PRINT*,'E/N=',EPN,'  R=',RFM,'  = ',D,'  1/cm**3'

      AUX2NV(I)=EMEVFM
      AUX3NV(I)=PMEVFM
      AUX4NV(I)=PMEVFM*UF4
      AUX5NV(I)=EMEVFM*UF6
   20 CONTINUE




!__ COPY EQUATION OF STATE TO ARRAYS  ENGD, PRESS
      DO 30 J=1,3
      DO 30 I=1,I100
      IF (J == 1) THEN
      II=I
      ENGD(II)=AUX2HW(I)
      PRESS(II)=AUX3HW(I)
      ELSE IF (J == 2) THEN
      II=I+I100
      ENGD(II)=AUX2NV(I)
      PRESS(II)=AUX3NV(I)
      ELSE


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!          HADRONIC EQUATION OF STATE: LEPTONS AND BARYONS
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      IF (I > L2) GO TO 30
      II=I+2*I100
      READ ( 50 , 5177 ) ENGD(II),PRESS(II),PDEN(II)
 5177 FORMAT (6X,E14.8,4X,E14.8,4X,F6.4)
      IEP=I100+I100+L2
      END IF
   30 CONTINUE

      READ ( 50 , 518 ) CHEOD
 518  FORMAT (A11)
      if (CHEOD /= 'END OF DATA') then
      WRITE(6,'(/,''Last line in prns_eos.dat file not read'')')         
      STOP 'ERROR STOP in SEOS1'         
      end if

      if((typeos=='E_SH93I                                           ')
     + .or. 
     +  (typeos=='E_SH93II                                           '))
     + then
!  compute pressure from thermodynamic relation  P = n *d/dn (E/N)
!  note: eos in this case is given in input file in the form 
!        "k_F, E/N , 0"    

      DO 1402 K=1,L2
!  note: array "engd"  contains the Fermi momentum  
!        array "press" contains the energy per nucleon
      AUX1(K)= ENGD(I100+I100+K)**3 / (3 * piq)
 1402 AUX3(K)=PRESS(I100+I100+K)

      DO 1401 K=1,L2
      XXX=AUX1(K)
      EASS='ORG: CALL SELSPL   d/dn E/N   '
!                                    *****         RETURN:  d/dn (E/N)
      CALL SELSPL(AUX1,AUX3,L2,SC,10,DEDNN,XXX,20,EASS)
      AUX4(K)=DEDNN
 1401 CONTINUE

      DO 1403 K=1,L2
      ENGD(I100+I100+K) = (AUX3(K)+AMNEUT) * AUX1(K)
      PRESS(I100+I100+K)= AUX1(K)**2 * AUX4(K)
 1403 PDEN(I100+I100+K) = AUX1(K)
      end if


!__ INTERPOLATE BARYON DENSITY, n_e  (LOW-DENSITY REGION)
      KDO123=KP52
      L200  =2*I100
      EASS  ='SEOS1: CALL SELSPL 123        '

! 05/08/2015: remove following Do loop
!      DO 1233  I=1,L200
! 1233 AUX41(1,I)=AXNR1(I)
      DO 123 I=1,L200
      XXX     =ENGD(I)
!                                          ***         RETURN: YYY=N_E
      CALL SELSPL(AXED1,AXNR1,KDO123,SC,10,YYY,XXX,10,EASS)
      PDEN(I)=YYY
! NO ACCESS
!                                        ***         RETURN: ERG=n_P(N-1)
!      CALL SELINT(AXED1,AUX41,1,1,KDO123,ERG,XXX,1)
!      PDEN(I)=ERG(1)
!      print*,' diff=',abs(yyy-erg(1))
! END OF NO ACCESSD
  123 CONTINUE




!__ COMPUTE INTERNAL ENERGY
      DO 124  I=1,IEP
      EBET     =ENGD(I) - AMNEUT*PDEN(I)
      IF ((EBET < 0.) .AND. (EBET >= (-1.E-14))) EBET=-EBET
      EINTNL(I)=EBET
  124 CONTINUE


      DO 517 I=1,IEP
      PRINT*,' I=',I,'  E=',ENGD(I),'  P=',PRESS(I),'  N=',
     +PDEN(I)
  517 CONTINUE
      END IF



      IF (IQM < 0) THEN
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                    QUARK MATTER EQUATION OF STATE
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


      ICOMB=IEP+1
      IR1=2*I100
      IR1P1=IR1+1

!__ INTERPOLATE HYPERONIC EOS
 8888 CONTINUE
      AA=ENGD(IR1P1)
      BB=ENGD(IEP)
      ARGINT='SEOS2: ENERGY DENSITY QUARK   '
!                          ****
      CALL INTV(AA,BB,I100,AUX4,0,ARGINT)

      DO 40 I=1,L2
      AUX3(I)=ENGD(IR1+I)
      AUX31(1,I)=PRESS(IR1+I)
   40 CONTINUE

      IOPT=2
      DO 50 I=1,I100
      E=AUX4(I)
!                                            ***
!>>>>>     CALL INTER(IOPT,AUX3,AUX31,1,1,L2,ERG,E)
      CALL SELINT(AUX3,AUX31,1,1,L2,ERG,E,IOPT)
!     COPY INTERPOLATED HADRONIC EOS
      AUX1HW(I)=E
      AUX2HW(I)=ERG(1)
   50 CONTINUE

      IF (I777 /= 0) THEN
      ICOMB=10*I100
      GO TO 7777
      END IF


!__ COMPUTE QUARK MATTER EOS AND INTERPOLATE

      RBYAA=0.5
      RBYBB=10.0
!     (RANGE OF NUCLEAR DENSITY IN  1/FM**3)
!                                ****

      ISECC=0
 1111 ISECC=ISECC+1
      IF (ISECC > 100) THEN
      WRITE (6,'(//'' _SEOS1_ SECURITY COUNTER LIMIT EXCEEDED''//)')
      WRITE (6,'('' ISECC='',I4,'' RBYAA='',E13.6,'' RBYBB='',E13.6)')
     +       ISECC,RBYAA,RBYBB
      WRITE (6,'('' E(HYP;1)='',E13.6,''E(QUARK;1)='',E13.6,/,'' E(H''
     +      ,''YP;I100)='',E13.6,'' E(QUARK;I100)='',E13.6)') AUX1HW(1
     +       ),AUX3(1),AUX1HW(I100),AUX3(I100)
      STOP 'ERROR STOP'
      END IF

      ARGINT='SEOS1: BARYON DENSITY - QUARK '
      CALL INTV(RBYAA,RBYBB,I100,AUX4,0,ARGINT)

      ED3=1./3.
      BAGC4=BAGCT**4

      DO 60 I=1,I100
      RDY=AUX4(I)*UF13
      CHPL=(PIQ*RDY)**ED3
      HMEV4=3.*CHPL**4/(4.*PIQ) - BAGC4
      AUX31(1,I)=HMEV4/UF13
      AUX3(I)=(3.*HMEV4 + 4.*BAGC4)/UF13

!  TEST: ENERGY DENSITY RANGE OF THE CALCULATED QUARK MATTER EOS MUST COVER
!        THE ONE OF ORDINARY HYPERON MATTER (NECESSARY FOR LATER CALCULATIONS)
      IF ( (I == 1) .AND. (AUX3(I) > AUX1HW(I)) ) THEN
      RBYAA=RBYAA/5.
      GO TO 1111
      ELSE
     +IF ( (I == I100) .AND. (AUX3(I) < AUX1HW(I)) ) THEN
      RBYBB=RBYBB+.5
      GO TO 1111
      END IF

   60 CONTINUE

      IOPT=2
      DO 70 I=1,I100
      E=AUX1HW(I)
!                                     ***
      CALL SELINT(AUX3,AUX31,1,1,I100,ERG,E,IOPT)
!     COPY INTERPOLATED QUARK MATTER EOS
      AUX3HW(I)=ERG(1)
   70 CONTINUE


!__ INFORMATION SECTOR
      IF (AUX2HW(1)  >=  AUX3HW(1)) THEN
      WRITE (6,'(//,'' _SEOS1_ P(HADRONIC EOS;E(1)) > P(QUARK EOS;'',
     +       ''E(1)'',//)')
      ELSE
     +IF (AUX2HW(I100)  >=  AUX3HW(I100)) THEN
      WRITE (6,'(//,'' _SEOS1_ P(HADRONIC EOS;E(100)) > P(QUARK '',
     +       ''EOS;E(I100)'',//)')
      END IF
!__ END OF INFORMATION


!__ SEARCH FOR THE VALUE OF THE ENERGY DENSITY AT WHICH HYPERONIC AND
!   QUARK MATTER EOS WILL BE COMBINED

      D=AUX2HW(1)-AUX3HW(1)
      DO 80 I=2,I100
      ICOMB=I
      DD=AUX2HW(I)-AUX3HW(I)
      DSIGN=D*DD
      IF (DSIGN <= 0.) GO TO 9999
      IF (I == I100) THEN
      WRITE (6,'(//,'' _SEOS1_ NO CROSSING OF HADRONIC AND QUARK '',
     +       ''MATTER EOS FOUND'',//)')
      WRITE (6,NAME2)
      WRITE (6,'(//''USE ONLY HADRONIC MATTER EOS FOR THE CALCULATIONS''
     +       //)')
      I777=1
      GO TO 8888
      END IF
   80 CONTINUE


 9999 CONTINUE
!__ COPY EOS TO ARRAYS ENGD, PRESS;
!                                   1<=I< ICOMB  ><  HADRONIC EOS,
!                               ICOMB<=I<=I100   ><  QUARK MATTER EOS

 7777 DO 90 I=1,I100
      IRUN=IR1+I
      ENGD(IRUN)=AUX1HW(I)
      IF (I < ICOMB) THEN
      PRESS(IRUN)=AUX2HW(I)
      ELSE
      ENGD(IRUN)=AUX3HW(I)
      END IF
   90 CONTINUE

      IEP=3*I100
      END IF



      IF (IQM == 10) THEN
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                QUARK MATTER EQUATION OF STATE ONLY
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      ED3=1./3.
      ED4=1./4.
      BAGC4=BAGCT**4
      RADK=BAGC4*VPI*PI/3.
      CHMIN=RADK**ED4
      CHMIN=CHMIN*(1.-1./100.)

      RBYAA=CHMIN**3/PIQ
      RBYAA=RBYAA
      RBYBB=RCQMAX*UF13
      ARGINT='SEOS1: BARYON DENSITY - QUARK '
!                                ****     RANGE OF NUCLEAR DENSITY (MeV**3)
      CALL INTV(RBYAA,RBYBB,I100,AUX4,0,ARGINT)

      DO 316 I=1,I100
      CHPL=(PIQ*AUX4(I))**ED3
      HMEV4=3.*CHPL**4/(4.*PIQ) - BAGC4
!__ P
      AUX1HW(I)=HMEV4/UF13
!__ E
      AUX2HW(I)=(3.*HMEV4+4.*BAGC4)/UF13
  316 CONTINUE

!__ INTERPOLATE  EOS  AT SMALL PRESSURE VALUES
      STGG=(AUX1HW(2)-AUX1HW(1))/(AUX2HW(2)-AUX2HW(1))
      I42=I100+2
      ARGINT='SEOS1: ENERGY DENSITY - QUARK '
!                                       ****    ENERGY DENSITY
      CALL INTV(AUX2HW(1),AUX2HW(2),I42,AUX3,0,ARGINT)


      DO 416 I=1,I42
      AUX4(I)=STGG * (AUX3(I)-AUX3(1)) + AUX1HW(1)
  416 CONTINUE

      DO 417 I=1,2*I100
      IF (I <= I42) THEN
      ENGD(I) =AUX3(I)
      PRESS(I)=AUX4(I)
      ELSE
      ENGD(I) =AUX2HW(I-I100)
      PRESS(I)=AUX1HW(I-I100)
      END IF
      CHI4   = 4.*PIQ * (PRESS(I)*UF13 + BAGC4) / 3.
      CHI    = CHI4**ED4
      PDEN(I)= CHI**3 / (PIQ * UF13)
  417 CONTINUE

      IEP=2*I100
      ICOMB=1
      END IF





!__ CHECK ORDERING OF ENERGY DENSITY AND PRESSURE
!                              *****
 5555 CALL SORTOR(ENGD,IEP,-10,IERSO)
      IERSO1=IERSO
!                               *****
      CALL SORTOR(PRESS,IEP,-10,IERSO)
      IERSO2=IERSO

      IF (IERSO2 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _SEOS1_ PRESSURE NOT MONOTONOUSLY'',
     +       '' INCREASING'')')
      WRITE (6,505)
      STOP '_SEOS1_ pressure not monotonously increasing'
      END IF

      IF (IERSO1 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _SEOS1_ ENERGY DENSITY NOT MONOTONOUSLY'',
     +       '' INCREASING'',//)')
      DO 510 I=1,IEP
  510 WRITE (6,5155) I,ENGD(I),PRESS(I)
 5155 FORMAT (' I=',I4,2X,'E=',E13.6,2X,'P=',E13.6)
      STOP '_SEOS1_ energy density not monotonously increasing'
      END IF
  505 FORMAT (//,'**************************************************',
     +        '************************',////)

 
              CLOSE ( 50 )

      RETURN

      END


!SEOS2                 ***       *****           ***** ***** *****
      SUBROUTINE SEOS2(IEP,BAGCT,ICOMB,SC,IVARDM,EDRIP,PDRIP,IDRIA,
!                      ***** ******        ******
     +                 IDRIB,ITOTAL,ECQMAX,EQDRIP,EJGCM3,AHW,BHW,AN
     +                 V,BNV)


! ----------------------------------------------------------------------
!
!  PURPOSE:
!
!     COMPUTATION OF THE EQUATION OF STATE OF A QUARK STAR HAVING A
!     CRUST CONSISTING OF HADRONIC MATTER AT DENSITIES LOWER OR EQUAL
!     TO THE NEUTRON DRIP DENSITY, i.e.
!
!
!     1) CORE (DENSITY LARGER THAN RHO(DRIP)): QUARK MATTER
!
!     2) CRUST (DENSITY SMALLER THAN RHO(DRIP)):
!         A) HARRISON AND WHEELER (1965), (7.86<MASS DENSITY/(g cm**-3)<1E11),
!                                         (4.4E-12<EPS/MeV FM**-3<5.6E-02),
!         B) NEGELE AND VAUTHERIN (1973), (1E11<     .  .   .   .      <1E13),
!                                         (5.6E-02<  .  .   .   .      <5.6 )
!     OR
!
!         C) BAYM-PETHICK-SUTHERLAND
!
!
!  RETURN: ENGD = ENERGY DENSITY (IN MeV/fm**3),
!          PRESS = PRESSURE (IN MeV/fm**3)
!
! ----------------------------------------------------------------------

      REAL SC(4,IVARDM)
      REAL * 8 EE18,EEM18,EE14,EE34,EE44
      real * 8 ee03,ee13,ee15
      
      PARAMETER (NP100=100,N7=7)
      REAL AUX1HW(NP100),AUX2HW(NP100),AUX3HW(NP100),AUX4HW(NP100)
      REAL AUX1NV(NP100),AUX2NV(NP100),AUX3NV(NP100),AUX4NV(NP100),
     +     AUX5NV(NP100)

      COMMON/PARNV1/CNVEOS(N7),A1NV,C0NV

      PARAMETER (N27=50,NPVAR=3*NP100)
      COMMON/ZZZZ/ZEDY(N27),ZMNS(N27),ZRNS(N27),ZTM(N27),ZKEPF(N27),
     +            ZFRP(N27),FIXMOV(N27)
      COMMON/ZZZ1/ZEDYOV(N27),ZEDYNW(N27),ZEDYHT(N27)
      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)

      PARAMETER (KP52=52)
      COMMON/CAXNR1/AXNR1(KP52),AXED1(KP52)

      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/ExpoMaxMin/EXMAX,EXMIN
      COMMON/B12345/B1,B2,B3,B4,B5
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/CTYPE/TYPEOS,TEOSCB

      CHARACTER * 50 TYPEOS
      CHARACTER * 45 TEOSCB
      CHARACTER * 30 ARGINT
      CHARACTER * 30 EASS

      NAMELIST/NAME1/I100,I7,B1,B2,B3,B4,B5,UF1,C,IVDY,IVAR
      NAMELIST/NAME2/I100,BAGC4,AUX1HW,AUX2HW,AUX3HW


      IF ( (NPVAR /= (3*I100)) .OR. (IVAR /= IVARDM) ) THEN
      IVDY=IVARDM
      WRITE (6,'(4X,'' _SEOS2_ ERROR STOP  '')')
      WRITE (6,NAME1)
      STOP
      END IF

!__ CONSTANTS
      ED3=1./3.
      ED4=1./4.
      BAGC4 =BAGCT**4
      BAGCMF=BAGC4/UF13
      RADK  =BAGC4*VPI*PI/3.



      IF (ICRUST == 10) THEN
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                        HARRISON-WHEELER EOS
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      A=AHW
      B=BHW
      AL=ALOG10(A)
      BL=ALOG10(B)
      NAB=60
      ARGINT='SEOS2: MASS DENSITY - HW EOS  '
      CALL INTV(AL,BL,NAB,AUX1HW,0,ARGINT)
!__ RETURN: AUX1HW, i.e. LOG(MASS DENSITY, IN  g/cm**3)

      DO 324  LO=1,NAB
      AUX1HW(LO)=10**AUX1HW(LO)
  324 CONTINUE


      DO 10 I=1,NAB
      AUX2HW(I)=AUX1HW(I)/UF6
      EMEVFM   =AUX2HW(I)
      CALL HWEOS(EMEVFM,PMEVFM)
!__ RETURN: PRESSURE  (MeV/fm**3)
      AUX3HW(I)=PMEVFM
      AUX4HW(I)=PMEVFM*UF4
   10 CONTINUE
      IHW=NAB

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                          NEGELE-VAUTHERIN EOS
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!      A=4.0E-04; CORRESPONDS TO  E=0.375 MeV/fm**3
!      B=1.0E-03; CORRESPONDS TO  E=0.94  MeV/fm**3
      A=1.02E-04
      B=2.20E-04
      NAB=20
      ARGINT='SEOS2: BARYON DENSITY -NV EOS '
      CALL INTV(A,B,NAB,AUX1NV,0,ARGINT)
!__ RETURN: AUX1NV  (BARYON DENSITY IN  1/fm**3)
      DO 20 I=1,NAB
      RFM=AUX1NV(I)
      CALL NVEOS(AMNEUT,RFM,EMEVFM,PMEVFM)
!__ RETURN: ENERGY DENSITY  EMEVFM  (MeV/fm**3)
!           PRESSURE        PMEVFM  (MeV/fm**3)

!__ COMPARISON: ENERGY PER NUCLEON
!     EPN=EMEVFM/RFM-AMNEUT
!     D=RFM*1.E39
!     PRINT*,'E/N=',EPN,'  R=',RFM,'  = ',D,'  1/cm**3'

      AUX2NV(I)=EMEVFM
      AUX3NV(I)=PMEVFM
      AUX4NV(I)=PMEVFM*UF4
      AUX5NV(I)=EMEVFM*UF6
   20 CONTINUE
      INV=NAB


!__ COPY
      ITOTAL=IHW+INV
      IF (ITOTAL > I100) THEN
      WRITE(6,'(///,'' _SEOS2_  FATAL ERROR: ITOTAL > I100'',/ ,
     +          '' ITOTAL='',I4,'' I100='',I4,///)') ITOTAL,I100
      STOP
      END IF
      DO 2910 I=1,ITOTAL
      IF (I <= IHW) THEN
      ENGD(I)=AUX2HW(I)
      PRESS(I)=AUX3HW(I)
      ELSE
      ENGD(I)=AUX2NV(I-IHW)
      PRESS(I)=AUX3NV(I-IHW)
      END IF
 2910 CONTINUE


!__ INTERPOLATE BARYON DENSITY  N(EPS)  (LOW-DENSITY REGION)
      KDO123=KP52
      EASS='SEOS2: CALL SELSPL 123        '
      DO 123 I=1,ITOTAL
      XXX=ENGD(I)
!                                          ***         RETURN: YYY=N_E
      CALL SELSPL(AXED1,AXNR1,KDO123,SC,10,YYY,XXX,10,EASS)
      PDEN(I)=YYY
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      SET  PARTICLE DENSITY  N(R)=0 (NOT NEEDED)
       PDEN(I)=0.
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  123 CONTINUE



      ELSE IF (ICRUST == 20) THEN
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!               BAYM-PETHICK-SUTHERLAND EQUATION OF STATE
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



!   COMPUTE FEYNMAN-METROPOLIS-TELLER EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL FMT49(AUX3,AUX4,AUX1,IDIM,IFMT)
      DO 1010 I=1,IFMT
      ENGD(I) =AUX3(I)
      PDEN(I) =AUX1(I)
 1010 PRESS(I)=AUX4(I)

!   COMPUTE BAYM-PETHICK-SUTHERLAND EOS
      IDIM=IVAR
!                **** **** ****      ******
      CALL BPSCR(AUX3,AUX4,AUX1,IDIM,IBPSCR)
      ITOTAL=IFMT+IBPSCR
      IF (IBPSCR > I100) THEN
!__ ERROR STOP
      WRITE (6,'(//,'' DIMENSION OF AUXILIARY ARRAYS AUX3 AND AUX4 IS'',
     +          '' TOO SMALL!'',/,'' ITOTAL='',I4,'' IDIM='',I4,'' IB'',
     +          ''PSCR='',I4,'' IVAR='',I4,'' I100='',I4)') ITOTAL,IDIM,
     +      IBPSCR,IVAR,I100
      STOP
      END IF

      DO 200 I=1,IBPSCR
      ENGD(I+IFMT) =AUX3(I)
      PDEN(I+IFMT) =AUX1(I)
  200 PRESS(I+IFMT)=AUX4(I)



      ELSE



      WRITE(6,'(///,'' _SEOS2_ FATAL ERROR: ICRUST NOT CORRECTLY'',
     +          '' INITIALIZED'',/,'' ICRUST='',I4,///)') ICRUST
      STOP

      END IF



      NZERO=10
      IF(PDEN(2) == 0.) NZERO=0


      IF (EJGCM3 /= 0.) THEN
!   CHECK: Is the inner density of the hadronic crust (input) equal or
!          below the neutron drip density? In the latter case PRESS, ENGD,
!          and PDEN, each of length ITOTAL, must be modified.

!   check whether  EJGCM3=e_join (in g/cm**3)  is properly specified
      EJMFM3 = EJGCM3/UF6
      L10    = 10

      IF ( (EJMFM3 < ENGD(L10)) .OR. (EJMFM3 >= ENGD(ITOTAL)) )  THEN
      WRITE (6,'(///,'' e_join  out of range'',/,'' hadronic matter'',
     +               '' eos cannot be joined'',//,'' EJGCM3='',E12.6,
     +               '' g/ccm;  EJMFM3='',E12.6,'' MeV/fm**3'',//,
     +               '' ENGD_L10='',E12.6,''; ENGD_ITOTAL='',E12.6,
     +               '' MeV/fm**3; ITOTAL='',I4,''; L10='',I3,//)')
     +      EJGCM3,EJMFM3,ENGD(L10),ENGD(ITOTAL),ITOTAL,L10
      STOP 'e_join out of range _SEOS2_'
      END IF


      V = EJMFM3
      DO 600 I=1,ITOTAL-1
      II = I
      IF ( (ENGD(I) <= V) .AND. (ENGD(I+1) >= V) ) THEN
      ITOTNW = II+1
      IF (ENGD(I) == V) ITOTNW = II
      GO TO 610
      END IF
  600 CONTINUE
  610           CONTINUE


!  COPY EOS
      DO 620 I=1,ITOTAL
      AUX1(I) = ENGD(I)
      AUX3(I) = PRESS(I)
      AUX4(I) = PDEN(I)
  620 CONTINUE


!  COMPUTE PRESSURE  P_join(e_join)
      EASS='SEOS2: CALL SELSPL P_join     '
      XXX=V
!                                        ***           RETURN: YYY=P_e
      CALL SELSPL(AUX1,AUX3,ITOTAL,SC,10,YYY,XXX,10,EASS)
      PRESS(ITOTNW) = YYY

!  COMPUTE PARTICLE DENSITY  n_join(e_join)
      IF(NZERO /= 0) THEN
      EASS='SEOS2: CALL SELSPL n_join     '
      XXX=V
!                                        ***           RETURN: YYY=n_e
      CALL SELSPL(AUX1,AUX4,ITOTAL,SC,10,YYY,XXX,10,EASS)
      ELSE
      YYY=0.
      END IF
      PDEN(ITOTNW) = YYY

      ENGD(ITOTNW) = V


!  COPY EOS,  now: 0 < e   <= e_join
!                  0 < P_e <= P_join
!                  0 < n_e <= n_join

      DO 630  I=1,ITOTNW-1
      ENGD(I)  = AUX1(I)
      PRESS(I) = AUX3(I)
      PDEN(I)  = AUX4(I)
      IF (NZERO == 0) PDEN(I) = 0.
  630 CONTINUE


      ITOTAL = ITOTNW

      END IF
!     ^^^^^^-------------------------- end of  e_join # e_drip



      ITOP1=ITOTAL+1
      IDRIA=ITOTAL


!__ DRIP PRESSURE  P_drip, MeV/fm**3
      PDRIP =PRESS(ITOTAL)
      PDRIPU=PDRIP*UF13
!__ E(hadronic phase,P_drip), MeV/fm**3
      EDRIP = ENGD(ITOTAL)
!__ E(quark,P_drip), MeV/fm**3
      EQDRIP=3.*PDRIP + 4.*BAGCMF
!__ N(hadronic,P_drip), 1/fm**3
      RNDRIP= PDEN(ITOTAL)



      AA=EDRIP  * (1.+1.E-02)
      BB=EQDRIP * (1.-1.E-02)
      IF ((BB-AA) <= 0.) THEN
      WRITE (6,'(////,'' e(quark phase at P=P_drip) <'',
     +                '' e(hadronic phase at P=P_drip)'',//,
     +                '' choose a larger bag constant; for this'',
     +                '' calculation, B='',F8.3,'' MeV'',////)')
     +      BAGCT
      STOP 'BAG CONSTANT TOO SMALL _SEOS2_'
      END IF



      I57=I100-ITOTAL
      ARGINT='SEOS2: ENERGY DENSITY - BPS   '
!                         ****   COMPUTE ENERGY DENSITY GRID (MeV/fm**3)
      CALL INTV(AA,BB,I57,AUX4,0,ARGINT)

      DO 4116 I=1,I57
!__ E, MeV/fm**3
      ENGD(I+ITOTAL) =AUX4(I)
!__ P, MeV/fm**3
      PRESS(I+ITOTAL)=PDRIP
!__ BARYON DENSITY: SEE BELOW
 4116 CONTINUE

      ITOTRE = ITOTAL
      ITOTAL = ITOTAL + I57
      ITOP1  = ITOTAL + 1
      IDRIB  = ITOP1
      ICOMB  = IDRIB
      IF (ITOTAL /= I100) THEN
      WRITE (6,'(////,'' _SEOS2_ ITOTAL # I100'',/,'' ITOTAL='',
     +          I4,/,'' I100='',I4,////)') ITOTAL,I100
      STOP
      END IF




! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                   QUARK MATTER EQUATION OF STATE
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


!__ MU(P_drip), MeV
      CHDRIP=(VPI*PI*(PDRIPU+BAGC4)/3.)**ED4
      CHMIN =CHDRIP
!__ N(quark,P_drip), 1/fm**3
      RNQDRI=CHMIN**3 / (PIQ*UF13)


      AA=RNDRIP * (1.+1.E-02)
      BB=RNQDRI * (1.-1.E-02)
      ARGINT='SEOS2: ENERGY DENSITY - QUARK '
!                         ****   COMPUTE BARYON DENSITY GRID (1/fm**3)
      CALL INTV(AA,BB,I57,AUX4,0,ARGINT)

      DO 4216 I=1,I57
!__ N, 1/fm**3
      PDEN(I+ITOTRE) =AUX4(I)
!>>>>      PDEN(I+ITOTRE) =0.
 4216 CONTINUE


      AA=EQDRIP
      BB=ECQMAX
      ARGINT='SEOS2: ENERGY DENSITY; DP-ECMX'

!                          ****   COMPUTE ENERGY DENSITY GRID  (MeV/fm^3)
      CALL INTV(AA,BB,I100,AUX4,0,ARGINT)

      DO 316 I=1,I100
!__ E, MeV/fm**3
      AUX2HW(I)=AUX4(I)
!__ P, MeV/fm**3
      AUX1HW(I)=(AUX2HW(I) - 4.*BAGCMF) / 3.
  316 CONTINUE


!>> CHECK ACCURACY
      DZ=(AUX2HW(1)-EQDRIP)/EQDRIP
      DZ=ABS(DZ)
      EPSDZ=1.E-06
      IF (DZ > EPSDZ) THEN
      WRITE (6,'(////,'' WARNING (_SEOS2_): E(QUARK;P=P_drip) '',
     +           ''DEVIATES FROM ITS NOMINAL VALUE'',/,'' E(QUAR'',
     +           ''K;P=P_drip)='',E14.8,'' (NOMINAL VALUE)'',/,1X ,
     +           ''E(QUARK;P=P_drip)='',E14.8,'' (CALCULATED!)'',/,
     +           '' DIFFERENCE: '',E12.6,//,'' CONTINUE CALCULAT'',
     +           ''ION'',////)') EQDRIP,AUX2HW(1),DZ
      END IF
!>> END OF ACCURACY CHECK




!__ INTERPOLATE  BAG EOS  AT PRESSURE VALUES SLIGHTLY LARGER THAN P_join

      AA=AUX2HW(1)
      BB=AUX2HW(2)
      ARGINT='SEOS2: ENERGY DENSITY - BAG I '
!                          ****          ENERGY DENSITY  (MeV/fm^3)
      CALL INTV(AA,BB,I100,AUX4,0,ARGINT)

      DO 3116 I=1,I100
!__ E, MeV/fm**3
      AUX4HW(I)=AUX4(I)
!__ P, MeV/fm**3
      AUX3HW(I)=(AUX4HW(I) - 4.*BAGCMF) / 3.
      
!cc  BEGIN OF CHANGE
      IF(AUX3HW(I) < PDRIP) AUX3HW(I)=PDRIP
!cc  END OF CHANGE
!__ N, 1/fm**3
      CHMEV4   =(AUX3HW(I)*UF13 + BAGC4)*VPI*PI/3.
      CHMEV    =CHMEV4**ED4
      CHFM     =CHMEV/UF1
      AUX5NV(I)=CHFM**3 / PIQ
 3116 CONTINUE


      AA=AUX2HW(3)
      BB=AUX2HW(I100)
      ARGINT='SEOS2: ENERGY DENSITY - BAG II'
!                          ****          ENERGY DENSITY  (MeV/fm^3)
      CALL INTV(AA,BB,I100,AUX4,0,ARGINT)

      DO 3117 I=1,I100
!__ E, MeV/fm**3
      AUX4NV(I)=AUX4(I)
!__ P, MeV/fm**3
      AUX3NV(I)=(AUX4NV(I) - 4.*BAGCMF) / 3.
!__ N, 1/fm**3
      CHMEV4   =(AUX3NV(I)*UF13 + BAGC4)*VPI*PI/3.
      CHMEV    =CHMEV4**ED4
      CHFM     =CHMEV/UF1
      AUX1NV(I)=CHFM**3 / PIQ
 3117 CONTINUE




!__ C O P Y   QUARK MATTER EOS TO ARRAYS ENGD, PRESS, PDEN:
!   NOTE:       1 <= I <  ITOP1    ><   HADRONIC MATTER EOS (CRUST)
!   NOW:    IDRIB <= I <= 3*I100   ><   QUARK MATTER EOS

      I143 = ITOP1+I100-1
      I144 = I143+1
      I243 = I144+I100-1
      IEP  = I243

      DO 90 I=ITOP1,I243
      IF (I >= ITOP1. AND. I <= I143) THEN
      ENGD(I)  = AUX4HW(I-ITOTAL)
      PRESS(I) = AUX3HW(I-ITOTAL)
      PDEN(I)  = AUX5NV(I-ITOTAL)
      ELSE
      ENGD(I)  = AUX4NV(I-I143)
      PRESS(I) = AUX3NV(I-I143)
      PDEN(I)  = AUX1NV(I-I143)
      END IF
   90 CONTINUE




!__ CHECK ORDERING OF ENERGY DENSITY AND PRESSURE
!                              *****
 5555 CALL SORTOR(ENGD,IEP,-10,IERSO)
      IERSO1=IERSO
!                               *****
      CALL SORTOR(PRESS,IEP,-10,IERSO)
      IERSO2=IERSO

      IF (IERSO2 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _SEOS2_ PRESSURE NOT MONOTONOUSLY'',
     +       '' INCREASING'')')
      WRITE (6,505)
      END IF

      IF (IERSO1 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _SEOS2_ ENERGY DENSITY NOT MONOTONOUSLY'',
     +       '' INCREASING'',//)')
      DO 510 I=1,IEP
  510 WRITE (6,515) I,ENGD(I),PRESS(I)
  515 FORMAT (' I=',I4,2X,'E=',E13.6,2X,'P=',E13.6)
      STOP 'ERROR STOP'
      END IF
  505 FORMAT (//,'**************************************************',
     +        '************************',//)


      RETURN

      END


!SEOSWhiteDwarf                 ***  Return  /ep1/engd,press,pden,e_intnl
      subroutine SEOSWhiteDwarf(IEP)


! ----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Input the EOS of white dwarf matter from an external data file
!
!  RETURN: ENGD  = energy density (in MeV/fm**3),
!          PRESS = pressure (in MeV/fm**3),
!          PDEN  = particle number density (in 1/fm**3)      
!
! ----------------------------------------------------------------------

      PARAMETER (NP100=100,NPVAR=3*NP100)
      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)

      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/CTYPE/TYPEOS,TEOSCB

      CHARACTER * 50 TYPEOS
      CHARACTER * 45 TEOSCB
      CHARACTER * 11 CHEOD

      OPEN (UNIT=50, FILE='prns_eos.dat', status='old')


!__ Read eos data from external data file
!   Check header
      READ ( 50 , 515 ) TYPEOS
 515  FORMAT (A50)
      READ ( 50 , 516 ) L2
 516  FORMAT (I3)
         IF ( L2 > (3*I100) ) THEN
         WRITE (6,'(///,4X,'' _SEOSWhiteDwarf_  ERROR STOP '')')
         STOP 'ERROR STOP in _SEOSWhiteDwarf_'
         END IF

!   Read data         
      DO I = 1, L2
      READ ( 50 , 5177 ) ENGD(I), PRESS(I), PDEN(I)
      END DO
 5177 FORMAT (6X,E14.8,4X,E14.8,4X,E12.6)      

!   Check trailer      
      READ ( 50 , 518 ) CHEOD
 518  FORMAT (A11)
      if (CHEOD /= 'END OF DATA') then
      WRITE(6,'(/,''Last line in prns_eos.dat file not read'')')
      STOP 'ERROR STOP in _SEOSWhiteDwarf_'
      end if

! Number of data       
      IEP = L2

      return
      end subroutine SEOSWhiteDwarf

      

!NVEOS                            **** *****
      SUBROUTINE NVEOS(AMNEUT,BYD,ENGD,PRESS)


! ----------------------------------------------------------------------
!
! PURPOSE:
!     CALCULATION OF THE EQUATION OF STATE OF NEGELE AND VAUTHERIN
!
!     INPUT: BARYON DENSITY        BYD (IN 1/fm**3)
!                            1E-04<BYD<1E-01  1/fm**3
!
!     RETURN: ENERGY DENSIY   ENGD (IN MeV/fm**3)
!             PRESSURE        PRESS (IN MeV/fm**3)
!
! ----------------------------------------------------------------------

      PARAMETER (NP100=100,N7=7)
      COMMON/PARNV1/CNVEOS(N7),A1NV,C0NV

      NAMELIST/NAME1/RO,N7NL,A1NL
      N7NL = N7
      A1NL = A1NV


!__ BARYON DENSITY IN  1/fm**3
      RO=BYD
      IF (RO < 1.E-04 .OR. RO > 1.E-01) THEN
      WRITE (6,'(4X,'' ++ _NVEOS_ BARYON DENSITY OUT OF RANGE ++'')')
      WRITE (6,NAME1)
      STOP
      END IF
      SU=0.
      SU2=0.

      DO 2 J=1,N7
      X=ALOG(RO*A1NV)
      IM1=J-1
      SU=SU+CNVEOS(J)*X**IM1
      IM2=J-2
      IF (J > 1) SU2=SU2+CNVEOS(J)*FLOAT(IM1)*X**IM2
    2 CONTINUE

!__ ENERGY DENSITY IN  MeV/fm**3
      EPS=RO*(AMNEUT+C0NV+EXP(SU))
      ENGD=EPS


!__ PRESSURE IN  MeV/fm**3
      DEDRO=EXP(SU)*SU2/RO
      PRESS=RO*RO*DEDRO


      RETURN

      END

!HWEOS                      *****
      SUBROUTINE HWEOS(ENGD,PRESS)


! ----------------------------------------------------------------------
!
! PURPOSE:
!     CALCULATION OF THE EQUATION OF STATE OF HARRISON AND WHEELER
!
!     INPUT:  ENERGY DENSITY  ENGD (IN MeV/fm**3)
!     RETURN: PRESSURE        PRESS (IN MeV/fm**3)
!
! ----------------------------------------------------------------------


      COMMON/B12345/B1,B2,B3,B4,B5
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10

      NAMELIST/NAME1/RGCM3,R,B1,B2,B3,B4,B5,E1,E2,E3,E4,I


      E1=1./3.
      E2=5./4.
      E3=-.5
      E4=-5./6.



!__ MASS DENSITY  IN  g/cm**3
      RGCM3=ENGD*UF6
      R=RGCM3

      IF (R <= 1.E06) THEN
!__ REGIME  I: 7.86 < MASS DENSITY < 1 E06 g/cm**3
!              4.41 E-12< ..... <5.61 E-02 MeV/fm**3

!__ PRESSURE IN g/(cm*s**2)=dyn/cm**2
      PGCMS2=B1*(R**E1-B2)**5-B3

      ELSE IF (R > 1.E06 .AND. R <= 1.E11) THEN

!__ REGIME  II: 1 E06 < MASS DENSITY < 3.22 E11 g/cm**3
!               5.6 E-07<  .....  <0.181    MeV/fm**3

      PGCMS2=B4*(R**E2)*((1.+B5*R**E3)**E4)

      ELSE
      WRITE (6,'(4X,'' ++ _HWEOS_ MASS DENSITY OUT OF RANGE ++''
     +         )')
      WRITE (6,NAME1)
      STOP
      END IF

!__ PRESSURE: dyn/cm**2 --> MeV/fm**3
      PMEVFM=PGCMS2/UF4
      PRESS=PMEVFM


      RETURN

      END

!HWHT68                 *** *****
      SUBROUTINE HWHT68(IEP,ICOMB)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATE THE EOS USED BY HARTLE AND THORNE IN [HT68]
!
!
! RETURN: EOS, AVAILABLE VIA  COMMON/EP1/
!
!-----------------------------------------------------------------------


      DOUBLE PRECISION GUF,P,E,EI 


      PARAMETER (NP100=100,NPVAR=3*NP100)

      DIMENSION P(NPVAR),E(NPVAR),EI(NPVAR),YYY(1)

      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST

      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT

      CHARACTER * 30 ARGINT


      DATA GUF/0.74214976E-28/
!     X (1/cm**2) / GUF = X (g/cm^3),    [GUF] = cm/g



!__ PRESSURE   [in  1/cm**2]
      DATA (P(I), I=1,36)/
     +8.31D-41,4.17D-40,8.31D-40,4.17D-39,8.31D-39,2.79D-38,2.38D-37,
     +1.37D-36,7.00D-36,6.96D-35,4.79D-34,1.74D-33,5.95D-33,1.56D-32,
     +4.62D-32,2.67D-31,9.63D-31,4.83D-30,2.32D-29,5.19D-29,1.65D-28,
     +8.23D-28,2.37D-27,7.19D-27,2.69D-26,7.86D-26,1.93D-25,6.60D-25,
     +1.65D-24,4.18D-24,1.35D-23,3.29D-23,8.07D-23,2.67D-22,6.53D-22,
     +1.21D-21/
      DATA (P(I), I=37,72)/
     +2.73D-21,6.49D-21,1.10D-20,1.88D-20,3.05D-20,4.58D-20,6.59D-20,
     +9.55D-20,1.50D-19,2.54D-19,4.49D-19,9.14D-19,1.88D-18,6.09D-18,
     +2.63D-17,8.23D-17,2.60D-16,1.09D-15,3.25D-15,9.71D-15,3.93D-14,
     +9.71D-14,2.42D-13,7.34D-13,1.65D-12,3.60D-12,1.01D-11,2.08D-11,
     +4.28D-11,1.10D-10,2.26D-10,4.64D-10,1.19D-09,2.43D-09,4.91D-09,
     +1.23D-08/

!__ ENERGY   [in  1/cm**2]
      DATA (E(I), I=1,36)/
     +5.82D-28,5.84D-28,5.86D-28,5.97D-28,6.04D-28,6.32D-28,8.52D-28,
     +1.23D-27,2.34D-27,6.54D-27,1.56D-26,2.95D-26,5.22D-26,8.52D-26,
     +1.53D-25,3.71D-25,7.41D-25,1.86D-24,4.68D-24,7.41D-24,1.48D-23,
     +3.71D-23,7.41D-23,1.48D-22,3.71D-22,7.41D-22,1.48D-21,3.71D-21,
     +7.41D-21,1.48D-20,3.71D-20,7.41D-20,1.48D-19,3.71D-19,7.41D-19,
     +1.17D-18/
      DATA (E(I), I=37,72)/
     +2.34D-18,4.68D-18,7.41D-18,1.17D-17,1.86D-17,2.95D-17,4.68D-17,
     +7.41D-17,1.17D-16,1.86D-16,2.95D-16,4.68D-16,7.41D-16,1.48D-15,
     +3.71D-15,7.41D-15,1.48D-14,3.71D-14,7.41D-14,1.48D-13,3.71D-13,
     +7.41D-13,1.48D-12,3.71D-12,7.41D-12,1.48D-11,3.71D-11,7.41D-11,
     +1.48D-10,3.71D-10,7.41D-10,1.48D-09,3.71D-09,7.41D-09,1.48D-08,
     +3.71D-08/

!__ INTERNAL ENERGY   [in  1/cm**2]
      DATA (EI(I), I=1,36)/
     +1.00D-45,7.11D-43,2.77D-42,4.16D-41,1.12D-40,8.64D-40,3.36D-38,
     +3.21D-37,3.47D-36,5.01D-35,3.76D-34,1.52D-33,5.18D-33,1.45D-32,
     +4.73D-32,2.72D-31,1.05D-30,5.82D-30,3.03D-29,6.82D-29,2.27D-28,
     +1.11D-27,3.59D-27,1.12D-26,4.84D-26,1.42D-25,4.03D-25,1.53D-24,
     +4.07D-24,1.07D-23,3.77D-23,9.57D-23,2.41D-22,8.17D-22,2.04D-21,
     +3.72D-21/
      DATA (EI(I), I=37,72)/
     +9.21D-21,2.25D-20,4.05D-20,7.21D-20,1.28D-19,2.25D-19,3.88D-19,
     +6.60D-19,1.11D-18,1.88D-18,3.17D-18,5.39D-18,9.28D-18,2.18D-17,
     +7.27D-17,1.90D-16,5.16D-16,2.02D-15,5.70D-15,1.61D-14,6.30D-14,
     +1.69D-13,4.34D-13,1.43D-12,3.37D-12,7.71D-12,2.23D-11,4.87D-11,
     +1.04D-10,2.82D-10,5.89D-10,1.22D-09,3.19D-09,6.52D-09,1.33D-08,
     +3.41D-08/



      I70=72
      IV =IVAR

      IF (IV < I70) THEN
      WRITE (6,'(////,'' _HWHT68_ FATAL ERROR DETECTED'',/,'' IV <'',
     +       '' I70:  IV='',I3,'' I70='',I3)') IV,I70
      STOP
      END IF


      IWT=0
      SCALE=GUF*UF6

      DO 1000 I=1,I70

!  X (1/cm^2) / SCALE = X (MeV/fm^3)
       ENGD(I)   = SNGL( E(I) / SCALE )
       PRESS(I)  = SNGL( P(I) / SCALE )
       EINTNL(I) = SNGL( EI(I) / SCALE )

      AUX3(I)   = ENGD(I)
      AUX31(1,I)= PRESS(I)

!  DENSITY OF BARYONS PER  1/fm**3
      PDEN(I)   = (ENGD(I)-EINTNL(I)) / AMNEUT

      IF (IWT /= 0) THEN
      IF (I == 1) WRITE (6,332)
      WRITE (6,333) I,PRESS(I),ENGD(I),EINTNL(I),PDEN(I)
      END IF
 1000 CONTINUE

  333 FORMAT (1X,I3,2X,E9.3,2X,E9.3,2X,E9.3,2X,E9.3)
  332 FORMAT (/,' HARRISON-WHEELER EOS AS GIVEN BY [HT68]',/)
      IEP=I70
      ICOMB=-1



!>>>> INTERPOLATE   EOS
      IIPOL=0
!....................   CHOOSE  KI < 29
      KI29=29
!....................
      KI29P1=KI29+1
      IF (IIPOL /= 0) THEN

         XBORDA=AUX3(KI29P1)
         XBORDB=AUX3(I70)
         INAB=IVAR-KI29
         IF (INAB < (2*KI29)) THEN
         WRITE (6,'(///,'' NOTE: INTERPOLATION IN SUBROUTINE  HWHT68'',
     +   '' DOES NOT MAKE'',/,'' SENSE FOR TOO SMALL VALUES OF  IVAR.'',
     +   /,'' CHOOSE LARGER VALUE ( < 300 ); PRESENT VALUE: IVAR='',I3)'
     +   ) IVAR
         END IF
      ARGINT='HWHT68: ENERGY DENSITY        '
!                                     ****    ENERGY DENSITY
         CALL INTV(XBORDA,XBORDB,INAB,AUX4,0,ARGINT)

         DO 111 I=1,INAB
         XXX=AUX4(I)
!                                       ***       PRESSURE
         CALL SELINT(AUX3,AUX31,1,1,I70,YYY,XXX,2)
         AUX41(1,I)=YYY(1)
  111    CONTINUE

         DO 222 I=KI29P1,IVAR
         ENGD(I) = AUX4(I-KI29)
         PRESS(I)= AUX41(1,I-KI29)
  222    CONTINUE
         IEP=IVAR
         ICOMB=-1

         IF (IWT /= 0) THEN
         DO 444 I=1,IVAR
 444     PRINT*,' I=',I,' E=',ENGD(I),'  P=',PRESS(I)
         END IF


      END IF



      RETURN

      END


!BPS71                 *** *****
      SUBROUTINE BPS71(IEP,ICOMB)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATE THE COMBINED EOS AS USED BY  BAYM, PETHICK, AND
!          SUTHERLAND  IN [BPS71], i.e.,
!          - FEYNMAN-METROPOLIS-TELLER (OUTER CRUST)
!          - BAYM-PETHICK-SUTHERLAND (OUTER CRUST)
!          - BAYM-BETHE-PETHICK (INNER REGION)
!          - PANDHARIPANDE 'C' (INNER REGION)
!
!
! RETURN: EOS, AVAILABLE VIA  COMMON/EP1/
!
!-----------------------------------------------------------------------


      PARAMETER (NP100=100,NPVAR=3*NP100)
      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)
      COMMON/AX2/AUX100(NPVAR)

      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST

      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10


!>> 1) COMPUTE FEYNMAN-METROPOLIS-TELLER EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL FMT49(AUX3,AUX4,AUX1,IDIM,IFMT)
      DO 100 I=1,IFMT
      ENGD(I) =AUX3(I)
      PDEN(I) =AUX1(I)
  100 PRESS(I)=AUX4(I)


!>> 2) COMPUTE BAYM-PETHICK-SUTHERLAND EOS
      IDIM=IVAR
!                **** **** ****      ******
      CALL BPSCR(AUX3,AUX4,AUX1,IDIM,IBPSCR)
      ITOTAL=IFMT+IBPSCR
      IF (ITOTAL > IVAR) GOTO 9999
      DO 200 I=1,IBPSCR
      ENGD(I+IFMT) =AUX3(I)
      PDEN(I+IFMT) =AUX1(I)
  200 PRESS(I+IFMT)=AUX4(I)


!>> 3) COMPUTE BAYM-BETHE-PETHICK EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL BBP71(AUX3,AUX4,AUX1,IDIM,IBBP)
      ITOTAL=IFMT+IBPSCR+IBBP
      IF (ITOTAL > IVAR) GOTO 9999
      DO 300 I=1,IBBP
      ENGD(I+IFMT+IBPSCR) =AUX3(I)
      PDEN(I+IFMT+IBPSCR) =AUX1(I)
  300 PRESS(I+IFMT+IBPSCR)=AUX4(I)


!>> 4) COMPUTE PANDHARIPANDE EOS 'C'
      IDIM=IVAR
!                **** **** ****      ****
      CALL PAN71(AUX3,AUX4,AUX1,IDIM,IPAN)
      ITOTAL=IFMT+IBPSCR+IBBP+IPAN
      IF (ITOTAL > IVAR) GOTO 9999
      DO 400 I=1,IPAN
      ENGD(I+IFMT+IBPSCR+IBBP) =AUX3(I)
      PDEN(I+IFMT+IBPSCR+IBBP) =AUX1(I)
  400 PRESS(I+IFMT+IBPSCR+IBBP)=AUX4(I)

      IEP=ITOTAL
      ICOMB=-1


!>> CHECK ORDERING OF ENERGY DENSITY AND PRESSURE
!                              *****
      CALL SORTOR(ENGD,IEP,-10,IERSO)
      IERSO1=IERSO
!                               *****
      CALL SORTOR(PRESS,IEP,-10,IERSO)
      IERSO2=IERSO

      IF (IERSO2 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _BBP71_ PRESSURE NOT MONOTONOUSLY'',
     +       '' INCREASING'')')
      WRITE (6,505)
      END IF

      IF (IERSO1 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _BPS_ ENERGY DENSITY NOT MONOTONOUSLY'',
     +       '' INCREASING'',//)')
      DO 510 I=1,IEP
  510 WRITE (6,515) I,ENGD(I),PRESS(I)
  515 FORMAT (' I=',I4,2X,'E=',E13.6,2X,'P=',E13.6)
      STOP
      END IF
  505 FORMAT (//,'**************************************************',
     +        '************************',//)
      RETURN


 9999 CONTINUE
!__ ERROR STOP
      WRITE (6,'(//,'' DIMENSION OF AUXILIARY ARRAYS AUX3 AND AUX4 IS'',
     +          '' TOO SMALL!'',/,'' ITOTAL='',I4,'' IFMT='',I4,'' IB'',
     +          ''PSCR='',I4,'' IBBP='',I4,'' IPAN='',I4,'' IVAR='',I4,/
     +          )') ITOTAL,IFMT,IBPSCR,IBBP,IPAN,IVAR
      STOP

      END

!MISHQ                 *** *****
      SUBROUTINE MISHQ(IEP,ICOMB)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          compute BPS eos at sub-nuclear densities; combine it with
!          Mishra's quark matter eos.
!          Note: The phase transition is determined inconsistently (only
!          one conserved charge rather than two!). For that reason a
!          energy gap occurs: 226 MeV/fm^3  -->  284 MeV/fm^3, which 
!          corresponds to a particle density gap of
!                             0.24 1/fm^3   -->  0.33 MeV/fm^3.
!          The pressure's value is  4 MeV/fm^3 (=constant)
!
! RETURN: EOS, AVAILABLE VIA  COMMON/EP1/
!
!-----------------------------------------------------------------------


      PARAMETER (NP100=100,NPVAR=3*NP100)
      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)

      PARAMETER (m651=160501)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA3/AUX61(1,M651),AUX6(M651),DPDRAX(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)

      COMMON/AX2/AUX100(NPVAR)
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT

      CHARACTER * 30 ARGINT

      parameter (ipemi=62)
      dimension emi(ipemi),pmi(ipemi),dmi(ipemi)

!  energy density  (MeV/fm^3)
      data (emi(i),i=1,ipemi)/
     +.26550E+03,.28860E+03,.31370E+03,.34088E+03,
     +.37025E+03,.40196E+03,.43604E+03,.47267E+03,.51197E+03,.55400E+03,
     +.59885E+03,.64685E+03,.69785E+03,.75215E+03,.81005E+03,.87125E+03,
     +.93605E+03,.10047E+04,.10776E+04,.11541E+04,.12351E+04,.13206E+04,
     +.14103E+04,.15048E+04,.16044E+04,.17088E+04,.18183E+04,.19332E+04,
     +.20535E+04,.21795E+04,.23115E+04,.24492E+04,.25935E+04,.27438E+04,
     +.29007E+04,.30645E+04,.32349E+04,.34128E+04,.35988E+04,.37908E+04,
     +.39888E+04,.41988E+04,.44148E+04,.46398E+04,.48708E+04,.51138E+04,
     +.53628E+04,.56238E+04,.58908E+04,.61698E+04,.64608E+04,.67578E+04,
     +.70668E+04,.73848E+04,.77148E+04,.80568E+04,.84078E+04,.87708E+04,
     +.91458E+04,.95328E+04,.99318E+04,.10343E+05/

!  pressure  (MeV/fm^3)
      data (pmi(i),i=1,ipemi)/
     +.88340E+00,.85850E+01,.16950E+02,.26010E+02,
     +.35800E+02,.46370E+02,.57730E+02,.69940E+02,.83040E+02,.97050E+02,
     +.11200E+03,.12800E+03,.14500E+03,.16310E+03,.18240E+03,.20280E+03,
     +.22440E+03,.24730E+03,.27160E+03,.29710E+03,.32410E+03,.35260E+03,
     +.38250E+03,.41400E+03,.44720E+03,.48200E+03,.51850E+03,.55680E+03,
     +.59690E+03,.63890E+03,.68290E+03,.72880E+03,.77690E+03,.82700E+03,
     +.87930E+03,.93390E+03,.99070E+03,.10500E+04,.11120E+04,.11760E+04,
     +.12420E+04,.13120E+04,.13840E+04,.14590E+04,.15360E+04,.16170E+04,
     +.17000E+04,.17870E+04,.18760E+04,.19690E+04,.20660E+04,.21650E+04,
     +.22680E+04,.23740E+04,.24840E+04,.25980E+04,.27150E+04,.28360E+04,
     +.29610E+04,.30900E+04,.32230E+04,.33600E+04/

!  baryon density  (1/fm^3)
      data (dmi(i),i=1,ipemi)/
     +.32020E+00,.34860E+00,.37850E+00,.40990E+00,
     +.44300E+00,.47760E+00,.51390E+00,.55190E+00,.59170E+00,.63320E+00,
     +.67650E+00,.72170E+00,.76870E+00,.81770E+00,.86860E+00,.92160E+00,
     +.97650E+00,.10340E+01,.10930E+01,.11540E+01,.12170E+01,.12830E+01,
     +.13510E+01,.14210E+01,.14940E+01,.15690E+01,.16460E+01,.17260E+01,
     +.18080E+01,.18930E+01,.19800E+01,.20700E+01,.21630E+01,.22580E+01,
     +.23560E+01,.24570E+01,.25610E+01,.26670E+01,.27760E+01,.28880E+01,
     +.30040E+01,.31220E+01,.32430E+01,.33670E+01,.34940E+01,.36240E+01,
     +.37580E+01,.38940E+01,.40340E+01,.41780E+01,.43240E+01,.44740E+01,
     +.46270E+01,.47840E+01,.49440E+01,.51080E+01,.52750E+01,.54460E+01,
     +.56210E+01,.57990E+01,.59800E+01,.61660E+01/


      iemi=ipemi-1

! copy Mishra's eos
      do 601 i=1,iemi
      emi(i)  = emi(i+1)
      pmi(i)  = pmi(i+1)
      dmi(i)  = dmi(i+1)
  601 continue


!>> 1) COMPUTE FEYNMAN-METROPOLIS-TELLER EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL FMT49(AUX3,AUX4,AUX1,IDIM,IFMT)
      DO 100 I=1,IFMT
      ENGD(I) =AUX3(I)
      PDEN(I) =AUX1(I)
  100 PRESS(I)=AUX4(I)


!>> 2) COMPUTE BAYM-PETHICK-SUTHERLAND EOS
      IDIM=IVAR
!                **** **** ****      ******
      CALL BPSCR(AUX3,AUX4,AUX1,IDIM,IBPSCR)
      ITOTAL=IFMT+IBPSCR
      IF (ITOTAL > IVAR) GOTO 9999
      DO 200 I=1,IBPSCR
      ENGD(I+IFMT) =AUX3(I)
      PDEN(I+IFMT) =AUX1(I)
  200 PRESS(I+IFMT)=AUX4(I)


!>> 3) COMPUTE BAYM-BETHE-PETHICK EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL BBP71(AUX3,AUX4,AUX1,IDIM,IBBP)
      ITOTAL=IFMT+IBPSCR+IBBP
      IF (ITOTAL > IVAR) GOTO 9999
      DO 300 I=1,IBBP
      ENGD(I+IFMT+IBPSCR) =AUX3(I)
      PDEN(I+IFMT+IBPSCR) =AUX1(I)
  300 PRESS(I+IFMT+IBPSCR)=AUX4(I)


      pc1=3.
      do 710 i=1,itotal
      ii=i
      if(press(i) >= pc1) goto 720
  710 continue
      write (6,'(/,'' fatal error: P(BPS) < P_const'',/,
     +             '' error stop in subroutine MISHQ'')')
      stop 'error stop in MISHQ'
  720 im1    = ii - 1


!  compute energy density grid
      i210= 50
      a   = engd(im1)
      b   = emi(1)
      argint='MISHQ: eng. grid points Mishra'
!                        ****
      call intv(a,b,i210,aux7,0,argint)
      do 730 i=1,i210
  730 engd(im1-1+i) = aux7(i)
      itotal = im1 - 1 + i210

!  pressure
      aux6(1)    = press(im1)
!  p=const(!); allow for a small increase of pressure, delp (MeV/fm^3)
      delp       = 0.08
      pmi(1)     = press(im1) + delp
      aux6(i210) = pmi(1)

! linear interpolation of P(eps)
      xxm      = (aux6(i210)-aux6(1)) / (aux7(i210)-aux7(1))
      do 762  i=1,i210
  762 press(im1-1+i) = xxm * (aux7(i)-aux7(1)) + press(im1)


!  baryon density
      a = pden(im1)
      b = dmi(1)*(1.-0.0001)
      argint='MISHQ: baryon density   Mishra'
!                        ****
      call intv(a,b,i210,aux5,0,argint)
      do 763 i=1,i210
  763 pden(im1-1+i) = aux5(i)


!  internal energy
      a =engd(im1)        - amneut*pden(i210)
      b =engd(im1-1+i210) - amneut*pden(im1-1+i210)
      b = b*(1.-0.0001)
!     argint='MISHQ: internal energy  Mishra'
!                        ****
!     call intv(a,b,i210,aux4,0,argint)
      do 764 i=1,i210
! 764 eintnl(im1-1+i) = aux4(i)
  764 eintnl(im1-1+i) = 0.


!>> 4) ADD MISHRA'S QUARK MATTER EOS
      IF (ITOTAL > IVAR) GOTO 9999
      DO 400 I=1,IEMI
      ENGD(I+ITOTAL-1) =EMI(I)
      PDEN(I+ITOTAL-1) =DMI(I)
  400 PRESS(I+ITOTAL-1)=PMI(I)

      ITOTAL = ITOTAL + IEMI - 1
      IEP    = ITOTAL
      ICOMB  = -1



!>> CHECK ORDERING OF ENERGY DENSITY AND PRESSURE
!                              *****
      CALL SORTOR(ENGD,IEP,-10,IERSO)
      IERSO1=IERSO
!                               *****
      CALL SORTOR(PRESS,IEP,-10,IERSO)
      IERSO2=IERSO

      IF (IERSO2 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _BBP71_ PRESSURE NOT MONOTONOUSLY'',
     +       '' INCREASING'')')
      WRITE (6,505)
      END IF

      IF (IERSO1 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _BPS_ ENERGY DENSITY NOT MONOTONOUSLY'',
     +       '' INCREASING'',//)')
      DO 510 I=1,IEP
  510 WRITE (6,515) I,ENGD(I),PRESS(I)
  515 FORMAT (' I=',I4,2X,'E=',E13.6,2X,'P=',E13.6)
      STOP
      END IF
  505 FORMAT (//,'**************************************************',
     +        '************************',//)
      RETURN


 9999 CONTINUE
!__ ERROR STOP
      WRITE (6,'(//,'' DIMENSION OF AUXILIARY ARRAYS AUX3 AND AUX4 IS'',
     +          '' TOO SMALL!'',/,'' ITOTAL='',I4,'' IFMT='',I4,'' IB'',
     +          ''PSCR='',I4,'' IBBP='',I4,'' IVAR='',I4,/)') ITOTAL,IFM
     +          T,IBPSCR,IBBP,IVAR
      STOP

      END




!Can75                 *** *****
      SUBROUTINE Can75(IEP,ICOMB)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          Compute the most favored eos of Canuto, 1975, consisting of
!          the following ones:
!          - FEYNMAN-METROPOLIS-TELLER (OUTER CRUST)
!          - BAYM-PETHICK-SUTHERLAND (OUTER CRUST)
!          - BAYM-BETHE-PETHICK (INNER REGION)
!          - BETHE AND JOHNSON, pure neutron gas (INNER REGION)
!
!
! RETURN: EOS, AVAILABLE VIA  COMMON/EP1/
!
!-----------------------------------------------------------------------


      PARAMETER (NP100=100,NPVAR=3*NP100)
      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)
      COMMON/AX2/AUX100(NPVAR)

      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST

      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10


!>> 1) COMPUTE FEYNMAN-METROPOLIS-TELLER EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL FMT49(AUX3,AUX4,AUX1,IDIM,IFMT)
      DO 100 I=1,IFMT
      ENGD(I) =AUX3(I)
      PDEN(I) =AUX1(I)
  100 PRESS(I)=AUX4(I)


!>> 2) COMPUTE BAYM-PETHICK-SUTHERLAND EOS
      IDIM=IVAR
!                **** **** ****      ******
      CALL BPSCR(AUX3,AUX4,AUX1,IDIM,IBPSCR)
      ITOTAL=IFMT+IBPSCR
      IF (ITOTAL > IVAR) GOTO 9999
      DO 200 I=1,IBPSCR
      ENGD(I+IFMT) =AUX3(I)
      PDEN(I+IFMT) =AUX1(I)
  200 PRESS(I+IFMT)=AUX4(I)


!>> 3) COMPUTE BAYM-BETHE-PETHICK EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL BBP71(AUX3,AUX4,AUX1,IDIM,IBBP)
      IBBP8 =IBBP-8
      ITOTAL=IFMT+IBPSCR+IBBP8
      IF (ITOTAL > IVAR) GOTO 9999
      DO 300 I=1,IBBP8
      ENGD(I+IFMT+IBPSCR) =AUX3(I)
      PDEN(I+IFMT+IBPSCR) =AUX1(I)
  300 PRESS(I+IFMT+IBPSCR)=AUX4(I)


!>> 4) COMPUTE BETHE-JOHNSON (PURE NEUTRON GAS)
      IDIM=IVAR
!               **** **** ****      *****
      CALL BJNG(AUX3,AUX4,AUX1,IDIM,IBJNG)
      ITOTAL=IFMT+IBPSCR+IBBP8+IBJNG
      IF (ITOTAL > IVAR) GOTO 9999
      DO 400 I=1,IBJNG
      ENGD(I+IFMT+IBPSCR+IBBP8) =AUX3(I)
      PDEN(I+IFMT+IBPSCR+IBBP8) =AUX1(I)
  400 PRESS(I+IFMT+IBPSCR+IBBP8)=AUX4(I)

      IEP=ITOTAL
      ICOMB=-1


!>> CHECK ORDERING OF ENERGY DENSITY AND PRESSURE
!                              *****
      CALL SORTOR(ENGD,IEP,-10,IERSO)
      IERSO1=IERSO
!                               *****
      CALL SORTOR(PRESS,IEP,-10,IERSO)
      IERSO2=IERSO

      IF (IERSO2 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _CAN75_ PRESSURE NOT MONOTONOUSLY'',
     +       '' INCREASING'')')
      WRITE (6,505)
      END IF

      IF (IERSO1 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _CAN75_ ENERGY DENSITY NOT MONOTONOUSLY'',
     +       '' INCREASING'',//)')
      WRITE (6,505)
      DO 510 I=1,IEP
  510 WRITE (6,515) I,ENGD(I),PRESS(I),PDEN(I)
  515 FORMAT (' I=',I4,2X,'e=',E13.6,2X,'P=',E13.6,2X,'n=',E13.6)
      STOP
      END IF
  505 FORMAT (//,'**************************************************',
     +        '************************',//)
      RETURN


 9999 CONTINUE
!__ ERROR STOP
      WRITE (6,'(//,'' DIMENSION OF AUXILIARY ARRAYS AUX3 AND AUX4 IS'',
     +          '' TOO SMALL!'',/,'' ITOTAL='',I4,'' IFMT='',I4,'' IB'',
     +          ''PSCR='',I4,'' IBBP='',I4,'' IPAN='',I4,'' IVAR='',I4,/
     +          )') ITOTAL,IFMT,IBPSCR,IBBP,IBJNG,IVAR
      STOP

      END


!BJNG                 * * ** ****
      SUBROUTINE BJNG(E,P,RN,IDIM,IOTCH)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          COMPUTE THE EOS OF  BETHE AND JOHNSON, PURE NEUTRON GAS, 1974
!
!
! RETURN: (1) NUMBER OF PRESSURE AND ENERGY DENSITY VALUES (=IOTCH)
!
!         (2) EQUATION OF STATE, i.e.,
!             - E(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - P(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - N(I)  CONTAINS THE BARYON DENSITY IN  1/fm**3
!
!-----------------------------------------------------------------------

      REAL * 8 AP


      PARAMETER (IP20=20)
      DIMENSION P(IDIM),E(IDIM),RN(IDIM),AP(IP20),AE(IP20),AN(IP20)

      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10




!__ PRESSURE   (in  dynes/cm^2)
      DATA (AP(I), I=1,IP20)/1.19E33,2.93E33,6.00E33,1.09E34,
     +                       1.83E34,4.09E34,7.61E34,1.26E35,
     +                       1.99E35,2.85E35,3.71E35,4.92E35,
     +                       6.23E35,8.58E35,1.14E36,1.34E36,
     +                       1.85E36,2.76E36,4.83E36,7.62E36/

!__ ENERGY DENSITY   (in  MeV/fm^3)
      DATA (AE(I), I=1,IP20)/12.6,16.6,21.2,26.0,32.2,46.9,64.4,
     +                       83.7,109.0,135.0,160.0,190.0,224.0,
     +                       274.0,327.0,360.0,442.0,560.0,788.0,
     +                       1040.0/

!__ BARYON DENSITY   (in  1/fm^3)
      DATA (AN(I), I=1,IP20)/0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,
     +                       0.8,0.9,1.0,1.1,1.25,1.4,1.5,1.7,2.0,
     +                       2.5,3.0/



      II8=IP20
      IF (II8 > IDIM) THEN
      WRITE (6,'(//,'' FATAL DIMENSIONING ERROR IN _BJNG_'',/,
     +           '' II8='',I4,'' <!  IDIM='',I4)') II8,IDIM
      STOP
      END IF

      IWT=0

!  COPY INPUT DATA AND CONVERT UNITS
      DO 1000 I=1,II8
      RN(I) = AN(I)
!     CONVERT  dynes/cm**2   -->    MeV/fm**3
      P(I)  = AP(I) / UF4
!     CONVERT  E/B  IN MeV   -->    e  IN MeV/fm**3
      E(I)  = ( AE(I) + AMNEUT ) * RN(I)

      IF (IWT /= 0) THEN
      IF (I == 1) WRITE (6,444)
      WRITE (6,333) I,P(I),E(I),RN(I)
      END IF
 1000 CONTINUE

  333 FORMAT (' I=',I4,' P=',E12.6,' E=',E12.6,' n=',E12.6)
  444 FORMAT (/,' Bethe-Johnson eos (pure neutron gas)',//)

      IOTCH=II8


      RETURN

      END


!FP81                 *** *****
      SUBROUTINE FP81(IEP,ICOMB,SC,IVAR,I100,ICCORR)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATE THE COMBINED EOS ACCORDING TO  [FP81].
!          NOTE: CHOOSE BETWEEN THE USE OF  BBP  OR  NV!
!          THE DIFFERENT EOS ARE:
!          - FEYNMAN-METROPOLIS-TELLER (OUTER CRUST)
!          - BAYM-PETHICK-SUTHERLAND (OUTER CRUST)
!
!                  - BAYM-BETHE-PETHICK \
!            OR                           (INNER REGION)
!                  - NEGELE-VAUTHERIN   /
!
!          - FRIEDMAN-PANDHARIPANDE (CORE REGION): V14+TNI MODEL
!            (CONTAINS TWO- AND THREE-NUCLEON INTERACTIONS)
!
!
! RETURN: EOS, AVAILABLE VIA  COMMON/EP1/
!
!-----------------------------------------------------------------------

      REAL SC(4,IVAR)

      PARAMETER (NP100=100,NPVAR=3*NP100)
      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)

      PARAMETER (m651=160501)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)


      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT

      CHARACTER * 30 ARGINT

      DATA ESFP81/789.736/,EPFP81/182.215/


!>> 1) COMPUTE FEYNMAN-METROPOLIS-TELLER EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL FMT49(AUX3,AUX4,AUX1,IDIM,IFMT)
      DO 100 I=1,IFMT
      ENGD(I) =AUX3(I)
      PDEN(I) =AUX1(I)
      PRESS(I)=AUX4(I)
  100 CONTINUE

!>> 2) COMPUTE BAYM-PETHICK-SUTHERLAND EOS
      IDIM=IVAR
!                **** **** ****      ******
      CALL BPSCR(AUX3,AUX4,AUX1,IDIM,IBPSCR)
      ITOTAL=IFMT+IBPSCR
      IF (ITOTAL > IVAR) GOTO 9999
      DO 200 I=1,IBPSCR
      ENGD(I+IFMT) =AUX3(I)
      PDEN(I+IFMT) =AUX1(I)
      PRESS(I+IFMT)=AUX4(I)
  200 CONTINUE


      IBBPE=-10
      IF (IBBPE > 0) THEN
!>> 3A) COMPUTE BAYM-BETHE-PETHICK EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL BBP71(AUX3,AUX4,AUX1,IDIM,IBBP)
      ITOTAL=IFMT+IBPSCR+IBBP
      IF (ITOTAL > IVAR) GOTO 9999
      DO 300 I=1,IBBP
      ENGD(I+IFMT+IBPSCR) =AUX3(I)
      PDEN(I+IFMT+IBPSCR) =AUX1(I)
  300 PRESS(I+IFMT+IBPSCR)=AUX4(I)

      ELSE

!>> 3B) COMPUTE NEGELE-VAUTHERIN EOS
      IDIM=IVAR
      A=4.0E-04
      B=0.1E-02
!      B=5.0E-02
!      B=1.0E-02
      NAB=I100
      ARGINT='FP81: ENERGY DENSITY GRID     '
      CALL INTV(A,B,NAB,AUX3,0,ARGINT)
!__ RETURN: AUX3  (BARYON DENSITY IN  1/fm**3)
      DO 500 I=1,NAB
      RFM=AUX3(I)
      CALL NVEOS(AMNEUT,RFM,EMEVFM,PMEVFM)
!__ RETURN: ENERGY DENSITY  EMEVFM  (MeV/fm**3)
!           PRESSURE        PMEVFM  (MeV/fm**3)

      AUX31(1,I)=EMEVFM
      AUX41(1,I)=PMEVFM
  500 CONTINUE
      ITOTAL=IFMT+IBPSCR+NAB
      IF (ITOTAL > IVAR) GOTO 9999
      DO 301 I=1,NAB
      ENGD(I+IFMT+IBPSCR) =AUX31(1,I)
      PDEN(I+IFMT+IBPSCR) =0.
      PRESS(I+IFMT+IBPSCR)=AUX41(1,I)
  301 CONTINUE
      END IF


!>> 4) COMPUTE FRIEDMAN-PANDHARIPANDE EOS (V14+TNI MODEL)
      IDIM=IVAR
!                 **** ****      ****
      CALL FRPA81(AUX3,AUX4,IDIM,IPAN,SC,IVAR)


      IF (ICCORR /= 0) THEN
! CORRECT EOS FOR CAUSALITY

      ESW = ESFP81
      EPW = EPFP81
      ENW = 1.0

      DO 1 I=1,IPAN
      I1P   =I-1
    1 IF (AUX3(I) > ESW) GOTO 2
    2 CONTINUE
      IF (I1P == IPAN) THEN
      WRITE (6,'(///,'' _FP81_  I1P=IPAN --> CAUSALITY CORRECTION'',
     +               '' FAILS'',//,'' ERROR STOP'')')
      STOP
      END IF

! COMPUTE DENSITY GRID
      ID1P = I100 - I1P
      XNA  = ENW  * (1. + .001)
      XNB  = 2.   * ENW
      ARGINT='FP81: CAUSALITY CORRECTION    '
!                            ****      BARYON DENSITY GRID  n(i)
      CALL INTV(XNA,XNB,ID1P,AUX1,0,ARGINT)

! COMPUTE CAUSALITY CORRECTED PART OF THE EOS
      DO 3 I=1,ID1P
      X       = AUX1(I)/ENW
      PR      = (EPW-ESW + (EPW+ESW) * X * X ) / 2.
      AUX8(I) = PR
      EY      = ESW - EPW + PR
      AUX7(I) = EY
    3 CONTINUE

! COPY
      DO 4       I=1,ID1P
      AUX3(I+I1P) =  AUX7(I)
      AUX4(I+I1P) =  AUX8(I)
    4 CONTINUE

      IPAN = I100

! END OF CAUSALITY CORRECTION
                                 END IF



      IF (IBBPE > 0) THEN
      IV00=IBBP
      ELSE
      IV00=NAB
      END IF
      ITOTAL=IFMT+IBPSCR+IV00+IPAN
      IF (ITOTAL > IVAR) GOTO 9999
      DO 400 I=1,IPAN
      ENGD(I+IFMT+IBPSCR+IV00) =AUX3(I)
      PDEN(I+IFMT+IBPSCR+IV00) =0.
      PRESS(I+IFMT+IBPSCR+IV00)=AUX4(I)
  400 CONTINUE

      IEP=ITOTAL
      ICOMB=-1


!>> CHECK ORDERING OF ENERGY DENSITY AND PRESSURE
!                              *****
      CALL SORTOR(ENGD,IEP,-10,IERSO)
      IERSO1=IERSO
!                               *****
      CALL SORTOR(PRESS,IEP,-10,IERSO)
      IERSO2=IERSO

      IF (IERSO2 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _BBP71_ PRESSURE NOT MONOTONOUSLY'',
     +       '' INCREASING'')')
      WRITE (6,505)
      END IF

      IF (IERSO1 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _BPS_ ENERGY DENSITY NOT MONOTONOUSLY'',
     +       '' INCREASING'',//)')
      DO 510 I=1,IEP
  510 WRITE (6,515) I,ENGD(I),PRESS(I)
  515 FORMAT (' I=',I4,2X,'E=',E13.6,2X,'P=',E13.6)
      STOP
      END IF
  505 FORMAT (//,'**************************************************',
     +        '************************',//)
      RETURN



 9999 CONTINUE
!__ ERROR STOP
      WRITE (6,'(//,'' DIMENSION OF AUXILIARY ARRAYS AUX3 AND AUX4 IS'',
     +          '' TOO SMALL!'',/,'' ITOTAL='',I4,'' IFMT='',I4,'' IB'',
     +          ''PSCR='',I4,'' IBBP='',I4,'' IPAN='',I4,'' IVAR='',I4,/
     +          )') ITOTAL,IFMT,IBPSCR,IBBP,IPAN,IVAR
      STOP

      END


!FMT49                 * * **      *****
      SUBROUTINE FMT49(E,P,RN,IDIM,IOTCH)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATE THE EOS OF  FEYNMAN-METROPOLIS-TELLER  (1949)
!
!
! RETURN:  NUMBER OF PRESSURE RESPECTIVELY ENERGY DENSITY VALUES (IOTCH)
!
!          AND EQUATION OF STATE, i.e.,
!             - E(I),  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - P(I),  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!
!-----------------------------------------------------------------------

      REAL * 8 AN  
      DOUBLE PRECISION BD


      PARAMETER (IP8=8)
      DIMENSION AP(IP8),AE(IP8),AN(IP8)

      DIMENSION P(IDIM),E(IDIM),RN(IDIM)


      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10

      DATA BD/1.0D-39/  



!__ PRESSURE   (in  dynes/cm^2)
      DATA (AP(I), I=1,IP8)/
     +1.01E09,1.01E10,1.01E11,1.21E12,1.40E13,1.70E14,5.82E15,
     +1.90E17/


!__ MASS DENSITY   (in  g/cm^3)
      DATA (AE(I), I=1,IP8)/
     +7.86,7.90,8.15,11.6,16.4,45.1,212.,1150./


!__ BARYON DENSITY   (in  1/cm^3)
      DATA (AN(I), I=1,IP8)/4.73E24,4.76E24,4.91E24,6.99E24,
     +                      9.90E24,2.72E25,1.27E26,6.93E26/



      II8=IP8
      IF (II8 > IDIM) THEN
      WRITE (6,'(//,'' FATAL DIMENSIONING ERROR IN _FMT49_'',/,
     +           '' II8='',I4,'' <!  IDIM='',I4)') II8,IDIM
      STOP
      END IF


      IWT=0
      DO 1000 I=1,II8

!  CONVERT  dynes/cm**2  -->  MeV/fm**3
!  CONVERT  g/cm**3      -->  MeV/fm**3
!  CONVERT  1/cm**3      -->    1/fm**3
      P(I) = AP(I) / UF4
      E(I) = AE(I) / UF6
      RN(I)= SNGL ( AN(I) ) * BD

      IF (IWT /= 0) THEN
      IF (I == 1) WRITE (6,444)
      WRITE (6,333) I,P(I),E(I),RN(I)
      END IF

 1000 CONTINUE

  333 FORMAT (' I=',I4,' P=',E12.6,' E=',E12.6,' N=',E12.6)
  444 FORMAT (/,' FEYNMAN-METROPOLIS-TELLER  EQUATION OF STATE',//)

      IOTCH=II8

      RETURN

      END

!BPSCR                 * * **      *****
      SUBROUTINE BPSCR(E,P,RN,IDIM,IOTCH)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATE THE EOS OF  BAYM-PETHICK-SUTHERLAND  (1971)
!
!
! RETURN: NUMBER OF PRESSURE RESPECTIVELY ENERGY DENSITY VALUES (IOTCH)
!
!         AND EQUATION OF STATE, i.e.,
!             - E(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - P(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - N(I)  CONTAINS THE BARYON DENSITY IN    1/fm**3
!
!-----------------------------------------------------------------------

      DOUBLE PRECISION BD

      PARAMETER (IP43=43)
      REAL * 8 AP(IP43),AE(IP43),AN(IP43)

      DIMENSION P(IDIM),E(IDIM),RN(IDIM)


      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10

      DATA BD/1.0D-39/


!__ PRESSURE   (in  dynes/cm^2)
      DATA (AP(I), I=1,IP43)/
     +9.744E18,4.968E19,2.431E20,1.151E21,5.266E21,2.318E22,9.755E22,
     +3.911E23,5.259E23,1.435E24,3.833E24,1.006E25,2.604E25,6.676E25,
     +8.738E25,1.629E26,3.029E26,4.129E26,5.036E26,6.860E26,1.272E27,
     +2.356E27,4.362E27,5.662E27,7.702E27,1.048E28,1.425E28,1.938E28,
     +2.503E28,3.404E28,4.628E28,5.949E28,8.089E28,1.100E29,1.495E29,
     +2.033E29,2.597E29,2.892E29,3.290E29,4.473E29,5.816E29,7.538E29,
     +7.805E29/


!__ MASS DENSITY   (in  g/cm^3)
      DATA (AE(I), I=1,IP43)/
     +1.044E04,2.622E04,6.587E04,1.654E05,4.156E05,1.044E06,2.622E06,
     +6.588E06,8.293E06,1.655E07,3.302E07,6.589E07,1.315E08,2.624E08,
     +3.304E08,5.237E08,8.301E08,1.045E09,1.316E09,1.657E09,2.626E09,
     +4.164E09,6.601E09,8.312E09,1.046E10,1.318E10,1.659E10,2.090E10,
     +2.631E10,3.313E10,4.172E10,5.254E10,6.617E10,8.332E10,1.049E11,
     +1.322E11,1.664E11,1.844E11,2.096E11,2.640E11,3.325E11,4.188E11,
     +4.299E11/


!__ PARTICLE DENSITY   (in  1/cm^3)
      DATA (AN(I), I=1,IP43)/
     +6.295E27,1.581E28,3.972E28,9.976E28,2.506E29,6.294E29,1.581E30,
     +3.972E30,5.000E30,9.976E30,1.990E31,3.972E31,7.924E31,1.581E32,
     +1.990E32,3.155E32,5.000E32,6.294E32,7.924E32,9.976E32,1.581E33,
     +2.506E33,3.972E33,5.000E33,6.294E33,7.924E33,9.976E33,1.256E34,
     +1.581E34,1.990E34,2.506E34,3.155E34,3.972E34,5.000E34,6.294E34,
     +7.924E34,9.976E34,1.105E35,1.256E35,1.581E35,1.990E35,2.506E35,
     +2.572E35/



      II8=IP43
      IF (II8 > IDIM) THEN
      WRITE (6,'(//,'' FATAL DIMENSIONING ERROR IN _BPSCR_'',/,
     +           '' II8='',I4,'' <!  IDIM='',I4)') II8,IDIM
      STOP
      END IF


      IWT=0
      DO 1000 I=1,II8
!  CONVERT  dynes/cm**2  -->  MeV/fm**3
!  CONVERT  g/cm**3      -->  MeV/fm**3
!  CONVERT  1/cm**3      -->    1/fm**3
       P(I) = SNGL ( AP(I) ) / UF4
       RN(I)= SNGL ( AN(I) * BD )

      E(I) = AE(I) / UF6
      IF (IWT /= 0) THEN
      IF (I == 1) WRITE (6,444)
      WRITE (6,333) I,P(I),E(I),RN(I)
      END IF

 1000 CONTINUE

  333 FORMAT (' I=',I4,' P=',E12.6,' E=',E12.6,' P=',E12.6)
  444 FORMAT (/,' BAYM-PETHICK-SUTHERLAND  EQUATION OF STATE',/,
     +          ' (INNER CRUST REGION)',//)

      IOTCH=II8


      RETURN

      END

!BBP71                 * * **      *****
      SUBROUTINE BBP71(E,P,RN,IDIM,IOTCH)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATE THE EOS OF  BAYM-BETHE-PETHICK (1971)
!
!
! RETURN: NUMBER OF PRESSURE RESPECTIVELY ENERGY DENSITY VALUES (IOTCH)
!
!         AND EQUATION OF STATE, i.e., 
!             - E(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - P(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!
!-----------------------------------------------------------------------

      REAL * 8 AP,AN
      DOUBLE PRECISION BD 


      PARAMETER (IP35=35)
      DIMENSION P(IDIM),E(IDIM),RN(IDIM),AP(IP35),AE(IP35),AN(IP35)


      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10

      DATA BD/1.0D-39/


!__ PRESSURE   (in  dynes/cm^2)
      DATA (AP(I), I=1,IP35)/
     +7.890E29,8.352E29,9.098E29,9.831E29,1.083E30,1.218E30,1.399E30,
     +1.638E30,1.950E30,2.592E30,3.506E30,4.771E30,6.481E30,8.748E30,
     +1.170E31,1.695E31,2.209E31,2.848E31,3.931E31,6.178E31,8.774E31,
     +1.386E32,1.882E32,2.662E32,3.897E32,5.861E32,8.595E32,1.286E33,
     +1.900E33,2.242E33,2.751E33,3.369E33,4.286E33,6.103E33,7.391E33/


!__ MASS DENSITY   (in  g/cm^3)
      DATA (AE(I), I=1,IP35)/
     +4.460E11,5.228E11,6.610E11,7.964E11,9.728E11,1.196E12,1.471E12,
     +1.805E12,2.202E12,2.930E12,3.833E12,4.933E12,6.248E12,7.801E12,
     +9.611E12,1.246E13,1.496E13,1.778E13,2.210E13,2.988E13,3.767E13,
     +5.081E13,6.193E13,7.732E13,9.826E13,1.262E14,1.586E14,2.004E14,
     +2.520E14,2.761E14,3.085E14,3.433E14,3.885E14,4.636E14,5.094E14/


!__ BARYON DENSITY   (in  1/cm^3)
      DATA (AN(I), I=1,IP35)/
     +2.670E35,3.126E35,3.951E35,4.759E35,5.812E35,7.143E35,8.786E35,
     +1.077E36,1.314E36,1.748E36,2.287E36,2.942E36,3.726E36,4.650E36,
     +5.728E36,7.424E36,8.907E36,1.059E37,1.315E37,1.777E37,2.239E37,
     +3.017E37,3.675E37,4.585E37,5.821E37,7.468E37,9.371E37,1.182E38,
     +1.484E38,1.625E38,1.814E38,2.017E38,2.280E38,2.715E38,2.979E38/



      II8=IP35
      IF (II8 > IDIM) THEN
      WRITE (6,'(//,'' FATAL DIMENSIONING ERROR IN _BBP71_'',/,
     +           '' II8='',I4,'' <!  IDIM='',I4)') II8,IDIM
      STOP
      END IF


      IWT=0
      DO 1000 I=1,II8

!  CONVERT  dynes/cm**2  -->  MeV/fm**3
!  CONVERT  g/cm**3      -->  MeV/fm**3
!  CONVERT  1/cm**3      -->    1/fm**3
      P(I) = SNGL ( AP(I) ) / UF4
      RN(I)= SNGL ( AN(I) * BD )

      E(I) = AE(I) / UF6
      IF (IWT /= 0) THEN
      IF (I == 1) WRITE (6,444)
      WRITE (6,333) I,P(I),E(I),RN(I)
      END IF

 1000 CONTINUE

  333 FORMAT (' I=',I4,' P=',E12.6,' E=',E12.6)
  444 FORMAT (/,' BAYM-BETHE-PETHICK  EQUATION OF STATE',//)

      IOTCH=II8


      RETURN

      END

!PAN71                 * * ** ****
      SUBROUTINE PAN71(E,P,RN,IDIM,IOTCH)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATE THE EOS OF  PANDHARIPANDE, MODEL 'C', FOR
!          HYPERONIC MATTER (1971)
!
!
! RETURN: (1) NUMBER OF PRESSURE AND ENERGY DENSITY VALUES (=IOTCH)
!
!         (2) EQUATION OF STATE, i.e., 
!             - E(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - P(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - N(I)  CONTAINS THE BARYON DENSITY IN  1/fm**3
!
!-----------------------------------------------------------------------

      REAL * 8 AN
      DOUBLE PRECISION BD


      PARAMETER (IP29=31)
      DIMENSION P(IDIM),E(IDIM),RN(IDIM),AP(IP29),AE(IP29),AN(IP29)


      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10


      DATA BD/1.0D-39/ 


!__ PRESSURE   (in  MeV/fm^3)
      DATA (AP(I), I=1,IP29)/9.465,12.87,
     +16.98,19.36,23.65,28.32,34.19,46.72,69.26,98.09,134.4,178.6,
     +256.4,352.7,469.5,613.4,849.0,1139.,1485.,2114.,2881.,3793.,
     +5266.,6974.,9516.,12583.,16067.,20121.,24611.,26997.,35238./

!__ ENERGY DENSITY   (in  MeV/fm^3)
      DATA (AE(I), I=1,IP29)/310.2,350.4,
     +390.9,442.0,493.5,545.4,597.8,704.1,812.3,924.1,1040.,1159.,
     +1346.,1544.,1753.,1975.,2292.,2637.,3012.,3637.,4342.,5134.,
     +6335.,7711.,9701.,11995.,14624.,17589.,20916.,24605.,28680./

!__ BARYON DENSITY   (in  1/cm^3)
!   DO NOT USE THE FOLLOWING VALUES. THESE DO NOT CORRESPOND TO THE
!   ABOVE GIVEN ENERGY AND PRESSURE VALUES
!      DATA (AN(I), I=1,IP29)/4.000E38,5.000E38,6.000E38,8.000E38,
!     +                       1.000E39,1.250E39,1.550E39,1.900E39,
!     +                       2.300E39,2.900E39,3.600E39,4.500E39,
!     +                       5.500E39,6.500E39,7.500E39/



      II8=IP29
      IF (II8 > IDIM) THEN
      WRITE (6,'(//,'' FATAL DIMENSIONING ERROR IN _PAN71_'',/,
     +           '' II8='',I4,'' <!  IDIM='',I4)') II8,IDIM
      STOP
      END IF

      IWT=0

!  COPY AND CONVERT UNITS
      DO 1000 I=1,II8
      P(I) = AP(I)
      E(I) = AE(I)
!     CONVERT  1/cm**3   -->    1/fm**3
      AN(I)=0.
      RN(I)= SNGL ( AN(I) ) * BD
      IF (IWT /= 0) THEN
      IF (I == 1) WRITE (6,444)
      WRITE (6,333) I,P(I),E(I)
      END IF

 1000 CONTINUE

  333 FORMAT (' I=',I4,' P=',E12.6,' E=',E12.6)
  444 FORMAT (/,' PANDHARIPANDE  EQUATION OF STATE  "C" ',//)

      IOTCH=II8


      RETURN

      END

!FRPA81                 * * ****
      SUBROUTINE FRPA81(E,P,IDIM,IOTCH,SC,IVAR)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATE THE EOS OF  FRIEDMAN-PANDHARIPANDE, V14+TNI MODEL
!          [FP81].
!
!
! RETURN: (1) NUMBER OF PRESSURE AND ENERGY DENSITY VALUES (=IOTCH),
!
!         (2) EQUATION OF STATE, i.e., 
!             - E(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - P(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!
!-----------------------------------------------------------------------

      REAL SC(4,IVAR)

      PARAMETER (IP24=28,NP100=100,NPVAR=3*NP100)
      DIMENSION P(IDIM),E(IDIM),AP(IP24),AE(IP24),AR(IP24)


      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)

      CHARACTER * 30 EASS


!__ DENSITY   (in  1/fm^3)
      DATA (AR(I), I=1,IP24)/4.222E-03,7.295E-03,1.158E-02,
     +1.729E-02,2.462E-02,3.377E-02,4.495E-02,5.836E-02,7.420E-02,
     +9.267E-02,1.140E-01,1.383E-01,1.659E-01,1.921E-01,2.223E-01,
     +2.574E-01,2.978E-01,3.451E-01,3.992E-01,4.622E-01,5.353E-01,
     +6.193E-01,7.382E-01,8.305E-01,1.000,1.250,1.500,2.000/


!__ ENERGY PER NUCLEON   (in  MeV)
      DATA (AE(I), I=1,IP24)/2.057,2.747,3.503,4.322,5.177,6.027,
     +6.873,7.744,8.699,9.822,11.22,13.01,15.26,17.69,20.73,24.67,
     +29.59,36.04,44.16,54.94,69.15,87.96,112.8,147.1,225.4,338.8,
     +491.7,921.5/


      II8=IP24
      IF (II8 > IDIM) THEN
      WRITE (6,'(//,'' FATAL DIMENSIONING ERROR IN _FRPA81_'',/,
     +           '' II8='',I4,'' <!  IDIM='',I4)') II8,IDIM
      STOP
      END IF

!  COMPUTE  D [E/N] / D [RHO]
      EASS='FRPA81: CALL SELSPL 123       '
      DO 123 I=1,II8
      XXX=AR(I)
      CALL SELSPL(AR,AE,II8,SC,10,YYY,XXX,100,EASS)
      AP(I)= AR(I)**2 * YYY
  123 CONTINUE

      IWT=0
      DO 1000 I=1,II8
      E(I) = (AE(I) + AMNEUT) * AR(I)
      P(I) = AP(I)
      IF (IWT /= 0) THEN
      IF (I == 1) WRITE (6,444)
      WRITE (6,333) I,P(I),E(I)
      END IF
 1000 CONTINUE


  333 FORMAT (' I=',I4,' P=',E12.6,' E=',E12.6)
  444 FORMAT (/,' FRIEDMAN-PANDHARIPANDE  EOS  V14+TNI',//)

      IOTCH=II8


      RETURN

      END

!WFF88                 *** *****
      SUBROUTINE WFF88(IEP,ICOMB,SC,IVAR,I100,KST1,ICCORR)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATE THE COMBINED EOS ACCORDING TO  [WFF88]
!          NOTE: CHOOSE BETWEEN THE USE OF  BBP  OR  NV!
!          THE DIFFERENT EOS ARE:
!          -  FEYNMAN-METROPOLIS-TELLER (OUTER CRUST)
!          -  BAYM-PETHICK-SUTHERLAND (OUTER CRUST)
!
!                  - BAYM-BETHE-PETHICK \
!            OR                           (INNER REGION)
!                  - NEGELE-VAUTHERIN   /
!
!           CORE REGION: WIRINGA, FIKS, AND FABROCINI [WFF88]
!           -  AV14+UVII   ><   KST1=-70,
!           -  UV14+UVII   ><   KST1=-60,
!           -  UV14+TNI    ><   KST1=-50
!
! PARAMETER:
!            ICCORR >< CORRECTION OF EOS FOR CAUSALITY IF ICCORR  /=  0
!
! RETURN: EOS, AVAILABLE VIA  COMMON/EP1/
!
!-----------------------------------------------------------------------

      REAL SC(4,IVAR)

      PARAMETER (NP100=100,NPVAR=3*NP100)
      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)

      PARAMETER (m651=160501)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)


      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT


      CHARACTER * 30 ARGINT


      DATA ESWUU/901./ ,EPWUU/326./,ENWUU/0.8088/
      DATA ESWAU/1009./,EPWAU/421./,ENWAU/0.9004/


!>> 1) COMPUTE FEYNMAN-METROPOLIS-TELLER EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL FMT49(AUX3,AUX4,AUX1,IDIM,IFMT)
      DO 100 I=1,IFMT
      ENGD(I) =AUX3(I)
      PDEN(I) =AUX1(I)
      PRESS(I)=AUX4(I)
  100 CONTINUE

!>> 2) COMPUTE BAYM-PETHICK-SUTHERLAND EOS
      IDIM=IVAR
!                **** **** ****      ******
      CALL BPSCR(AUX3,AUX4,AUX1,IDIM,IBPSCR)
      ITOTAL=IFMT+IBPSCR
      IF (ITOTAL > IVAR) GOTO 9999
      DO 200 I=1,IBPSCR
      ENGD(I+IFMT) =AUX3(I)
      PDEN(I+IFMT) =AUX1(I)
      PRESS(I+IFMT)=AUX4(I)
  200 CONTINUE


      IBBPE=-10
      IF (IBBPE > 0) THEN
!>> 3A) COMPUTE BAYM-BETHE-PETHICK EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL BBP71(AUX3,AUX4,AUX1,IDIM,IBBP)
      ITOTAL=IFMT+IBPSCR+IBBP
      IF (ITOTAL > IVAR) GOTO 9999
      DO 300 I=1,IBBP
      ENGD(I+IFMT+IBPSCR) =AUX3(I)
      PDEN(I+IFMT+IBPSCR) =AUX1(I)
  300 PRESS(I+IFMT+IBPSCR)=AUX4(I)

      ELSE

!>> 3B) COMPUTE NEGELE-VAUTHERIN EOS
      IDIM=IVAR
      A=4.0E-04
      B=0.1E-02
!      B=5.0E-02
!      B=1.0E-02
      NAB=I100
      ARGINT='WFF81: ENERGY DENSITY GRID    '
      CALL INTV(A,B,NAB,AUX3,0,ARGINT)
!__ RETURN: AUX3  (BARYON DENSITY IN  1/fm**3)
      DO 500 I=1,NAB
      RFM=AUX3(I)
      CALL NVEOS(AMNEUT,RFM,EMEVFM,PMEVFM)
!__ RETURN: ENERGY DENSITY  EMEVFM  (MeV/fm**3)
!           PRESSURE        PMEVFM  (MeV/fm**3)

      AUX31(1,I)=EMEVFM
      AUX41(1,I)=PMEVFM
  500 CONTINUE
      ITOTAL=IFMT+IBPSCR+NAB
      IF (ITOTAL > IVAR) GOTO 9999
      DO 301 I=1,NAB
      ENGD(I+IFMT+IBPSCR) =AUX31(1,I)
      PDEN(I+IFMT+IBPSCR) =0.
      PRESS(I+IFMT+IBPSCR)=AUX41(1,I)
  301 CONTINUE
      END IF


!>> 4) COMPUTE  WIRINGA, FIKS, AND FABROCINI'S EOS
      IDIM=IVAR
!                 **** ****      ****              ****
      CALL WFF123(AUX3,AUX4,IDIM,IPAN,SC,IVAR,KST1,AUX5)



      IF (ICCORR /= 0) THEN
! CORRECT EOS FOR CAUSALITY

      IF (KST1 == (-60)) THEN
      ESW = ESWUU
      EPW = EPWUU
      ENW = ENWUU
      ELSE IF (KST1 == (-70)) THEN
      ESW = ESWAU
      EPW = EPWAU
      ENW = ENWAU
      ELSE
      WRITE (6,'(///,'' _WFF88_  KST1 OUT OF RANGE'',
     +            //,'' ERROR STOP'')')
      STOP
             END IF

      DO 1 I=1,IPAN
      I1P   =I-1
    1 IF (AUX3(I) > ESW) GOTO 2
    2 CONTINUE
      IF (I1P == IPAN) THEN
      WRITE (6,'(///,'' _WFF88_  I1P=IPAN --> CAUSALITY CORRECTION'',
     +               '' FAILS'',//,'' ERROR STOP'')')
      STOP
      END IF

! COMPUTE DENSITY GRID
      ID1P = I100 - I1P
      XNA  = ENW  * (1. + .001)
      XNB  = 2.   * ENW
      ARGINT='WFF88: CAUSALITY CORRECTION   '
!                            ****      BARYON DENSITY GRID  N(I)
      CALL INTV(XNA,XNB,ID1P,AUX1,0,ARGINT)

! COMPUTE CAUSALITY CORRECTED PART OF THE EOS
      DO 3 I=1,ID1P
      X       = AUX1(I)/ENW
      PR      = (EPW-ESW + (EPW+ESW) * X * X ) / 2.
      AUX8(I) = PR
      EY      = ESW - EPW + PR
      AUX7(I) = EY
    3 CONTINUE

! COPY
      DO 4       I=1,ID1P
      AUX3(I+I1P) =  AUX7(I)
      AUX4(I+I1P) =  AUX8(I)
      AUX5(I+I1P) =  AUX1(I)
    4 CONTINUE

      IPAN = I100

! END OF CAUSALITY CORRECTION
                                 END IF



      IF (IBBPE > 0) THEN
      IV00=IBBP
      ELSE
      IV00=NAB
      END IF
      ITOTAL=IFMT+IBPSCR+IV00+IPAN
      IF (ITOTAL > IVAR) GOTO 9999
      DO 400 I=1,IPAN
      ENGD(I+IFMT+IBPSCR+IV00) =AUX3(I)
      PDEN(I+IFMT+IBPSCR+IV00) =AUX5(I)
      PRESS(I+IFMT+IBPSCR+IV00)=AUX4(I)
  400 CONTINUE

      IEP=ITOTAL
      ICOMB=-1


!>> CHECK ORDERING OF ENERGY DENSITY AND PRESSURE
!                              *****
      CALL SORTOR(ENGD,IEP,-10,IERSO)
      IERSO1=IERSO
!                               *****
      CALL SORTOR(PRESS,IEP,-10,IERSO)
      IERSO2=IERSO

      IF (IERSO2 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _BBP71_ PRESSURE NOT MONOTONOUSLY'',
     +       '' INCREASING'')')
      WRITE (6,505)
      END IF

      IF (IERSO1 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _BPS_ ENERGY DENSITY NOT MONOTONOUSLY'',
     +       '' INCREASING'',//)')
      DO 510 I=1,IEP
  510 WRITE (6,515) I,ENGD(I),PRESS(I)
  515 FORMAT (' I=',I4,2X,'E=',E13.6,2X,'P=',E13.6)
      STOP
      END IF
  505 FORMAT (//,'**************************************************',
     +        '************************',//)
      RETURN



 9999 CONTINUE
!__ ERROR STOP
      WRITE (6,'(//,'' DIMENSION OF AUXILIARY ARRAYS AUX3 AND AUX4 IS'',
     +          '' TOO SMALL!'',/,'' ITOTAL='',I4,'' IFMT='',I4,'' IB'',
     +          ''PSCR='',I4,'' IBBP='',I4,'' IPAN='',I4,'' IVAR='',I4,/
     +          )') ITOTAL,IFMT,IBPSCR,IBBP,IPAN,IVAR
      STOP

      END

!WFF123                 * *      ****               ***
      SUBROUTINE WFF123(E,P,IDIM,IOTCH,SC,IVAR,KST1,RHO)

!-----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATE THE EOS'S OF WIRINGA, FIKS, AND FABROCINI [WFF88].
!
!
!
! RETURN: (1) NUMBER OF PRESSURE AND ENERGY DENSITY VALUES (=IOTCH),
!
!         (2) EQUATION OF STATE, i.e., 
!             - E(I)    CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - P(I)    CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - RHO(I)  CONTAINS THE DENSITY OF PARTICLES PER 1/fm**3
!
!-----------------------------------------------------------------------

      REAL SC(4,IVAR),RHO(IVAR)

      PARAMETER (IP24=18,NP100=100,NPVAR=3*NP100)
      DIMENSION P(IDIM),E(IDIM)
      DIMENSION AP(IP24),AE(IP24),AR(IP24)
      DIMENSION RAV14U(IP24),EAV14U(IP24)
      DIMENSION RUV14U(IP24),EUV14U(IP24)
      DIMENSION RUV14T(IP24),EUV14T(IP24)


      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT

      CHARACTER * 30 EASS


!__ DENSITY   (in  1/fm^3)
!     MODEL: AV14+UVII
      DATA (RAV14U(I), I=1,IP24)/0.07,0.08,0.10,0.125,0.15,0.175,
     +0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.70,0.80,1.00,1.25,1.50/

!     MODEL: UV14+UVII
      DATA (RUV14U(I), I=1,IP24)/0.07,0.08,0.10,0.125,0.15,0.175,
     +0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.70,0.80,1.00,1.25,1.50/

!     MODEL: UV14+TNI
      DATA (RUV14T(I), I=1,IP24)/0.07,0.08,0.10,0.125,0.15,0.175,
     +0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.70,0.80,1.00,1.25,1.50/


!__ ENERGY PER NUCLEON   (in  MeV)
!     MODEL: AV14+UVII
      DATA (EAV14U(I), I=1,IP24)/7.35,7.94,8.97,10.18,11.43,12.74,
     +14.12,16.96,20.48,24.98,30.44,45.15,66.40,93.60,132.10,233.,
     +410.0,635.0/

!     MODEL: UV14+UVII
      DATA (EUV14U(I), I=1,IP24)/8.13,8.66,9.70,11.06,12.59,14.18,
     +15.92,20.25,25.78,32.60,40.72,61.95,90.2,126.2,170.5,291.1,
     +501.0,753.0/

!     MODEL: UV14+TNI
      DATA (EUV14T(I), I=1,IP24)/5.95,6.06,6.40,7.17,8.27,9.70,
     +11.55,16.29,22.19,28.94,36.60,56.00,79.2,106.1,135.5,200.9,
     +294.0,393.0/



      II8=IP24
      IF (II8 > IDIM) THEN
      WRITE (6,'(//,'' FATAL DIMENSIONING ERROR IN _WFF123_'',/,
     +           '' II8='',I4,'' <!  IDIM='',I4)') II8,IDIM
      STOP
      END IF


!  COPY
      DO 23 I=1,II8
      IF (KST1 == (-70)) THEN
         AR(I)=RAV14U(I)
         AE(I)=EAV14U(I)
         ELSE IF (KST1 == (-60)) THEN
         AR(I)=RUV14U(I)
         AE(I)=EUV14U(I)
         ELSE IF (KST1 == (-50)) THEN
         AR(I)=RUV14T(I)
         AE(I)=EUV14T(I)
         ELSE
         WRITE (6,'(//,1X,''KST1 NOT CORRECTLY INITIATLIZED'',/,
     +              '' FATAL ERROR IN  WFF123, KST1='',I6)') KST1
         STOP
      END IF
   23 CONTINUE



!  COMPUTE  D [E/N] / D [RHO]
      EASS='WFF123: CALL SELSPL 123       '
      DO 123 I=1,II8
      XXX=AR(I)
      RHO(I)=XXX
      CALL SELSPL(AR,AE,II8,SC,10,YYY,XXX,100,EASS)
      AP(I)= AR(I)**2 * YYY
  123 CONTINUE

      IWT=0
      DO 1000 I=1,II8
      E(I) = (AE(I) + AMNEUT) * AR(I)
      P(I) = AP(I)
      IF (IWT /= 0) THEN
      IF (I == 1 .AND. KST1 == (-70)) WRITE (6,444)
      IF (I == 1 .AND. KST1 == (-60)) WRITE (6,445)
      IF (I == 1 .AND. KST1 == (-50)) WRITE (6,446)
      WRITE (6,333) I,P(I),E(I)
      END IF
 1000 CONTINUE


  333 FORMAT (' I=',I4,' P=',E12.6,' E=',E12.6)
  444 FORMAT (/,' WIRINGA-FIKS-FABROCINI EOS:  AV14+UVII',//)
  445 FORMAT (/,' WIRINGA-FIKS-FABROCINI EOS:  UV14+UVII',//)
  446 FORMAT (/,' WIRINGA-FIKS-FABROCINI EOS:  UV14+TNI',//)

      IOTCH=II8


      RETURN

      END



!RD84                 *** *****
      SUBROUTINE RD84(IEP,ICOMB)


! ----------------------------------------------------------------------
!
!  PURPOSE:
!
!  COMPUTATION OF THE EQUATIONS OF STATE OF:
!
!       -  FEYNMAN-METROPOLIS-TELLER (1949)    \   TAKEN FROM
!       -  BAYM-PETHICK-SUTHERLAND (1971)      /   [BPS71]
!
!       -  NEGELE-VAUTHERIN (1973)            \
!       -  PANDHARIPANDE-PINES-SMITH (1976)    \   TAKEN FROM
!                                              /   [PPS76]
!       -  BETHE-JOHNSON (1974), MODEL  I     /
!
!
!  RETURN: ENGD  = ENERGY DENSITY (IN MeV/fm**3),
!          PRESS = PRESSURE (IN MeV/fm**3),
!          IEP   = NUMBER OF P VS. E  VALUES
!
! ----------------------------------------------------------------------


      REAL * 8 EE18,EEM18,EE14,EE34,EE44
      real * 8 ee03,ee13,ee15

      PARAMETER (NP100=100,NPVAR=3*NP100)
      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)

      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/ExpoMaxMin/EXMAX,EXMIN
      COMMON/B12345/B1,B2,B3,B4,B5
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST

      CHARACTER * 30 ARGINT

      DATA AL1/236./,AL2/364./



!>> 1) COMPUTE FEYNMAN-METROPOLIS-TELLER EOS
      IDIM=IVAR
!                **** **** ****      ****
      CALL FMT49(AUX3,AUX4,AUX1,IDIM,IFMT)
      DO 10 I=1,IFMT
      ENGD(I) =AUX3(I)
      PDEN(I) =AUX1(I)
   10 PRESS(I)=AUX4(I)


!>> 2) COMPUTE BETHE-PETHICK-SUTHERLAND EOS (CRUST REGION!)
      IDIM=IVAR
!                **** **** ****      ******
      CALL BPSCR(AUX3,AUX4,AUX1,IDIM,IBPSCR)
      ITOTAL=IFMT+IBPSCR
      IF (ITOTAL > IVAR) GOTO 9999
      DO 20 I=1,IBPSCR
      ENGD(I+IFMT) =AUX3(I)
      PDEN(I+IFMT) =AUX1(I)
   20 PRESS(I+IFMT)=AUX4(I)


!>> 3) COMPUTE NEGELE-VAUTHERIN EOS   A N D   EOS AS OBTAINED BY GRADUALLY
!      AVERAGING THE NV NUCLEAR ENERGIES WITH THE NEUTRON MATTER ENERGIES
!      CALCULATED IN THE BJ(I) MODEL (SEE [PPS76])
!      MASS DENSITY:  3.7 E11 < M.D. < 1 E13   (NV),
!                     1   E13 < M.D. < 1 E14   (PPS-AVERAGING)
      RMD1=ENGD(ITOTAL)
      RMD2=4.0E14/UF6
      IDIM=IVAR
      CALL BJNVAV(RMD1,RMD2,AUX3,AUX4,IDIM,IOBJ)
      ITOTAL=IFMT+IBPSCR+IOBJ
      IF (ITOTAL > IVAR) GOTO 9999
      DO 30 I=1,IOBJ
      ENGD(I+IFMT+IBPSCR) =AUX3(I)
      PDEN(I+IFMT+IBPSCR) =0.
   30 PRESS(I+IFMT+IBPSCR)=AUX4(I)


!>> 4) COMPUTE BETHE-JOHNSON EQUATION OF STATE (MODEL I)
      BORD1=0.15
      BORD2=2.50
      NAB=I100
      ARGINT='RD84: ENERGY DENSITY GRID     '
      CALL INTV(BORD1,BORD2,NAB,AUX3,0,ARGINT)

      DO 100 I=1,NAB
      R=AUX3(I)
      E= AL1 * R**1.54 + AMNEUT
      AUX3(I)= R * E

      IF (AUX3(1)  <=  ENGD(ITOTAL)) THEN
      WRITE (6,'(///,'' COMBINATION OF BJ(I) EOS WITH PPS76  FAILED'',/,
     +           '' CHOOSE A VALUE LARGER THAN  '',E14.8,'' MeV/fm**3'',
     +           ///)') ENGD(ITOTAL)
      STOP 'ERROR STOP IN  _RD84_'
      END IF

      AUX4(I)= AL2 * R**2.54
  100 CONTINUE

      DO 200 I=1,NAB
      ENGD(I+ITOTAL) =AUX3(I)
      PDEN(I+ITOTAL) =0.
  200 PRESS(I+ITOTAL)=AUX4(I)

      IEP=ITOTAL+NAB
      ICOMB=-1


!>>> CHECK ORDERING OF ENERGY DENSITY AND PRESSURE
!                              *****
      CALL SORTOR(ENGD,IEP,-10,IERSO)
      IERSO1=IERSO
!                               *****
      CALL SORTOR(PRESS,IEP,-10,IERSO)
      IERSO2=IERSO

      IF (IERSO2 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _RD84_ PRESSURE NOT MONOTONOUSLY'',
     +       '' INCREASING'')')
      WRITE (6,505)
      END IF

      IF (IERSO1 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _RD84_ ENERGY DENSITY NOT MONOTONOUSLY'',
     +       '' INCREASING'',//)')
      DO 510 I=1,IEP
  510 WRITE (6,515) I,ENGD(I),PRESS(I)
  515 FORMAT (' I=',I4,2X,'E=',E13.6,2X,'P=',E13.6)
      STOP 'ERROR STOP'
      END IF
  505 FORMAT (//,'**************************************************',
     +        '************************',//)


      RETURN


 9999 CONTINUE
!__ ERROR STOP
      WRITE (6,'(//,'' DIMENSION OF AUXILIARY ARRAYS AUX3 AND AUX4 IS'',
     +          '' TOO SMALL!'',/,'' ITOTAL='',I4,'' IFMT='',I4,'' IB'',
     +          ''PSCR='',I4,'' IVAR='',I4,/)') ITOTAL,IFMT,IBPSCR,IVAR
      STOP 'ERROR STOP IN _RD84_'

      END

!BJNVAV                           * *      *****
      SUBROUTINE BJNVAV(RMD1,RMD2,E,P,IDIM,IOTCH)


!-----------------------------------------------------------------------
!
! PURPOSE:
!         -  CALCULATE THE EOS OF NEGELE AND VAUTHERIN (FOR MASS DENSI-
!            TIES:  3.7 E11 < M.D. < E13) ([PPS76]),
!         -  CALCULATED THE EOS OF [PPS76], i.e., GRADUALLY AVERAGING
!            OF THE NV NUCLEAR ENERGIES WITH THE NEUTRON MATTER ENERGIES
!            CALCULATED IN THE BJ MODEL (E13 < M.D. < E14)
!         -  FOR M.D. > E14  APPROXIMATE  P(E) BY PURE NEUTRON MATTER
!            RESULTS ([PPS76]).
!
! INPUT: RMD1=MINIMUM ENERGY DENSITY (IN MeV/fm**3),
!        RMD2=MAXIMUM ENERGY DENSITY (IN MeV/fm**3); THE EOS WILL BE
!        CALCULATED ONLY FOR ENERGY DENSITIES  RMD1 < E.D. < RMD2!
!
! RETURN:  NUMBER OF PRESSURE RESPECTIVELY ENERGY DENSITY VALUES (IOTCH)
!
!          AND EQUATION OF STATE, i.e., 
!             - E(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!             - P(I)  CONTAINS THE ENERGY DENSITY IN  MeV/fm**3
!
!-----------------------------------------------------------------------


      REAL * 8 AP

      PARAMETER (IP19=19)
      DIMENSION P(IDIM),E(IDIM),AP(IP19),AE(IP19)


      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10



!__ PRESSURE   (in  dynes/cm^2)
      DATA (AP(I), I=1,IP19)/
     +1.402E29,3.134E29,7.157E29,1.036E30,1.257E30,2.122E30,3.780E30,
     +8.527E30,1.162E31,3.262E31,9.407E31,4.888E32,8.584E32,3.879E33,
     +1.549E34,7.009E34,1.167E35,5.434E35,1.896E36/


!__ MASS DENSITY   (in  g/cm^3)
      DATA (AE(I), I=1,IP19)/
     +1.000E11,2.000E11,4.000E11,8.000E11,1.000E12,2.000E12,4.000E12,
     +8.000E12,1.000E13,2.000E13,4.000E13,8.000E13,1.000E14,2.000E14,
     +4.000E14,8.000E14,1.000E15,2.000E15,4.000E15/



      II8=IP19
      IF (II8 > IDIM) THEN
      WRITE (6,'(//,'' FATAL DIMENSIONING ERROR IN _BJNVAV_'',/,
     +           '' II8='',I4,'' <!  IDIM='',I4)') II8,IDIM
      STOP
      ELSE IF (RMD1 >= RMD2) THEN
      WRITE (6,'(//,'' FATAL ARGUMENT ERROR IN _BJNVAV_'',/,
     +           '' RMD1='',E12.6,'' MeV/fm**3;  RMD2='',E12.6,
     +           '' MeV/fm**3'')') RMD1,RMD2
      STOP
      END IF


      IWT=0
      DO 1000 I=1,II8

!  CONVERT  dynes/cm**2  -->  MeV/fm**3
!  CONVERT  g/cm**3      -->  MeV/fm**3
      P(I) = SNGL ( AP(I) ) / UF4
      E(I) = AE(I) / UF6
      IF (IWT /= 0) THEN
      IF (I == 1) WRITE (6,444)
      WRITE (6,333) I,P(I),E(I)
      END IF

 1000 CONTINUE

  333 FORMAT (' I=',I4,' P=',E12.6,' E=',E12.6)
  444 FORMAT (/,' PANDHARIPANDE-PINES-SMITH EQUATION OF STATE',//)


      IF ((RMD1 <= E(1)) .AND. (RMD2 >= E(II8))) THEN
      RETURN

      ELSE

      DO 222 I=1,II8
      EEE=E(I)
      PPP=P(I)

      IF (EEE <= RMD1) THEN
      KK=I
      ELSE IF ((EEE > RMD1) .AND. (EEE < RMD2)) THEN
      E(I-KK)=EEE
      P(I-KK)=PPP
      IOTCH=I-KK
      END IF

  222 CONTINUE

      END IF


      RETURN

      END


!PTROP                 *** *****
      SUBROUTINE PTROP(IEP,ICOMB)


! ----------------------------------------------------------------------
!
!  PURPOSE:
!
!  COMPUTATION OF THE n=3/2 RELATIVISTIC POLYTROPE EQUATION OF STATE OF
!  TROOPER (1965) (USE BY AND HARTLE (1973)).
!  HERE: USE ONLY FOR THE PURPOSE OF TEST CALCULATIONS.
!
!
!
!  RETURN: ENGD  = ENERGY DENSITY (IN MeV/fm**3),
!          PRESS = PRESSURE (IN MeV/fm**3),
!          IEP   = NUMBER OF  P(E)  VALUES.
!
! ----------------------------------------------------------------------


      REAL * 8 EE18,EEM18,EE14,EE34,EE44
      real * 8 ee03,ee13,ee15

      PARAMETER (NP100=100,NPVAR=3*NP100)
      COMMON/EP1/ENGD(NPVAR),PRESS(NPVAR),EINTNL(NPVAR),EINTH(NPVAR),
     +           PDEN(NPVAR)
      COMMON/AX1/AUX1(NPVAR),AUX3(NPVAR),AUX4(NPVAR),AUX31(1,NPVAR),
     +           AUX41(1,NPVAR)

      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/ExpoMaxMin/EXMAX,EXMIN
      COMMON/B12345/B1,B2,B3,B4,B5
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST

      CHARACTER * 30 ARGINT



!>> CONSTANTS OF THE CALCULATION
      VV15=1.5
      VV68=68.09138
      E35 =3./5.


      IEP=IVAR
      APRESS=1.0E-24
      BPRESS=3.0E+03

!   COMPUTE PRESSURE  P,  [P]=MeV/fm^3
      AXA=ALOG10(APRESS)
      BXB=ALOG10(BPRESS)
      ARGINT='PTROP: ENERGY DENSITY GRID    '
      CALL INTV(AXA,BXB,IEP,AUX3,0,ARGINT)
      DO 100 I=1,IEP
      PRESS(I)=10.**AUX3(I)
  100 CONTINUE

!   COMPUTE ENERGY DENSITY  E,  [E]=MeV/fm^3
      DO 200 I=1,IEP
      ENGD(I) = VV15 * PRESS(I) + VV68 * PRESS(I)**E35
  200 CONTINUE

!   NO CALCULATION OF THE PARTICLE DENSITY  n,  [n]=1/fm
      DO 210 I=1,IEP
      PDEN(I)=0.
  210 CONTINUE

      ICOMB=-1


!>>  CHECK ORDERING OF ENERGY DENSITY AND PRESSURE
!                              *****
      CALL SORTOR(ENGD,IEP,-10,IERSO)
      IERSO1=IERSO
!                               *****
      CALL SORTOR(PRESS,IEP,-10,IERSO)
      IERSO2=IERSO

      IF (IERSO2 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _PTROP_ PRESSURE NOT MONOTONOUSLY'',
     +       '' INCREASING'')')
      WRITE (6,505)
      END IF

      IF (IERSO1 /= 0) THEN
      WRITE (6,505)
      WRITE (6,'('' _PTROP_ ENERGY DENSITY NOT MONOTONOUSLY'',
     +       '' INCREASING'',//)')
      DO 510 I=1,IEP
  510 WRITE (6,515) I,ENGD(I),PRESS(I)
  515 FORMAT (' I=',I4,2X,'E=',E13.6,2X,'P=',E13.6)
      STOP 'ERROR STOP'
      END IF
  505 FORMAT (//,'**************************************************',
     +        '************************',//)


      RETURN

      END

!SDWF
!                                                             **** ****
      SUBROUTINE SDWF  (EDDY,PRESS,IEP,EC,AUX3,AUX4,NHU,IFAIL,RRNS,RQOZ,
!                       **               ***
     +                  WI,C,AUX31,AUX41,RTH,NI,MLF,IREM,XIPEOS,EDRIP,EQ
     +                  DRIP,PDRIP,IDRIA,IDRIB,BORD0,BORD1,PDEN,EJGCM3,M
     +                  LFLL,MLFUL,R5KM,IQ2020,SSLT,SSGT,C8,PCUT,PCUTS,I
     +                  RUKOV,POPMFM,BAGCT,MLG)


!-----------------------------------------------------------------------
!
! PURPOSE:
!          CONSTRUCT STRANGE DWARF STAR MODELS WITH VANISHINGLY SMALL 
!          STRANGE CORE RADII.
!
! INPUT:   R_core (in meters)
!
! RETURN:
!          RRNS=R_NS  (in meters)
!          RQOZ=M_NS / M_sun
!          WI=LOG(MOMENT OF INERTIA OF THE STAR [g cm^2])
!
! ERROR RETURN: IFAIL<0 
!-----------------------------------------------------------------------

      PARAMETER (m651=160501)


      REAL * 8 DGRAV,A,AR0,EE18,EE55,P16,GCM2KM,EEM18,F896,F178,F288
      REAL * 8 RMSUN,RMNS,GRAV,GRAVK,GRAVKM,EE14,EE34,EE44
      real * 8 ee03,ee13,ee15
      REAL * 8 C8(4,NHU),DMN,F25116(M651)

      DIMENSION EDDY(NHU),PRESS(NHU),AUX3(NHU),AUX4(NHU),C(4,NHU)
      DIMENSION AUX31(1,NHU),AUX41(1,NHU),RTH(0:NI,4),PDEN(NHU)


      PARAMETER (KP1=1,KPIRF6=6)
      DIMENSION PIP(KP1),EIP(KP1),PDENS(KP1),YPD(KP1)

      PARAMETER (N27=50,NPOVDP=4)
      COMMON/TABCCC/TOVDP(N27,NPOVDP)

      COMMON/ExpoMaxMin/EXMAX,EXMIN
      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      COMMON/Exponent55/EE55
      common/ConversionFactor30/ufe30
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/GRAVITATIONAL/GRAVK,GRAV,GRAVKM,GCM2KM
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/MISHRA/IM05
      COMMON/SOLARPROPERTIES/RMSUN,RSUN,RSUNKM
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/MISCINT2/J651,J650,IGIVMX,IVONR,KC20
      COMMON/RADIALSCALINGFACTOR/R0R0
      COMMON/NUCLEARMATTERDENSITY/ENGNM0
      COMMON/IPLOT2/IOUTD4,IPOLXX,ICMI99

      COMMON/SPL1/ST(M651),YOUT(M651),XSTROB(M651),YSTROB(M651)
      COMMON/CNSTA1/F251(M651),FOUT(M651),FST(M651),FST251(M651),
     +              FENG1(M651),X251(M651),Y251(M651),F251PM(M651),
     +              FDRIP(M651)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA3/AUX61(1,M651),AUX6(M651),DPDRAX(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)
      COMMON/BN01/AXBN1(M651),AXBN2(M651),PDENEX(M651),EINTHX(M651),
     +            DEINDP(M651),DEHDRH(M651),PDENT(M651)

      PARAMETER (NTHOV1=12)
      COMMON/TABSPH/THOV1(N27,NTHOV1)

      parameter (kpc20=20,kpc5=6)
      common/tc/tcore(kpc20,kpc5),r_core(kpc20)


      CHARACTER * 25 XIPEOS
      CHARACTER * 30 ARGINT
      CHARACTER * 30 EASS
      CHARACTER *  5 ERK


      DATA FMKM/1.E18/,FMEVFM/5.61012/,F178/1.782531E15/
      DATA F896/8.961946E-16/,F288/2.883716E-06/

      NAMELIST/NAME1/RMNS,RMSUN,FOUT,FST251,J650,J651 



! ------------------------- AUXILIARY FUNCTIONS -------------------------

!   NEWTONIAN LIMIT OF OV
      FH(DEN) = - DGRAV * DEN * DMN / DRN2

! -----------------------------------------------------------------------



      WRITE (6,'(//,''  _SDWF_  ->  N O N - ROTATING STAR CAL'',
     +         ''CULATION;  MLF='',I3,'';  MLG='',I3)') MLF,MLG
      WRITE (6,'(''                  Range of MLF: '',I2,'' <= MLF <='',
     +           I2)') MLFLL,MLFUL
      RX=BORD1*R0R0/EE18
      WRITE (6,'(''  NOTE: e_c='',E12.6,'' MeV/fm^3;  R_max='',E12.6,
     +           '' km;  R_core='',E12.6,'' m'')') EC,RX,BORD0

      IFAIL=10
      IF (IEP > NHU) THEN
      IFAIL=-10
      WRITE (6,'('' ** _SDWF_ FIELD LENGTH OUT OF RANGE ** '')')
      WRITE (6,'(5X,'' IEP='',I4,5X,''NHU='',I4,///)') IEP,NHU
      RETURN
      END IF

      ISY1  =0
      LSTOP = 10
      IRF6  =KPIRF6
      IRF6M1=IRF6-1

      ERK='EULER'
      IF (IRUKOV /= 0) ERK='RU-KU'

      POPMFM=PDRIP
      EJGCM3=EC*UF6


!__ CONSTANTS
      IF (IKP1 /= KP1) THEN
      WRITE (6,'(//,'' ****  _SDWF_  IKP1 # KP1  ****'',//,
     +              '' IKP1='',I4,3X,'' KP1='',I4)') IKP1,KP1
      IFAIL=10
      RETURN
      END IF
      AR0   =R0R0
      A     =EC*AR0**3/(RMSUN*ufe30)
      P16   =AR0/EE18
      DGRAV =(GRAV/ufe30)*(RMSUN*ufe30)/(UF1*AR0)
      RMETER=R0R0/EE15
      DPDRIP=PDRIP/EC
      VBAG  =4.*BAGCT**4/UF1**3


!__ INITIALIZE
      DO 8720 I=1,J651
      AXBN1(I) =0.
      AXBN2(I) =0.
 8720 CONTINUE
      RLASPH   =0.
      RLACOR   =0.
      RLACRU   =0.




!__ ENERGY DENSITY AND PRESSURE (DIMENSIONLESS)
!   Note: only eos of white dwarf star matter is needed
      DO 10    I=1,IEP
      AUX3(I)   =EDDY(I)/EC
      AUX4(I)   =PRESS(I)/EC
      AUX31(1,I)=AUX3(I)
      AUX41(1,I)=AUX4(I)
!   BARYON DENSITY (1/fm**3)
      AUX61(1,I)=PDEN(I)
   10 CONTINUE

      DEC    =EC/EC

      FMP = 1.E-04
      IF ((DEC > 1.+FMP).OR.(DEC < 1.-FMP)) THEN
      WRITE (6,'(///,'' FATAL ERROR IN _SDWF_'',/,'' FAILED '',
     +           ''TO DETERMINE CENTRAL DENSITY'',/,'' EC='',E12
     +           .6,''  EJGCM3='',E12.6)') EC,EJGCM3
      STOP '_SDWF_ DRIP PRESSURE'
      END IF

!   terminate calculation if pressure  p < p_crust(Fe)
      DPCUTS =PCUTS/EC


!__ COMPUTE GRID FOR RADIAL INTEGRATION,  R=[0,1]

      WRITE (6,'(/,2x,''begin computation of radial grid:'')')
      BORD00=BORD0/RMETER
      R5KM  =BORD00*RMETER/1.E03
      BORD11=BORD1
      ARGINT='SDWF: RADIAL GRID - FST251    '
!                                  ******
      CALL INTV(BORD00,BORD11,J651,FST251,0,ARGINT)
      SSGT=(FST251(2)-FST251(1))*(R0R0/EE18)*1.E03
      SSLT=SSGT
      WRITE (6,'(2X,''step size for radial integration: '',e12.6,'' m'',
     +           /)') ssgt




      IOPTIP=IPOLXX
      IF (IOPTIP < 0) THEN
      XIPEOS='LOGARITHMIC INTERPOLATION'
      ELSE IF (IOPTIP == 1) THEN
      XIPEOS='LINEAR INTERPOLATION     '
      ELSE IF (IOPTIP == 2) THEN
      XIPEOS='POLYNOMIAL INTERPOLATION '
      ELSE IF (IOPTIP == 3) THEN
      XIPEOS='SPLINE INTERPOLATION     '
      ELSE
      XIPEOS=' ***** NOT SPECIFIED ****'
      END IF


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                 SIMULTANEOUS INTEGRATION OF OV EQUATIONS
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


!__ INITIAL VALUES AT R_core AND R_core+h

!     PRESSURE AT R_core
!                                    ***           RETURN: PIP=P_core
      CALL SELINT(AUX3,AUX41,1,1,IEP,PIP,DEC,1)
      FOUT(1)=PIP(IKP1)
      PN     =PIP(IKP1)
      PDRIP  =FOUT(1)*EC

!*** DRMC16: 2 lines have been moved down
      
!  MASS OF STRANGE CORE (in units of M_sun)
      DRMC16      =VPI * BORD0**3 * (3.*PDRIP+VBAG) * F896 / 3.
      DRMC        =DRMC16
      tcore(mlg,2)=DRMC
      PRINT*,' DRMC=',DRMC,' M_sun'

! DRMC16 is not initialized -> moved these 2 lines from !***DRMC16 to here
      DPDR00 =-A*DGRAV*DEC*DRMC16/FST251(1)**2
      FOUT(2)=FOUT(1) + (FST251(2)-FST251(1))*DPDR00
      

!     DENSITY AT R_core
      EIP(1)  =DEC
      FENG1(1)=DEC

!     MASS AT R_core
      F25116(1) =DRMC16
      F25116(2) =VPI*A*DEC * (FST251(2)**3-FST251(1)**3) / 3.
      F25116(2) =F25116(1) + F25116(2)

!     CRUST MASS
      FDRIP(1)=0.
      FDRIP(2)=VPI*A*DEC * (FST251(2)**3-FST251(1)**3) / 3.


!     BARYON DENSITY, BPS AT R_core
!                                    *****
      CALL SELINT(AUX3,AUX61,1,1,IEP,PDENS,DEC,1)
      PDENEX(1)    =PDENS(IKP1)
      THOV1(MLF,10)=PDENS(IKP1)
      AXBN2(2)=VPI*A* PDENEX(1) * (FST251(2)**3-FST251(1)**3) / 3.

      IF (EJGCM3 /= 0.)
     +PRINT*,' Note: e_join < e_drip (hadronic crust below n-drip)'
      IF (EJGCM3 == 0.)
     +PRINT*,' Note: e_join = e_drip (hadronic crust at n-drip)'
      PRINT*,'       P_join =',PDRIP,' MeV/fm^3'



! ******************************************************************
!           START INTEGRATING HYDROSTATIC EQUILIBRIUM EQS
! ******************************************************************
      INT997=997
      NXY3  =3
      INTCT1=NXY3+INT997

      DO 1000 N=NXY3,J650

      IREM     =N-1
      DMN      =F25116(N-1)
      DRN      =FST251(N-1)
      DRN2     =DRN*DRN
      DRN3     =DRN2*DRN
      DPN      =FOUT(N-1)
      PIP(IKP1)=DPN
      APRMEV   =DPN*EC


      IF (PIP(IKP1) <= AUX4(1)) THEN
      EIP(IKP1)=AUX3(1)


      ELSE



      IF (XIPEOS == 'SPLINE INTERPOLATION     ') THEN

      IF (APRMEV > PCUT) THEN

      EASS='SDWF: CALL SELSPL AAA         '
!                                    ****
      CALL SELSPL(AUX4,AUX3,IEP,C,10,YYYY,DPN,10,EASS)
      EIP(IKP1)=YYYY
      EASS='SDWF:   CALL SELSPL BBB       '
!                                    ****
      CALL SELSPL(AUX4,PDEN,IEP,C,10,YYYY,DPN,10,EASS)
      YPD(IKP1)=YYYY
      ELSE

!                                    ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,-10)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                    ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,  1)

! no access (logarithmic interpolation works better)
!      EASS='SDWF: CALL SELSPL AAADDD      '
!                                      ****
!      CALL SELSPL8(AUX4,AUX3,IEP,C8,10,YYYY,DPN,10,EASS)
!      EIP(IKP1)=YYYY
!      EASS='SDWF:   CALL SELSPL BBBDDD    '
!                                      ****
!      CALL SELSPL8(AUX4,PDEN,IEP,C8,10,YYYY,DPN,10,EASS)
!      YPD(IKP1)=YYYY
! end of no access

      END IF

      ELSE
!                                    ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,IOPTIP)
!                                    ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,IOPTIP)
      END IF


              END IF



      DEN       =EIP(IKP1)
      FENG1(N-1)=DEN

!__ n(p) (1/fm^3)
      PDENEX(N-1)=YPD(IKP1)


      HRK       = FST251(N)-FST251(N-1)
      HRKH      = HRK/2.


      IF (IRUKOV == 0) THEN


!__ integrate hydrostatic equilibrium equations via Euler method

      DPDR      = FH(DEN)
      DPDRAX(N) = DPDR

!............................................
!>>>> P(N)
      PNP1  = PN + HRK * DPDR
!............................................


      ELSE

!__ preliminary computations for solving the equations of hydrostatic
!   equilibrium via the fourth-order Runge-Kutta method
      DPDR      = FH(DEN)
      DPDRAX(N) = DPDR
      XK1       = HRK * DPDR

      DRN       = FST251(N-1)+HRKH
      DRN2      = DRN*DRN
      DRN3      = DRN2*DRN
      DPK1      = DPN+XK1/2.
      DPDR2     = FH(DEN)
      XK2       = HRK * DPDR2

      DPK2      = DPN+XK2/2.
      DPDR3     = FH(DEN)
      XK3       = HRK * DPDR3

      DRN       = FST251(N-1)+HRK
      DRN2      = DRN*DRN
      DRN3      = DRN2*DRN
      DPK3      = DPN+XK3
      DPDR4     = FH(DEN)
      XK4       = HRK * DPDR4

!................................................................
!>>>> P(N)
!     Compute advanced pressure via fourth order Runge-Kutta
      PNP1 = PN + (XK1 + XK4)/6. + (XK2 + XK3)/3.
!................................................................


      DRN       = FST251(N-1)
      DRN2      = DRN*DRN
      DRN3      = DRN2*DRN

      END IF



      PIP(IKP1)=PNP1
      FOUT(N)  =PNP1
      PN       =PNP1


      IF (N == INTCT1) THEN
      RKM1=FST251(N)*R0R0/EE18
      PMF3=PN*EC
      PRINT*,' N=',N,'  P=',PMF3,' MeV/fm^3','  r=',RKM1,' km','  ',ERK
      INTCT1=INTCT1+1000
      END IF



      PNLLIM=DPCUTS

      IF (PN <= PNLLIM) THEN
          THOV1(MLF,11)   =   PN*EC
          IF(PN > 0.) PN = - PNLLIM
          FOUT(N)         =   PN
          GO TO 2000
      END IF




!__ MASSES  (VIA SOLVING THE DIFFERENTIAL EQUATION):
!     GRAVITATIONAL MASS
      F25116(N)=F25116(N-1)+HRK*VPI*A*DRN2*DEN

!     PROPER STAR MASS
      F251PM(N)=0.

!     MASS OF THE STAR'S CRUST (DENSITY < DRIP DENSITY)
      FDRIP(N)=FDRIP(N-1)+HRK*VPI*A*DRN2*DEN

!     BARYON NUMBER of CRUST
      AXBN2(N)=AXBN2(N-1) + HRK*DRN2*PDENEX(N-1)
      ACRUST  =AXBN2(N)



 1000 CONTINUE
! ******************************************************************
!                        END OF INTEGRATION
! ******************************************************************

! NO ACCESS WRITE (6,'(5X,'' ** _SDWF_ ERROR DETECTED  **'',////)')
! NO ACCESS WRITE (6,NAME1)
      IFAIL=-20
      RRNS =-1.
      RQOZ =-1.
      RETURN



 2000 CONTINUE



!__ STAR'S RADIUS

      WRITE (6,'(''  surface of white dwarf/hyperon/hybrid/quark s'',
     +           ''tar encountered: irem='',I6,/)') irem



!   DETERMINE VALUE OF THE STAR'S RADIUS (in km)
      K44  =5
      KHIGH=IREM+1
      KLOW =KHIGH-K44

      DO 581      I=KLOW,KHIGH
      IRUN         =I-KLOW+1
      AUX3(IRUN)   =FST251(I)
      AUX31(1,IRUN)=FOUT(I)
  581 CONTINUE

      IDIM=KHIGH-KLOW+1
      IOUT=40
      IYY =1


      WRITE (6,'(''  --> Computing stars radius'',//)')
!                                                        ****
      CALL ZERO22(AUX3,AUX31,IYY,IDIM,AUX4,AUX41,ST,IOUT,YOUT)
      FST251(IREM)=YOUT(IYY)

      RRNS        =FST251(IREM)*RMETER
      tcore(mlg,4)=rrns/1.e03      
      RQOZ        =F25116(IREM)
      tcore(mlg,1)=RQOZ
      RMNS  =RQOZ*RMSUN*ufe30/5.61012E29



!__ COMPUTE NUMBER OF BARYONS CONTAINED IN THE STAR'S CRUST AND CORE
!   STRANGE CORE
      XN_C=4. * (PDRIP+BAGCT**4/UF1**3) / 3.
      XN_C=XN_C**0.75 / (SQRT(PI) * UF13**0.25)
      RLACOR=ALOG10(VPI/3.) + 3.*ALOG10(BORD0) + ALOG10(XN_C) + 45.
!   HADRONIC CRUST
      RLACRU=3.*DLOG10(AR0) + ALOG10(VPI*ACRUST)
      THOV1(MLF,8)=RLACOR
      THOV1(MLF,9)=RLACRU
      tcore(mlg,5)=RLACOR
      tcore(mlg,6)=RLACRU



!     AMOUNT OF THE STAR'S MASS BELOW NEUTRN DRIP
      RMDRIS      =FDRIP(IREM)*RMSUN*EE55
      TOVDP(MLF,2)=RMDRIS
      TOVDP(MLF,3)=FDRIP(IREM)
      TOVDP(MLF,4)=F251PM(IREM)
      tcore(mlg,3)=FDRIP(IREM)
!     R_core (in km)
      TOVDP(MLF,1)=BORD00*RMETER/1.E03


      DO 481 I=0,NI
  481 RTH(I,1)=RRNS
      FOUT(IREM)=(FOUT(IREM-1)+FOUT(IREM-2))/10.





! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                   MOMENT OF INERTIA OF THE STAR 
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      WI=0.
      IF (ICMI99 < 0) GO TO 6060

  800 DO 600 J=2,IREM
      DRN      =FST251(J)
      AUX8(J-1)=DRN
      DRN2     =DRN*DRN
      DRN4     =DRN2*DRN2
      DEN      =FENG1(J)
  600 FST(J-1) =DRN4*DEN 
      IREMM1=IREM-1

!__ PERFORM INTEGRATION
!
!                               R        4
!                   TM = const  I dr' (r') e(r')
!                               0                         

      EASS='_SDWF_ CALL SELSPL FFF      '
!                                              **
      CALL IGRN(AUX8,FST,IREMM1,X251,Y251,IREM,TM,C,NHU,10,1,AUX51,EASS)
!     RETURN: TM

!     I_core * const
      XICORE=BORD00**5*FENG1(1)/5.   

!__ COMPUTE  LOG_10 (I_total  (g cm^2) )
      TM=TM+XICORE
      WI= ALOG10(2.*VPI*EC/3.)+5.*ALOG10(R0R0)+ALOG10(TM)-ALOG10(FMEVFM)
     +   -52.



 6060 CONTINUE

      IEPRES=0
      IF (IEPRES /= 0) THEN
         OPEN (UNIT=66,FILE='USER.E_OUT',STATUS='unknown')
         OPEN (UNIT=67,FILE='USER.P_OUT',STATUS='unknown')
         DO 3799 I=1,IREM,10
         EU=FENG1(I)*EC/ENGNM0
         PU=FOUT(I)*EC
         P16=AR0/EE18
         RU=FST251(I)*P16
         WRITE (66,610) RU,EU
         WRITE (67,610) RU,PU
 3799    CONTINUE
  610    FORMAT (2X,E14.7,4X,E14.7)

         CLOSE ( 66 )
                      CLOSE ( 67 )

      END IF


      RETURN

      END



!OVNROT
      SUBROUTINE OVNROT(EDDY,PRESS,IEP,EC,AUX3,AUX4,NHU,IFAIL,RRNS,RQOZ,
     +                  WI,C,AUX31,AUX41,RTH,NI,MLF,IREM,XIPEOS,EDRIP,EQ
     +                  DRIP,PDRIP,IDRIA,IDRIB,BORD0,BORD1,PDEN,EJGCM3,M
     +                  LFLL,MLFUL,R5KM,IQ2020,SSLT,SSGT,C8,PCUT,PCUTS,I
     +                  RUKOV,POPMFM)


!-----------------------------------------------------------------------
!
! PURPOSE:
!          SOLVE THE OPPENHEIMER-VOLKOFF (OV) EQUATIONS OF A SPHERICAL,
!          NON-ROTATING STAR. GET RADIUS, GRAVITATIONAL MASS, PROPER
!          MASS, AND MOMENT OF INERTIA (FRAGE DRAGGING NEGLECTED) OF
!          THE STAR MODEL.
!          NOTE: START AT A GIVEN VALUE FOR THE CENTRAL ENERGY DENSITY
!                AND INTEGRATE THE OV EQUATIONS UNTILL THE PRESSURE
!                VANISHES, P(R)=0; THIS DEFINES THE RADIUS OF THE STAR.
!                THE STAR'S MASS FOLLOWS BY A SIMPLE VOLUME INTEGRATION
!                OVER THE ENERGY DENSITY FUNCTION.
!
! INPUT:
!        EOS P(EPSILON)  (IN MeV/fm**3): EDDY,PRESS
!        CENTRAL ENERGY DENSITY EC=EPSILON_C  (IN MeV/fm**3)
!        (150 MeV/fm**3 < E_C < 5000 MeV/fm**3)
!
!
! CONSTANTS:
!        GRAVITATIONAL CONSTANT: GRAV (IN fm**2), GRAVK (IN fm/MeV)
!        SCALING FACTOR: R0R0 (IN fm)
!        MASS OF THE SUN RMSUN: (IN MeV)
!
!
! RETURN:
!         RRNS=R_NS  (IN METERS)
!         RTH=R_NS   (0<=THETA<=#/2)
!         RQOZ=M_NS / M_SUN
!         WI=LOG(MOMENT OF INERTIA OF THE STAR (g cm^2))
!
!         IFAIL<0  <=>  ERROR RETURN
!-----------------------------------------------------------------------


      REAL * 8 DGRAV,A,AR0,EE18,EE55,P16,GCM2KM,EEM18,EE14,EE34
      REAL * 8 RMSUN,RMNS,GRAV,GRAVK,GRAVKM,EE44
      real * 8 ee03,ee13,ee15
      REAL * 8 C8(4,NHU)

      DIMENSION EDDY(NHU),PRESS(NHU),AUX3(NHU),AUX4(NHU),C(4,NHU)
      DIMENSION AUX31(1,NHU),AUX41(1,NHU),RTH(0:NI,4),PDEN(NHU)

      PARAMETER (KPIRF6=6)
      REAL A3(KPIRF6),A4(KPIRF6),P4(KPIRF6),A31(1,KPIRF6),A61(1,KPIRF6)

      PARAMETER (KP1=1)
      DIMENSION PIP(KP1),EIP(KP1),PDENS(KP1),YPD(KP1)

      PARAMETER (N27=50,NPOVDP=4)
      COMMON/TABCCC/TOVDP(N27,NPOVDP)

      COMMON/ExpoMaxMin/EXMAX,EXMIN
      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      COMMON/Exponent55/EE55
      common/ConversionFactor30/ufe30
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/GRAVITATIONAL/GRAVK,GRAV,GRAVKM,GCM2KM
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/MISHRA/IM05
      COMMON/SOLARPROPERTIES/RMSUN,RSUN,RSUNKM
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/MISCINT2/J651,J650,IGIVMX,IVONR,KC20
      COMMON/RADIALSCALINGFACTOR/R0R0
      COMMON/NUCLEARMATTERDENSITY/ENGNM0
      COMMON/IPLOT2/IOUTD4,IPOLXX,ICMI99

      PARAMETER (m651=160501)
      COMMON/SPL1/ST(M651),YOUT(M651),XSTROB(M651),YSTROB(M651)
      COMMON/CNSTA1/F251(M651),FOUT(M651),FST(M651),FST251(M651),
     +              FENG1(M651),X251(M651),Y251(M651),F251PM(M651),
     +              FDRIP(M651)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA3/AUX61(1,M651),AUX6(M651),DPDRAX(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)
      COMMON/BN01/AXBN1(M651),AXBN2(M651),PDENEX(M651),EINTHX(M651),
     +            DEINDP(M651),DEHDRH(M651),PDENT(M651)

      PARAMETER (NTHOV1=12)
      COMMON/TABSPH/THOV1(N27,NTHOV1)




      CHARACTER * 25 XIPEOS
      CHARACTER * 30 ARGINT
      CHARACTER * 30 EASS


      DATA RMC/0./
      DATA FMKM/1.E18/,FMEVFM/5.61012/

      NAMELIST/NAME1/RMNS,RMSUN,FOUT,FST251,EM,J650,J651



! ------------------------- STATEMENT FUNCTIONS -------------------------

!   OPPENHEIMER-VOLKOFF EQUATION

!   D [P] / D [R]
      FH(DEN,DPN)=-DGRAV*(DEN+DPN)*(DMN+VPI*A*DRN3*DPN)/
     +             (DRN2*(1.-2.*DGRAV*DMN/DRN))

!   NEWTONIAN LIMIT OF OV
!       FHNEWT(DEN)=-DGRAV*DEN*DMN / DRN2


!   METRIC FUNCTION

!   D [PHI] / D[R]
      FN(    DPN)=2.*DGRAV*(DMN+VPI*A*DRN3*DPN)/
     +             (DRN2*(1.-2.*DGRAV*DMN/DRN))

! -----------------------------------------------------------------------



      WRITE (6,'(//,''   _OVNROT_   ->  N O N - ROTATING STAR CAL'',
     +         ''CULATION;  MLF='',I3)') MLF
      WRITE (6,'(''                  Range of MLF: '',I2,'' <= MLF <='',
     +           I2)') MLFLL,MLFUL
      RX=BORD1*R0R0/EE18
      WRITE (6,'(''  NOTE: e_c='',E14.6,'' MeV/fm^3;  R_max='',E12.6,
     +           '' km'')') EC,RX

      IFAIL=10
      IF (IEP > NHU) THEN
      IFAIL=-10
      WRITE (6,'('' **  _OVNROT_  FIELD LENGTH OUT OF RANGE ** '')')
      WRITE (6,'(5X,'' IEP='',I4,5X,''NHU='',I4,///)') IEP,NHU
      RETURN
      END IF

      ISY1  =0
      LSTOP = 10
      IRF6  =KPIRF6
      IRF6M1=IRF6-1
      IDPREM=-10


!__ CONSTANTS
      IF (IKP1 /= KP1) THEN
      WRITE (6,'(//,'' ****   _OVNROT_   IKP1 # KP1  ****'',//,
     +              '' IKP1='',I4,3X,'' KP1='',I4)') IKP1,KP1
      IFAIL=10
      RETURN
      END IF

      AR0   =R0R0
      A     =EC*AR0**3/(RMSUN*ufe30)
      P16   =AR0/EE18
      DGRAV =(GRAV/ufe30)*(RMSUN*ufe30)/(UF1*AR0)
      DRMC  =RMC/(RMSUN*ufe30)
      RMETER=R0R0/EE15
      DPDRIP=PDRIP/EC


!__ INITIALIZE
      DO 8720 I=1,J651
      FDRIP(I) =0.
      AXBN1(I) =0.
      AXBN2(I) =0.
 8720 CONTINUE
      RLASPH   =0.
      RLACOR   =0.
      RLACRU   =0.




!__ ENERGY DENSITY AND PRESSURE (DIMENSIONLESS)
      DO 10    I=1,IEP
      AUX3(I)   =EDDY(I)/EC
      AUX4(I)   =PRESS(I)/EC
      AUX31(1,I)=AUX3(I)
      AUX41(1,I)=AUX4(I)
!   BARYON DENSITY (1/fm**3)
      AUX61(1,I)=PDEN(I)
   10 CONTINUE

      DEC    =EC/EC
!__ stop calculation if pressure  p < p_crust(Fe)
      DPCUTS =PCUTS/EC


!__ COMPUTE GRID FOR RADIAL INTEGRATION,  R=[0,1]

      WRITE (6,'(/,2x,''begin computation of radial grid:'')')
      IF (EDRIP == (-1.)) THEN
!     WIDE-SPACED GRID FOR NEUTRON/HYBRID/DWARF STAR CALCULATIONS
      BORD00=BORD0
      BORD11=BORD1
      ARGINT='OVNROT: RADIAL GRID - FST251  '
!                                  ******
      CALL INTV(BORD00,BORD11,J651,FST251,0,ARGINT)
      SSLT=(FST251(3)-FST251(2))*(R0R0/EE18)*1.E03
      SSGT=0.
      WRITE (6,'(2X,''step size for radial integration: '',e12.6,'' m'',
     +           /)') sslt

      ELSE IF (EDRIP /= (-1.)) THEN
!     NARROWLY-SPACED GRID FOR QUARK MATTER STAR CALCULATIONS, NUCLEAR
!     SOLID CRUST BELOW NUCLEAR DRIP DENSITY

      XDIV  =R0R0/(EE18*R5KM)
      BQS0  =0.
      BQS1  =1./XDIV
      ARGINT='OVNROT: RADIAL GRID - ST      '
!                                **             narrow grid,  [0,1/x]
      CALL INTV(BQS0,BQS1,IQ2020,ST,0,ARGINT)
      R6KM  =R5KM+1.
      XDIV  =R0R0/(EE18*R6KM)
      BQS0  =1./XDIV
      BQS1  =1.
      IH8080=J651-IQ2020


      ARGINT='OVNROT: RADIAL GRID - YOUT    '
!                                ****           wide grid,  [1/x,1]
      CALL INTV(BQS0,BQS1,IH8080,YOUT,0,ARGINT)
      SSLT=(ST(3)-ST(2))    *(R0R0/EE18)*1.E03
      SSGT=(YOUT(3)-YOUT(2))*(R0R0/EE18)*1.E03
      WRITE(6,'(2X,''step size for  r < '',e12.6,'' km: '',f7.4,'' m'')
     +         ') r5km,sslt
      WRITE(6,'(2X,''step size for  r > '',e12.6,'' km: '',e12.6,'' m'',
     +           /)') r5km,ssgt

!__ COPY
      DO 2511 I=1,J651
      IF (I <= IQ2020) THEN
      FST251(I)=ST(I)
      ELSE
      FST251(I)=YOUT(I-IQ2020)
      END IF
 2511 CONTINUE

      END IF
      FST251(1)=FST251(2)/2.



      IOPTIP=IPOLXX
      IF (IOPTIP < 0) THEN
      XIPEOS='LOGARITHMIC INTERPOLATION'
      ELSE IF (IOPTIP == 1) THEN
      XIPEOS='LINEAR INTERPOLATION     '
      ELSE IF (IOPTIP == 2) THEN
      XIPEOS='POLYNOMIAL INTERPOLATION '
      ELSE IF (IOPTIP == 3) THEN
      XIPEOS='SPLINE INTERPOLATION     '
      ELSE
      XIPEOS=' ***** NOT SPECIFIED ****'
      END IF


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     TOV EQUATION
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


!   STARTING VALUES
      EIP(1)  =DEC
      FENG1(1)=DEC
      F251(1) =DRMC
      F251(2) =VPI*A*DEC*FST251(1)**3/3.
      FST(1)  =0.
 

!  PRESSURE AT STAR'S ORIGIN         ***           RETURN: PIP=P_C
      CALL SELINT(AUX3,AUX41,1,1,IEP,PIP,DEC,1)
      FOUT(1)=PIP(IKP1)
      PN     =PIP(IKP1)
      DPDR00 =-VPI*A*DGRAV*(DEC+PN)*(DEC/3.+PN)*FST251(1)
      FOUT(2)=FOUT(1) + (FST251(2)-FST251(1))*DPDR00



!  BARYON DENSITY AT THE STAR'S CENTER
!                                    *****
      CALL SELINT(AUX3,AUX61,1,1,IEP,PDENS,DEC,1)
      PDENEX(1)    =PDENS(IKP1)
      THOV1(MLF,10)=PDENS(IKP1)

!  initialize N(1) and N(2)
      AXBN1(1)=  VPI * FST251(1)**3 * PDENEX(1) / 3.
      AXBN1(2)=  VPI * FST251(2)**3 * PDENEX(1) / 3.


! ******************************************************************
!                   START INTEGRATING TOV EQUATION
! ******************************************************************
      INT997=997
      NXY3  =3
      INTCT1=NXY3+INT997

      DO 1000 N=NXY3,J650

      IREM     =N-1
      DMN      =F251(N-1)
      DRN      =FST251(N-1)
      DRN2     =DRN*DRN
      DRN3     =DRN2*DRN
      DPN      =FOUT(N-1)
      PIP(IKP1)=DPN
      APRMEV   =DPN*EC


      IF (PIP(IKP1) <= AUX4(1)) THEN
      EIP(IKP1)=AUX3(1)


      ELSE



      IF (EDRIP == (-1.)) THEN

      IF (XIPEOS == 'SPLINE INTERPOLATION     ') THEN

      IF (APRMEV > PCUT) THEN

      EASS='OVNROT: CALL SELSPL AAA       '
!                                    ****
      CALL SELSPL(AUX4,AUX3,IEP,C,10,YYYY,DPN,10,EASS)
      EIP(IKP1)=YYYY
      EASS='OVNROT: CALL SELSPL BBB       '
!                                    ****
      CALL SELSPL(AUX4,PDEN,IEP,C,10,YYYY,DPN,10,EASS)
      YPD(IKP1)=YYYY
      ELSE

!                                    ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,-10)
!                                    ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,1)

! no access (logarithmic interpolation works better)
!      EASS='OVNROT: CALL SELSPL AAADDD    '
!                                      ****
!      CALL SELSPL8(AUX4,AUX3,IEP,C8,10,YYYY,DPN,10,EASS)
!      EIP(IKP1)=YYYY
!      EASS='OVNROT: CALL SELSPL BBBDDD    '
!                                      ****
!      CALL SELSPL8(AUX4,PDEN,IEP,C8,10,YYYY,DPN,10,EASS)
!      YPD(IKP1)=YYYY
! end of no access

      END IF

      ELSE
!                                    ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,IOPTIP)
!                                    ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,IOPTIP)
      END IF


      ELSE IF (EDRIP /= (-1.)) THEN


      LDA=(IDRIA+1) - IRF6
      LDB=(IDRIB-1) + IRF6


      IF (DPN <= AUX4(LDA) .OR. DPN >= AUX4(LDB)) THEN

      IF (DPN <= AUX4(LDA) .AND. IOUTD4 /= 0)
     +   PRINT*,' OV:  P < P_drip:  N=',N,' P(r)/e(c)=',DPN
      IF (DPN >= AUX4(LDB) .AND. IOUTD4 /= 0)
     +   PRINT*,' OV:  P > P_drip:  N=',N,' P(r)/e(c)=',DPN

      IF (XIPEOS == 'SPLINE INTERPOLATION     ') THEN

      IF (APRMEV > PCUT) THEN
      EASS='OVNROT: CALL SELSPL BBB       '
!                                    ****
      CALL SELSPL(AUX4,AUX3,IEP,C,10,YYYY,DPN,10,EASS)
      EIP(IKP1)=YYYY
      EASS='OVNROT: CALL SELSPL EEE       '
!                                    ****
      CALL SELSPL(AUX4,PDEN,IEP,C,10,YYYY,DPN,10,EASS)
      YPD(IKP1)=YYYY

      ELSE

!                                    ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,-10)
!                                    ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,-10)

      END IF


      ELSE
!                                    ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,IOPTIP)
!                                    ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,IOPTIP)
      END IF

      ELSE IF (DPN > AUX4(LDA) .AND. DPN < AUX4(LDB)) THEN

      IF (DPN <= DPDRIP) THEN
      IF(IOUTD4 /= 0) PRINT*,' OV:  P < P_drip: N=',N ,' P(r)/e(c)=',DPN
      IA1920=LDA
      IB1920=IDRIA
      LCLC  =LDA-1
      ELSE IF (DPN > DPDRIP) THEN
      IA1920=IDRIB
      IB1920=LDB
      LCLC  =IDRIB-1
      END IF

      ID20=IB1920-IA1920+1
      IF (ID20 /= IRF6) THEN
      WRITE (6,'(///,'' _OVNROT_  ID20 # IRF6 --> ERROR STOP'',/,
     +           '' ID20='',I4,'' IRF6='',I4,/,'' IB1920='',I4,/,
     +           '' IA1920='',I4,'' IDRIA='',I4,/,'' LCLC='',I4,/
     +           ,'' LDA='',I4,/, '' LDB='',I4,///)') ID20,IRF6,I
     +           B1920,IA1920,IDRIA,LCLC,LDA,LDB
      STOP
      END IF

      DO 1920 LC=IA1920,IB1920
      LS        =LC-LCLC
      A3(LS)    =AUX3(LC)
      A4(LS)    =AUX4(LC)
      P4(LS)    =PDEN(LC)
      A31(1,LS) =AUX31(1,LC)
      A61(1,LS) =AUX61(1,LC)
 1920 CONTINUE

      IF (XIPEOS == 'SPLINE INTERPOLATION     ') THEN

      IF (APRMEV > PCUT) THEN
      EASS='OVNROT: CALL SELSPL CCC       '
!                                 ****
      CALL SELSPL(A4,A3,IRF6,C,10,YYYY,DPN,10,EASS)
      EIP(IKP1)=YYYY
      EASS='OVNROT: CALL SELSPL HHH       '
!                                 ****
      CALL SELSPL(A4,P4,IRF6,C,10,YYYY,DPN,10,EASS)
      YPD(IKP1)=YYYY

      ELSE

!                                ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(A4,A31,1,1,IEP,EIP,DPN,-10)
!                                ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(A4,A61,1,1,IEP,YPD,DPN,-10)

      END IF



      ELSE
!                                 ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(A4,A31,1,1,IRF6,EIP,DPN,IOPTIP)
!                                 ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(A4,A61,1,1,IRF6,YPD,DPN,IOPTIP)
      END IF

             END IF

                    END IF

                           END IF



      DEN       =EIP(IKP1)
      FENG1(N-1)=DEN

!__ n(p) (1/fm^3)
      PDENEX(N-1)=YPD(IKP1)


      HRK       = FST251(N)-FST251(N-1)
      HRKH      = HRK/2.


      IF (IRUKOV == 0) THEN

!__ integration of Oppenheimer-Volkoff equations via Euler method
      DPDR      = FH(DEN,DPN)
      DPDRAX(N) = DPDR

!............................................
!>>>> P(N)
      PNP1  = PN + HRK * DPDR
!............................................


      ELSE

!__ preliminary computations to set up integration of Oppenheimer-Volkoff
!   equations via the fourth-order Runge-Kutta method below
      DPDR      = FH(DEN,DPN)
      DPDRAX(N) = DPDR
      XK1       = HRK * DPDR

      DRN       = FST251(N-1)+HRKH
      DRN2      = DRN*DRN
      DRN3      = DRN2*DRN
      DPK1      = DPN+XK1/2.
      DPDR2     = FH(DEN,DPK1)
      XK2       = HRK * DPDR2

      DPK2      = DPN+XK2/2.
      DPDR3     = FH(DEN,DPK2)
      XK3       = HRK * DPDR3

      DRN       = FST251(N-1)+HRK
      DRN2      = DRN*DRN
      DRN3      = DRN2*DRN
      DPK3      = DPN+XK3
      DPDR4     = FH(DEN,DPK3)
      XK4       = HRK * DPDR4

!................................................................
!>>>> P(N)
!     Compute advanced pressure via fourth order Runge-Kutta
      PNP1 = PN + (XK1 + XK4)/6. + (XK2 + XK3)/3.
!................................................................


      DRN       = FST251(N-1)
      DRN2      = DRN*DRN
      DRN3      = DRN2*DRN

      END IF



      PIP(IKP1)=PNP1
      FOUT(N)  =PNP1
      PN       =PNP1


      IF (N == INTCT1) THEN
      RKM1=FST251(N)*R0R0/EE18
      PMF3=PN*EC
      PRINT*,' N=',N,'  P=',PMF3,' MeV/fm^3','  r=',RKM1,' km'
      INTCT1=INTCT1+1000
                       END IF


      IF (LSTOP > 0) THEN
         IF ( (PN <= DPDRIP) .AND. (EQDRIP >= EDRIP) ) THEN
      PNP1     =DPDRIP * (1.0 - 1.0E-05)
      PIP(IKP1)=PNP1
      FOUT(N)  =PNP1
      PN       =PNP1
      LSTOP    =-10
            END IF
                     END IF


      IF      ( (FOUT(N) <= DPDRIP) .AND. (EC <= EQDRIP) ) THEN
!  IMPORTANT ONLY IN THE CASE OF WHITE DWARF CALCULATIONS
              IDPREM = 1
              GO TO 1991

      ELSE IF ( (FOUT(N) > DPDRIP) .AND. (EC <= EQDRIP) ) THEN

      WRITE (6,'(///,'' FATAL ERROR IN  _OVNROT_ '',/,'' P > P_drip'',
     +           '' or e_c > e_drip; are not compatible with each '',
     +           ''other for this'',/,'' calculation'',/,'' p='',e14.8,
     +           /,'' e='',e14.8,/,'' e_drip='',e14.8,/,'' p_drip='',
     +           e14.8)') FOUT(N)*EC,FENG1(N-1)*EC,EDRIP,DPDRIP*EC
      STOP '_OVNROT_  DRIP PRESSURE II'

         END IF



      IF (DPDRIP <= FOUT(N-1) .AND. DPDRIP > FOUT(N)) THEN
      IF (IM05 >= 0) GOTO 1771
      IDPREM=N-1
      FOUT(IDPREM)=DPDRIP
      POPMFM      =DPDRIP*EC
      PRINT*,' IDPREM=',IDPREM,' Note: P_IDPREM   > P_join=',PDRIP,
     +       ' MeV/fm^3'
      PRINT*,'                      P_IDPREM+1 < P_join'
      IF (EJGCM3 /= 0.)
     +PRINT*,' Note: e_join < e_drip (hadronic crust below n-drip)'
      IF (EJGCM3 == 0.)
     +PRINT*,' Note: e_join = e_drip (hadronic crust at n-drip)'
      WRITE (6,'(//,''  -----------------  P < P_join encountered'',
     +              ''  ----------------'',/)')
      ISY1=ISY1+1
      IF (ISY1 /= 1) THEN
      WRITE (6,'(///,'' FATAL ERROR IN  _OVNROT_ '',/,'' FAILED '',
     +           ''TO DETERMINE THE DRIP PRESSURE'',//,'' STOP '',
     +           ''CALCULATION;  ISY1='',I3,///)') ISY1
      STOP '_OVNROT_ DRIP PRESSURE'
      END IF
             END IF



 1771 CONTINUE

      IF (IDPREM < 0 .AND. PN <= 0. .AND. EDRIP /= (-1.)) THEN
      WRITE (6,'(///,'' FATAL ERROR IN  _OVNROT_ '',/,'' DRIP '',
     +           ''PRESSURE OUT OF RANGE'',/,'' P(HAT,DRIP)='',E12.6,
     +           /,'' EOS: P(1)='',E12.6,/,'' EOS: P(C)='',E12.6,
     +           /,'' IDPREM='',I4,/,'' EDRIP='',E12.6,///)')
     +           DPDRIP,AUX4(1),AUX4(IEP),IDPREM,EDRIP
      STOP
      END IF



 1991 CONTINUE

      PNLLIM=DPCUTS

      IF (PN <= PNLLIM) THEN
          THOV1(MLF,11)   =   PN*EC
          IF(PN > 0.) PN = - PNLLIM
          FOUT(N)         =   PN
          GO TO 2000
      END IF




      FST(N-1)=EIP(IKP1)*FST251(N-1)**2

      IDIFF1= 1
      IF (IDIFF1 < 0) THEN

      IF (IREM < 4) THEN
      WT=DEC*FST251(N-1)**3/3.
      ELSE IF (IREM >= 4) THEN

!..INTEGRATE OVER ENERGY DENSITY TO DETERMINE STAR'S MASS
!
!                        R             2
!                   WT = I [D R']  (R')   EPSILON(R'*R )
!                        0                            0

      EASS='OVNROT: CALL SELSPL DDD       '
!                                              **
      CALL IGRN(FST251,FST,IREM,X251,Y251,IREM,WT,C,NHU,10,1,AUX51,EASS)
!     RETURN:  WT
      END IF

!__ M A S S E S   (INTEGRATION OF MASS EQUATION):
!     GRAVITATIONAL MASS M(N-1)
      F251(N-1)  =VPI*A*WT

!     PROPER MASS: WILL NOT BE CALCULATED FOR  IDIFF1 < 0
      F251PM(N-1)=0.

!     CRUST MASS:  WILL NOT BE CALCULATED FOR  IDIFF1 < 0
      FDRIP(N-1) =0.

      ELSE

!__ MASSES  (VIA SOLVING THE DIFFERENTIAL EQUATION):
!     GRAVITATIONAL MASS
      F251(N)=F251(N-1)+HRK*VPI*A*DRN2*DEN

!     PROPER STAR MASS
      F251PM(N)=F251PM(N-1)+HRK*VPI*A*DRN2*DEN /
     +                    SQRT(1.-2.*DMN*DGRAV/DRN)

!     MASS OF THE STAR'S CRUST (DENSITY < DRIP DENSITY)
      IF (IDPREM < 0) THEN
      FDRIP(N)=0.
      ELSE
      FDRIP(N)=FDRIP(N-1)+HRK*VPI*A*DRN2*DEN
      END IF

      RAD =1.-2.*DGRAV*DMN/DRN
      SRAD=SQRT(RAD)


!__ STAR'S BARYON NUMBER (TOTAL, CORE, AND CRUST)
      IF (ICRUST == 0) THEN

!   TOTAL BARYON NUMBER
      AXBN1(N)=AXBN1(N-1) + HRK*DRN2*PDENEX(N-1) / SRAD
      ASPH    =AXBN1(N)
      AXBN2(N)=0.

      ELSE


      IF (IDPREM < 0) THEN

!   NUMBER OF BARYONS CONTAINED IN THE CORE (i.e., P > P_drip)
      AXBN1(N)=AXBN1(N-1) + HRK*DRN2*PDENEX(N-1) / SRAD
      ACORE   =AXBN1(N)
      AXBN2(N)=0.
      ELSE

!   NUMBER OF BARYONS IN THE HADRONIC CRUST (i.e., P < P_drip)
      AXBN1(N)=AXBN1(N-1)
      AXBN2(N)=AXBN2(N-1) + HRK*DRN2*PDENEX(N-1) / SRAD
      ACRUST  =AXBN2(N)

      END IF

              END IF


      END IF


 1000 CONTINUE
! ******************************************************************
!                        END OF INTEGRATION
! ******************************************************************

! NO ACCESS WRITE (6,'(5X,'' ** _OVNROT_ ERROR DETECTED  **'',////)')
! NO ACCESS WRITE (6,NAME1)
      IFAIL=-20
      RRNS =-1.
      RQOZ =-1.
      RETURN


!__ RADIUS AND MASS OF THE STAR

 2000 CONTINUE
      WRITE (6,'(''  surface of white dwarf/hyperon/hybrid/quark s'',
     +           ''tar encountered: irem='',I6,/)') irem

! DETERMINE VALUE OF THE STAR'S RADIUS (in km)
      K44  =5
      KHIGH=IREM+1
      KLOW =KHIGH-K44

      DO 581      I=KLOW,KHIGH
      IRUN         =I-KLOW+1
      AUX3(IRUN)   =FST251(I)
      AUX31(1,IRUN)=FOUT(I)
  581 CONTINUE

      IDIM=KHIGH-KLOW+1
      IOUT=40
      IYY =1


      WRITE (6,'(''  --> Computing stars radius'',//)')
!                                                        ****
      CALL ZERO22(AUX3,AUX31,IYY,IDIM,AUX4,AUX41,ST,IOUT,YOUT)
      FST251(IREM)=YOUT(IYY)

      RRNS  =FST251(IREM)*RMETER
      RQOZ  =F251(IREM)
      RMNS  =RQOZ*RMSUN*ufe30/5.61012E29



! COMPUTE AND COPY NUMBER OF BARYONS CONTAINED IN THE STAR'S CORE AND
! HADRONIC CRUST
      IF (ICRUST == 0) THEN
      RLASPH=3.*DLOG10(AR0) + ALOG10(VPI*ASPH)
      THOV1(MLF,7)=RLASPH
      ELSE

      IF(ACORE == 0.) THEN
      RLACOR=0.
      ELSE
      RLACOR=3.*DLOG10(AR0) + ALOG10(VPI*ACORE)

      END IF

      RLACRU=3.*DLOG10(AR0) + ALOG10(VPI*ACRUST)

      THOV1(MLF,8)=RLACOR
      THOV1(MLF,9)=RLACRU
      END IF


! AMOUNT OF THE STAR'S MASS WHICH IS BELOW THE DRIP DENSITY
      RMDRIS      =FDRIP(IREM)*RMSUN*EE55
      TOVDP(MLF,2)=RMDRIS
      TOVDP(MLF,3)=FDRIP(IREM)

      TOVDP(MLF,4) =F251PM(IREM)
      TOVDP(MLF,1) =0.
      IF (IDPREM > 0) TOVDP(MLF,1)=FST251(IDPREM)*P16
      IF (IDPREM == 1) TOVDP(MLF,1)=0.


      DO 481 I=0,NI
  481 RTH(I,1)=RRNS
      FOUT(IREM)=(FOUT(IREM-1)+FOUT(IREM-2))/10.



! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   INTEGRATE DIFFERENTIAL EQUATION OF THE METRIC FUNCTION  PHI(R)
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      FST(1)=0.

      DO 500 N=2,IREM
      DMN=F251(N-1)
      DRN=FST251(N-1)
      DRN2=DRN*DRN
      DRN3=DRN2*DRN
      DPN=FOUT(N-1)

!..................................
!>>>  PHI(R)
      HRK=FST251(N)-FST251(N-1)
      FST(N)=FST(N-1) + HRK*FN(DPN)
!..................................

  500 CONTINUE


!   RENORMALIZATION OF THE METRIC FUNCTION
!   (ADDITIV CONSTANT)
      EM=1.-2.*F251(IREM)*DGRAV/FST251(IREM)
      IF (EM < 0.) THEN
      WRITE (6,'(//,4X,'' **  EM < 0 **'',///)')
      WRITE (6,NAME1)
      RENK=0.
      GO TO 800
      ELSE
      RENK=ALOG(EM)-FST(IREM)
      END IF

!__ RENORMALIZE  PHI(R)
      DO 700 N=1,IREM
  700 FST(N)=FST(N)+RENK



! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!              MOMENT OF INERTIA OF THE STAR (I)
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      WI=0.
      IF (ICMI99 < 0) GO TO 6060

  800 DO 600 J=1,IREM
      DRN=FST251(J)
      DRN2=DRN*DRN
      DRN4=DRN2*DRN2
      DMN=F251(J)
      DEN=FENG1(J)
      DPN=FOUT(J)
      W=FST(J)/2.
      IF (W > EXMAX) THEN
      W=EXMAX
      ELSE IF (W < EXMIN) THEN
      W=EXMIN
      END IF
      EXPOF=EXP(W)
      RAD=1.-2.*DGRAV*DMN/DRN
      IF (RAD <= 0.) RAD=1.E-08
  600 FST(J)=DRN4*(DEN+DPN)/(EXPOF*SQRT(RAD))

!__ PERFORM INTEGRATION
!
!                        R             4
!                   TM = I [D R']  (R')  FUNCTION(R'*R )
!                        0                            0

      EASS='OVNROT: CALL SELSPL FFF       '
!                                              **
      CALL IGRN(FST251,FST,IREM,X251,Y251,IREM,TM,C,NHU,10,1,AUX51,EASS)
!     RETURN: TM
!__ COMPUTE  LOG_10 (I (g cm^2) )
      IF (TM <= 0.) THEN
      WI=1.
      ELSE
      WI= ALOG10(2.*VPI*EC/3.)+5.*ALOG10(R0R0)+ALOG10(TM)-ALOG10(FMEVFM)
     +   -52.
      END IF


 6060 CONTINUE

      IEPRES=0
      IF (IEPRES /= 0) THEN
         OPEN (UNIT=66,FILE='USER.E_OUT',STATUS='unknown')
         OPEN (UNIT=67,FILE='USER.P_OUT',STATUS='unknown')
         DO 3799 I=1,IREM,10
         EU=FENG1(I)*EC/ENGNM0
         PU=FOUT(I)*EC
         P16=AR0/EE18
         RU=FST251(I)*P16
         WRITE (66,610) RU,EU
         WRITE (67,610) RU,PU
 3799    CONTINUE
  610    FORMAT (2X,E14.7,4X,E14.7)

         CLOSE ( 66 )
                       CLOSE ( 67 )

      END IF


      RETURN

      END



!HT1000
!                                                             **** ****
      SUBROUTINE HT1000(EDDY,PRESS,IEP,EC,AUX3,AUX4,NHU,IFAIL,RSPH,RQOZ,
!                       **               ***
     +                  WI,C,AUX31,AUX41,RTH,NI,RQOZSP,THT,NIP1,IQM,MLF,
     +                  IREM,SUM124,OGC,XIPEOS,EINTNL,EINTH,PDEN,IATOTY,
!                                   ******       ******
     +                  ICMOI,SUMW1,RMTOTH,I5651,OGCOMS,EDRIP,EQDRIP,PDR
     +                  IP,IDRIA,IDRIB,RJI,NJ,NJP1,AX1,AX2,BORD0,BORD1,E
     +                  JGCM3,MLFLL,MLFUL,R5KM,IQ2020,SSLT,SSGT,PCUT,PCU
     +                  TS,C8,IRUKHT,POPMFM,R_90,e_90,R_00,e_00,cthr,anu
     +                  efip,k_qmc,k_mix,iwhiteDwarf)


!-----------------------------------------------------------------------
!
! PURPOSE:
!          SOLVE THE HARTLE-THORNE (HT) EQUATIONS OF A ROTATING RELATI-
!          VISTIC STAR (AP. J. 150 (1967) 1005; AP. J. 153 (1968) 807).
!
!          NOTE: OBTAIN THE PROPERTIES OF A RAPIDLY ROTATING COMPACT STAR
!                BY MEANS OF COMPUTING THE PERTURBATION SOLUTION ABOUT A
!                NON-ROTATING OPPENHEIMER-VOLKOFF STAR.
!
! INPUT:
!        EOS P(EPSILON)  (IN MeV/fm**3): EDDY,PRESS
!        CENTRAL ENERGY DENSITY EC=EPSILON_C  (IN MeV/fm**3)
!        (150 MeV/fm**3 < E_C < 5000 MeV/fm**3)
!        /TABOG1/RLOMG(.): OMEGA_LARGE(STAR)  (IN  1/SEC)
!
!
! CONSTANTS:
!        GRAVITATIONAL CONSTANT: GRAV (IN fm**2), GRAVK (IN fm/MeV)
!        SCALING FACTOR: R0R0 (IN fm)
!        MASS OF THE SUN RMSUN: (IN MeV)
!
!
! RETURN:
!         RSPH=RADIUS OF THE NON-ROTATING SPHERICAL STAR (IN METERS)
!         RTH=R_NS(0<=THETA<=#/2)
!         RQOZSP=M(SPHERICAL STAR)/M_SUN
!         RQOZ=M(NON-SPERICAL STAR)/M_SUN
!         WI=LOG(MOMENT OF INERTIA OF THE STAR (g cm^2))
!
!         IFAIL<0  <=>  ERROR RETURN
!-----------------------------------------------------------------------


      REAL * 8 DGRAV,A,AR0,AGD,P16,EE18,EE55,GCM2KM
      REAL * 8 RMSUN,RMNS,GRAV,GRAVK,GRAVKM,EEM18,EE14,EE34,EE44
      real * 8 ee03,ee13,ee15
      REAL * 8 C8(4,NHU)

      DIMENSION EDDY(NHU),PRESS(NHU),AUX3(NHU),AUX4(NHU),C(4,NHU)
      DIMENSION AUX31(1,NHU),AUX41(1,NHU),RTH(0:NI,4),THT(0:NI)
      DIMENSION EINTNL(NHU),EINTH(NHU),PDEN(NHU),R_90(0:NJ)
      DIMENSION RJI(0:NJ,0:NI),AX1(0:NJ),AX2(0:NJP1),e_90(0:NJ)
      DIMENSION R_00(0:NJ),e_00(0:NJ)
      REAL OGDIT(2)

      PARAMETER (KPIRF6=6)
      REAL A3(KPIRF6),A4(KPIRF6),P4(KPIRF6),E4(KPIRF6)
      REAL A31(1,KPIRF6),A61(1,KPIRF6),A71(1,KPIRF6)

      PARAMETER (KP1=1)
      DIMENSION PIP(KP1),EIP(KP1),YPD(KP1),YEI(KP1),DDO(KP1),PDENS(KP1),
     +          EINTC(KP1)

      COMMON/ExpoMaxMin/EXMAX,EXMIN
      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      COMMON/Exponent55/EE55
      common/ConversionFactor30/ufe30
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/GRAVITATIONAL/GRAVK,GRAV,GRAVKM,GCM2KM
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/MISHRA/IM05
      COMMON/SOLARPROPERTIES/RMSUN,RSUN,RSUNKM
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/MISCINT2/J651,J650,IGIVMX,IVONR,KC20
      COMMON/RADIALSCALINGFACTOR/R0R0
      COMMON/REDSEP/H0RCBE,H0RCBP,PHIRCP,PHIRCE,RDCBEP,RDCBEE,ODCBEE,
     +              H2RCBE,H2RCBP,V2RCBE,V2RCBP


      PARAMETER (m651=160501)
      COMMON/SPL1/ST(M651),YOUT(M651),XSTROB(M651),YSTROB(M651)
      COMMON/CNSTA1/F251(M651),FOUT(M651),FST(M651),FST251(M651),
     +              FENG1(M651),X251(M651),Y251(M651),F251PM(M651),
     +              FDRIP(M651)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA3/AUX61(1,M651),AUX6(M651),DPDRAX(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)
      COMMON/BN01/AXBN1(M651),AXBN2(M651),PDENEX(M651),EINTHX(M651),
     +            DEINDP(M651),DEHDRH(M651),PDENT(M651)
      COMMON/HT68A/DNUEDR(M651),RLJ(M651),DOGDR(M651),DJ2DR(M651)
      COMMON/HT68B/DEDP(M651),RNUED(M651),OGDR(M651),DLNJ(M651)
      COMMON/HT68C/RM0D(M651),P0SD(M651),DM0DR(M651),DP0SDR(M651)
      COMMON/HT68D/DMSDR(M651),DPHDRH(M651),P2SD(M651),RM2H(M651)
      COMMON/HT68E/H0OUT(M651),H0IN(M651),DH0IN(M651)
      COMMON/HT68F/DV2DR(M651),DH2DR(M651),V2RIN(M651),H2RIN(M651)
      COMMON/HT68G/Z0AR(M651),Z2AR(M651)
      COMMON/HT68H/DV2HOM(M651),DH2HOM(M651),V2HOM(M651),H2HOM(M651)
      COMMON/HT68I/OGAUX(M651),DOGAUX(M651),D2OGAX(M651)
      COMMON/MRDS2/DRM2H(M651),D2OGDR(M651),D2JDR2(M651),D2J22(M651)
      COMMON/MASSFRACTIONS/RMquark(m651),RMMixed(m651),RMNucl(m651)
      COMMON/CONVERGENCYLIMITS/ETASC6,ETASC7,ISCMX7
      COMMON/IPLOT1/IRADOT,IOUT20,IOTOMG,IOUT33,IOUTES,IOTMP0,IOP30,IOTH
     +              V2,IOUTPE,IOUTSE
      COMMON/IPLOT2/IOUTD4,IPOLXX,ICMI99
      common/Thresholds/thresh(15,2),rmtotm,rmtotq,nthrcb
      common/GivenInteriorDensity/einter

      PARAMETER (N27=50,N17=50,NP4=26,NPAPP=10)
      COMMON/TABHT1/THT1(N27,N17),THT2(N27,NP4),TAPPX(N27,NPAPP)
      common/tabquad/tquad(n27,npapp)
      COMMON/TABSTR/TSTROB(8,N27)
      common/tabmm/tabmm(8,n27)

      PARAMETER (NTHOV1=12)
      COMMON/TABSPH/THOV1(N27,NTHOV1)


      PARAMETER (NPOMG=50,NP10OG=13)
      COMMON/TABOG1/RLOMG(NPOMG),TFREQ(NPOMG,NP10OG)


      CHARACTER * 25 XIPEOS
      CHARACTER * 30 ARGINT
      CHARACTER * 30 EASS
      character * 30 cthr


      DATA RMC/0./
      DATA FMKM/1.E18/,FMEVFM/5.61012/

! 06/23/2015: parameter ALN used to compute O_K(i+1) in term of O_K(i)
      DATA ALN/0.15/      ! Use 0.15 for rotating white dwarf calculations
                          ! and 0.7 for rotating neutron star calculations
                          ! Appropritate value is assigned to ALN automatically
      

      NAMELIST/NAME1/RMNS,RMSUN,FOUT,FST251,EM,J650,J651



!....................... AUXILIARY FUNCTIONS ...........................

!                        OPPENHEIMER-VOLKOFF EQUATION

!                           D [P] / D [R]
      FH(DEN,DPN)=-DGRAV*(DEN+DPN)*(DMN+VPI*A*DRN3*DPN)/
     +             (DRN2*(1.-2.*DGRAV*DMN/DRN))

!                             METRIC FUNCTION

!                              D [PHI] / D[R]
      FN(    DPN)=2.*DGRAV*(DMN+VPI*A*DRN3*DPN)/
     +             (DRN2*(1.-2.*DGRAV*DMN/DRN))



!     AUX FUNCTIONS FOR INTEGRATING DIFFERENTIAL EQ FOR OMEGA_BAR

      FC1R(X,Y)=(4.+X*Y)/X
      FC2R(X,Y)= 4.*Y/X



!                                                            M
!     ASSOCIATED LEGENDRE POLYNOMIALS OF THE SECOND KIND:  Q  (Y)
!                                                           N
!     (N=2, M=1,2)
      Q21(Y,Y2)=-SQRT(Y2-1.) * ( 3. * Y * ALOG((Y+1.)/(Y-1.)) / 2.
     +          - (3.*Y2-2.) / (Y2-1.) )
      Q22(Y,Y2)=3. * (Y2-1.) * ALOG((Y+1.)/(Y-1.)) / 2.
     +          - Y * (3.*Y2 - 5.) / (Y2-1.)
      DQ22(Y,Y2,Z)=3.*Y*Z * ALOG((Y+1.)/(Y-1.)) - 3.*Z
     +            +2.*Y2*Z * (3.*Y2-5.) / ((Y2-1.)**2)
     +            -(9.*Y2-5.) * Z / (Y2-1.)
      DQ21(Y,Y2,Z)=Y*Z * ((3.*Y2-2.)/(Y2-1.) - 3.*Y * ALOG((Y+1.)/
     +                    (Y-1.)) /2.) / SQRT(Y2-1.)
     +             - SQRT(Y2-1.) * Z * ( 2.*Y/(Y2-1.)**2
     +             + 3. * ALOG((Y+1)/(Y-1.)) / 2. - 3.*Y/(Y2-1.) )



! ....................... END OF AUX FUNCTIONS .........................

      if (iwhiteDwarf == 0) ALN=0.7
      
      ltc=0
      ISY1=0
      LSTOP= 10
      IFAIL=10
      IF (IEP > NHU) THEN
      IFAIL=-10
      WRITE (6,'('' HT1000: FIELD LENGTH OUT OF RANGE ** '')')
      WRITE (6,'(5X,'' IEP='',I4,5X,''NHU='',I4,///)') IEP,NHU
      RETURN
      END IF

! 08/09/2018
! Initialize nthr (used to determine amount of quark matter present at a given density)
      if(ioutse == 10) then
      nthr=1
      cthr='G_240^B180 ...................'
      else if(ioutse == 20) then
      nthr=2
      cthr='G_300^B180 ...................'
      else if(ioutse == 30) then
      nthr=3
      cthr='G_300^B180_m7_a0 .............'
      else if(ioutse == 40) then
      nthr=4
      cthr='G_350^B180_m73 ...............'
      else if(ioutse == 50) then
      nthr=5
      cthr='G_290^B180_m7 ................'
      else if(ioutse == 60) then
      nthr=6
      cthr='Alford CFL1 ..................'
      else if(ioutse == 70) then
      nthr=7
      cthr='epNL_GM1_Gv009................'
      else if(ioutse == 80) then
      nthr=8
      cthr='DD2_npemu_ms119_gv00_hybrid...'
      else          
      write (6,'(//, ''  _ht1000_  Warning: ioutse=0 -> particle thre'',
     +               ''holds not provided for this nuclear EoS'', //)')
      cthr='ioutse (input) = 0            '
      end if

      nthrcb=nthr

      WRITE (6,'(''  HT1000:  ->  BEGIN ROTATING STAR CALCULATION'',
     +           '';  MLF='',I3)') MLF

      WRITE (6,'(''                  Range of MLF: '',I2,'' <= MLF <='',
     +           I2)') MLFLL,MLFUL

      RX=BORD1*R0R0/EE18
      WRITE (6,'(''  NOTE: e_c='',E14.6,'' MeV/fm^3; R_(max;ht)='',E12.6,
     +           '' km'')') EC,RX

      IRF6  =KPIRF6
      IRF6M1=IRF6-1
      IDPREM=-10

!__ CONSTANTS AND PARAMETERS
      IF (IKP1 /= KP1 .AND. IRF6 /= KPIRF6) THEN
      WRITE (6,'(//,'' HT1000:  IKP1 /= KP1  ****'',//,
     +              '' IKP1='',I4,3X,'' KP1='',I4)') IKP1,KP1
      IFAIL=10
      RETURN
      END IF
      AR0   =R0R0
      A     =EC*AR0**3/(RMSUN*ufe30)
      P16   =AR0/EE18
      DGRAV =(GRAV/ufe30)*(RMSUN*ufe30)/(UF1*AR0)
      DRMC  =RMC/(RMSUN*ufe30)
      AGD   =A*DGRAV
      RMETER=R0R0/EE15
      DPDRIP=PDRIP/EC
      k_mix=0
      k_qmc=0
      rm_mix=0.
      rm_qmc=0.
      rm_nlc=0.
      rmtotq=0.
      rmtotm=0.
      rmtoth=0.
      rmtotc=0.
      delmh =0.
      delmix=0.
      delqmc=0.

!__ INITIALIZE
      DO 8720 I=1,J651
      FDRIP(I)  =0.
      PDENT(I)  =0.
      RMQuark(i)=0.
      RMMixed(i)=0.
      RMNucl(i) =0.
 8720 CONTINUE

      H0RCBE=0.
      H0RCBP=0.
      PHIRCP=0.
      PHIRCE=0.
      RDCBEP=0.
      RDCBEE=0.
      ODCBEE=0.
      H2RCBE=0.
      H2RCBP=0.
      V2RCBE=0.
      V2RCBP=0.


!__ COPY
      DO 10 I=1,IEP

!   ENERGY DENSITY AND PRESSURE (DIMENSIONLESS)
      AUX3(I)   =EDDY(I)/EC
      AUX4(I)   =PRESS(I)/EC
      AUX31(1,I)=AUX3(I)
      AUX41(1,I)=AUX4(I)

!   INTERNAL ENERGY (DIMENSIONLESS)
      EINTH(I)  =EINTNL(I)/EC
      AUX71(1,I)=EINTH(I)

!   BARYON DENSITY (1/fm**3)
      AUX61(1,I)=PDEN(I)
   10 CONTINUE

      DEC    =EC/EC
!     stop calculation if pressure  p < p_crust(Fe)
      DPCUTS =PCUTS/EC


!__ COMPUTE GRID FOR RADIAL INTEGRATION,  R=[0,1]

      WRITE (6,'(/,2x,''begin computation of radial grid:'')')
      IF (EDRIP == (-1.)) THEN
!     WIDE-SPACED GRID FOR NEUTRON/HYBRID/DWARF STAR CALCULATIONS
      BORD00=BORD0
      BORD11=BORD1
      ARGINT='OVNROT: RADIAL GRID - FST251  '
!                                  ******
      CALL INTV(BORD00,BORD11,J651,FST251,0,ARGINT)
      SSLT=(FST251(3)-FST251(2))*(R0R0/EE18)*1.E03
      SSGT=0.
      WRITE (6,'(2X,''step size for radial integration: '',e12.6,'' m'',
     +           /)') sslt

      ELSE IF (EDRIP /= (-1.)) THEN
!     NARROWLY-SPACED GRID FOR QUARK MATTER STAR CALCULATIONS, NUCLEAR
!     SOLID CRUST BELOW NUCLEAR DRIP DENSITY

      XDIV=R0R0/(EE18*R5KM)
      BQS0  =0.
      BQS1  =1./XDIV
      ARGINT='OVNROT: RADIAL GRID - ST      '
!                                **             narrow grid,  [0,1/x]
      CALL INTV(BQS0,BQS1,IQ2020,ST,0,ARGINT)
      R6KM  =R5KM+1.
      XDIV  =R0R0/(EE18*R6KM)
      BQS0  =1./XDIV
      BQS1  =1.
      IH8080=J651-IQ2020


      ARGINT='OVNROT: RADIAL GRID - YOUT    '
!                                ****           wide grid,  [1/x,1]
      CALL INTV(BQS0,BQS1,IH8080,YOUT,0,ARGINT)
      SSLT=(ST(3)-ST(2))    *(R0R0/EE18)*1.E03
      SSGT=(YOUT(3)-YOUT(2))*(R0R0/EE18)*1.E03
      WRITE(6,'(2X,''step size for  r < '',e12.6,'' km: '',f7.4,'' m'')
     +         ') r5km,sslt
      WRITE(6,'(2X,''step size for  r > '',e12.6,'' km: '',e12.6,'' m'',
     +           /)') r5km,ssgt

!__ COPY
      DO 2511 I=1,J651
      IF (I <= IQ2020) THEN
      FST251(I)=ST(I)
      ELSE
      FST251(I)=YOUT(I-IQ2020)
      END IF
 2511 CONTINUE

      END IF
      FST251(1)=FST251(2)/2.


!__ COMPUTE THETA GRID POINTS
      ABORD=0.
      BBORD=PID2
      ARGINT='HT1000: THETA GRID  [0,#/2]   '
      CALL INTV(ABORD,BBORD,NIP1,X251,0,ARGINT)
      DO 530 I=0,NI
  530 THT(I)  =X251(I+1)


      IOPTIP=IPOLXX
      IF (IOPTIP < 0) THEN
      XIPEOS='LOGARITHMIC INTERPOLATION'
      ELSE IF (IOPTIP == 1) THEN
      XIPEOS='LINEAR INTERPOLATION     '
      ELSE IF (IOPTIP == 2) THEN
      XIPEOS='POLYNOMIAL INTERPOLATION '
      ELSE IF (IOPTIP == 3) THEN
      XIPEOS='SPLINE INTERPOLATION     '
      ELSE
      XIPEOS=' ***** NOT SPECIFIED ****'
      END IF



! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!           SIMULTANEOUS INTEGRATION OF OV EQUATONS
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


!   STARTING VALUES
      EIP(1)   =DEC
      FENG1(1) =DEC
      F251(1)  =DRMC
      F251(2)  =VPI*A*DEC*FST251(1)**3/3.
      F251PM(1)=F251(1)
      F251PM(2)=F251(2)
      RNUED(1) =0.
      FST(1)   =0.
      DMSDR(1) =VPI*A*FST251(1)**2*FENG1(1)
      DNUEDR(1)=2.*VPI*A*DGRAV * (DEC/3.+FOUT(1)) * FST251(1)

!   PRESSURE AT THE STAR'S CENTER
!                                    ***
      CALL SELINT(AUX3,AUX41,1,1,IEP,PIP,DEC,1)
!     RETURN: PIP=P_C

      FOUT(1)  =PIP(IKP1)
      PN       =PIP(IKP1)
      DPHDRH(1)=-DGRAV*VPI*A*(DEC+PN)*(DEC/3.+PN)*FST251(1)
      FOUT(2)  =FOUT(1) + (FST251(2)-FST251(1))*DPHDRH(1)

      tstrob(2,MLF) = fout(1)*ec*uf4/ee34

!   BARYON DENSITY AT THE STAR'S CENTER
!                                    *****
      CALL SELINT(AUX3,AUX61,1,1,IEP,PDENS,DEC,1)
      PDENEX(1)   =PDENS(IKP1)
      THT1(MLF,25)=PDENS(IKP1)

!   INTERNAL ENERGY AT THE STAR'S CENTER
!                                    *****
      CALL SELINT(AUX3,AUX71,1,1,IEP,EINTC,DEC,1)
      EINTHX(1)   =EINTC(IKP1)
      THT1(MLF,26)=EINTC(IKP1)*EC




! ******************************************************************
!                 START INTEGRATING OV EQS
! ******************************************************************
      INT997=997
      NXY3  =3
      INTCT1=NXY3+INT997

      DO 1000 N=NXY3,J650

      IREM=N-1
      DMN =F251(N-1)
      DRN =FST251(N-1)
      DRN2=DRN*DRN
      DRN3=DRN2*DRN
      DPN =FOUT(N-1)
      PIP(IKP1)=DPN
      APRMEV   =DPN*EC


      IF (PIP(IKP1) <= AUX4(1)) THEN
      EIP(IKP1)=AUX3(1)

      ELSE

      IF (EDRIP == (-1.)) THEN

      IF (XIPEOS == 'SPLINE INTERPOLATION     ') THEN

      IF (APRMEV > PCUT) THEN

      EASS='HT1000: CALL SELSPL AAA       '
!                                    ****
      CALL SELSPL(AUX4,AUX3,IEP,C,10,YYYY,DPN,10,EASS)
      EIP(IKP1)=YYYY
      EASS='HT1000: CALL SELSPL BBB       '
!                                    ****
      CALL SELSPL(AUX4,PDEN,IEP,C,10,YYYY,DPN,10,EASS)
      YPD(IKP1)=YYYY
      EASS='HT1000: CALL SELSPL BBBB      '
!                                     ****
      CALL SELSPL(AUX4,EINTH,IEP,C,10,YYYY,DPN,10,EASS)
      YEI(IKP1)=YYYY

      ELSE

!                                    ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,-10)
!                                    ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,-10)
!                                    ***            RETURN: YEI=E(int)_P(N-1)
      CALL SELINT(AUX4,AUX71,1,1,IEP,YEI,DPN,  1)

! no access (logarithmic interpolation works better)
!      EASS='HT1000: CALL SELSPL AAADDD    '
!                                      ****
!      CALL SELSPL8(AUX4,AUX3,IEP,C8,10,YYYY,DPN,10,EASS)
!      EIP(IKP1)=YYYY
!      EASS='HT1000: CALL SELSPL BBBDDD    '
!                                      ****
!      CALL SELSPL8(AUX4,PDEN,IEP,C8,10,YYYY,DPN,10,EASS)
!      YPD(IKP1)=YYYY
!      EASS='HT1000: CALL SELSPL BBBBDDD   '
!                                       ****
!      CALL SELSPL8(AUX4,EINTH,IEP,C8,10,YYYY,DPN,10,EASS)
!      YEI(IKP1)=YYYY
! end of no access

      END IF

      ELSE


!                                    ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,IOPTIP)
!                                    ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,IOPTIP)
!                                    ***            RETURN: YEI=E(int)_P(N-1)
      CALL SELINT(AUX4,AUX71,1,1,IEP,YEI,DPN,IOPTIP)

      END IF



      ELSE IF (EDRIP /= (-1.)) THEN


      LDA=(IDRIA+1) - IRF6
      LDB=(IDRIB-1) + IRF6


      IF (DPN <= AUX4(LDA) .OR. DPN >= AUX4(LDB)) THEN

      IF (DPN <= AUX4(LDA) .AND. IOUTD4 /= 0)
     +   PRINT*,'HT:  P < P_drip:  N=',N,' P(r)/e(c)=',DPN

      IF (XIPEOS == 'SPLINE INTERPOLATION     ') THEN

      IF (APRMEV > PCUT) THEN
      EASS='HT1000: CALL SELSPL DDD       '
!                                    ****
      CALL SELSPL(AUX4,AUX3,IEP,C,10,YYYY,DPN,10,EASS)
      EIP(IKP1)=YYYY
      EASS='HT1000: CALL SELSPL EEE       '
!                                    ****
      CALL SELSPL(AUX4,PDEN,IEP,C,10,YYYY,DPN,10,EASS)
      YPD(IKP1)=YYYY
      EASS='HT1000: CALL SELSPL FFF       '
!                                     ****
      CALL SELSPL(AUX4,EINTH,IEP,C,10,YYYY,DPN,10,EASS)
      YEI(IKP1)=YYYY


      ELSE

!                                    ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,-10)
!                                    ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,-10)
!                                    ***            RETURN: YEI=E(int)_P(N-1)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL SELINT(AUX4,AUX71,1,1,IEP,YEI,DPN,  1)

      END IF


      ELSE
!                                    ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,IOPTIP)
!                                    ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,IOPTIP)
!                                    ***            RETURN: YEI=E(int)_P(N-1)
      CALL SELINT(AUX4,AUX71,1,1,IEP,YEI,DPN,IOPTIP)
      END IF

      ELSE IF (DPN > AUX4(LDA) .AND. DPN < AUX4(LDB)) THEN

      IF (DPN <= DPDRIP) THEN
      IF(IOUTD4 /= 0) PRINT*,'HT:  P < P_drip:  N=',N,' P(r)/e(c)=',DPN
      IA1920=LDA
      IB1920=IDRIA
      LCLC  =LDA-1
      ELSE IF (DPN > DPDRIP) THEN
      IA1920=IDRIB
      IB1920=LDB
      LCLC  =IDRIB-1
      END IF
      ID20=IB1920-IA1920+1
      IF (ID20 /= IRF6) THEN
      WRITE (6,'(///,'' HT1000:  ID20 # IRF6 --> ERROR STOP'',/,
     +           '' ID20='',I4,'' IRF6='',I4,/,'' IB1920='',I4,/,
     +           '' IA1920='',I4,'' IDRIA='',I4,/,'' LCLC='',I4,/
     +           ,'' LDA='',I4,/, '' LDB='',I4,///)') ID20,IRF6,I
     +           B1920,IA1920,IDRIA,LCLC,LDA,LDB
      STOP
      END IF

      DO 1920 LC=IA1920,IB1920
              LS=LC-LCLC
      A3(LS)   =AUX3(LC)
      A4(LS)   =AUX4(LC)
      P4(LS)   =PDEN(LC)
      E4(LS)   =EINTH(LC)
      A31(1,LS)=AUX31(1,LC)
      A61(1,LS)=AUX61(1,LC)
      A71(1,LS)=AUX71(1,LC)
 1920 CONTINUE


      IF (XIPEOS == 'SPLINE INTERPOLATION     ') THEN

      IF (APRMEV > PCUT) THEN
      EASS='HT1000: CALL SELSPL GGG       '
!                                 ****
      CALL SELSPL(A4,A3,IRF6,C,10,YYYY,DPN,10,EASS)
      EIP(IKP1)=YYYY
      EASS='HT1000: CALL SELSPL HHH       '
!                                 ****
      CALL SELSPL(A4,P4,IRF6,C,10,YYYY,DPN,10,EASS)
      YPD(IKP1)=YYYY
      EASS='HT1000: CALL SELSPL III       '
!                                 ****
      CALL SELSPL(A4,E4,IRF6,C,10,YYYY,DPN,10,EASS)
      YEI(IKP1)=YYYY

      ELSE

!                                    ***          RETURN: EIP=E_P(N-1)
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,-10)
!                                    ***          RETURN: YPD=n_P(N-1)
      CALL SELINT(AUX4,AUX61,1,1,IEP,YPD,DPN,-10)
!                                 ***             RETURN: YEI=E(int)_P(N-1)
      CALL SELINT(A4,A71,1,1,IRF6,YEI,DPN,-10)

      END IF


      ELSE
!                                 ***            RETURN: EIP=E_P(N-1)
      CALL SELINT(A4,A31,1,1,IRF6,EIP,DPN,IOPTIP)
!                                 ***            RETURN: YPD=n_P(N-1)
      CALL SELINT(A4,A61,1,1,IRF6,YPD,DPN,IOPTIP)
!                                 ***            RETURN: YEI=E(int)_P(N-1)
      CALL SELINT(A4,A71,1,1,IRF6,YEI,DPN,IOPTIP)
      END IF

             END IF

                    END IF

                           END IF

      DEN       =EIP(IKP1)
      FENG1(N-1)=DEN
 
! 06/29/05: begin 
      if (ltc == 0) then
      idprem = 0            ! >< integrate from [0,surface] in IDEFST
      cdens=einter/uf6/ec   ! [einter]=g/cm^3, [cdens]=1
      if ( (feng1(n-2) > cdens) .and. (feng1(n-1) < cdens) ) then 
         idprem = n
         print *, feng1(n-2)*ec, feng1(n-1)*ec, n, cdens*ec, einter
         print *, ' Moment of inertia of crust (<2.4e14 g/cm^3) is ',
     +            'computed! ', 'idprem=', idprem
         ltc=ltc+1
      end if
             end if
! 06/29/05: end


!__ n(P) (1/fm^3)
      PDENEX(N-1)=YPD(IKP1)

!__ E(int;P)  [1]
      EINTHX(N-1)=YEI(IKP1)

!__ D [P] / D [R]
      DPDR       =FH(DEN,DPN)
      DPDRAX(N)  =DPDR
      DPHDRH(N-1)=DPDR

!__ D [PHI] / D [R]
      DNUEDR(N-1)=FN(FOUT(N-1))

!__ D [M] / D [R]
      DMSDR(N-1)=VPI*A * FST251(N-1)**2 * FENG1(N-1)



      HRK       = FST251(N)-FST251(N-1)
      HRKH      = HRK/2.


      IF (IRUKHT == 0) THEN

!__ integration of Oppenheimer-Volkoff equations via Euler method

!............................................
!>>>> P(N)
      PNP1  = PN + HRK * DPDR
!............................................


      ELSE

!__ preliminary computations for integrating the Oppenheimer-Volkoff
!   equations via the fourth-order Runge-Kutta method below
      XK1       = HRK * DPDR

      DRN       = FST251(N-1)+HRKH
      DRN2      = DRN*DRN
      DRN3      = DRN2*DRN
      DPK1      = DPN+XK1/2.
      DPDR2     = FH(DEN,DPK1)
      XK2       = HRK * DPDR2

      DPK2      = DPN+XK2/2.
      DPDR3     = FH(DEN,DPK2)
      XK3       = HRK * DPDR3

      DRN       = FST251(N-1)+HRK
      DRN2      = DRN*DRN
      DRN3      = DRN2*DRN
      DPK3      = DPN+XK3
      DPDR4     = FH(DEN,DPK3)
      XK4       = HRK * DPDR4

!................................................................
!>>>> P(N)
!     Compute advanced pressure via fourth order Runge-Kutta
      PNP1 = PN + (XK1 + XK4)/6. + (XK2 + XK3)/3.
!................................................................


      DRN       = FST251(N-1)
      DRN2      = DRN*DRN
      DRN3      = DRN2*DRN

      END IF



      PIP(IKP1)=PNP1
      FOUT(N)  =PNP1
      PN       =PNP1


      IF (N == INTCT1) THEN
      RKM1=FST251(N)*R0R0/EE18
      PMF3=PN*EC
      PRINT*,' N=',N,'  P=',PMF3,' MeV/fm^3','  r=',RKM1,' km'
      INTCT1=INTCT1+1000
                       END IF


      IF (LSTOP > 0) THEN
         IF ( (PN <= DPDRIP) .AND. (EQDRIP >= EDRIP) ) THEN
      PNP1=DPDRIP * (1. - 1.E-05)
      PIP(IKP1)=PNP1
      FOUT(N)=PNP1
      PN=PNP1
      LSTOP=-10
             END IF
                      END IF


      IF      ( (FOUT(N) <= DPDRIP) .AND. (EC <= EQDRIP) ) THEN
! IMPORTANT ONLY FOR WHITE DWARF CALCULATIONS
              IDPREM = 1
              GO TO 1991


      ELSE IF ( (FOUT(N) > DPDRIP) .AND. (EC <= EQDRIP) ) THEN

      WRITE (6,'(///,'' FATAL ERROR IN  HT1000: '',/,'' P>P_drip'',
     +           '' or e_c>e_drip; are not compatible with each '',
     +           ''other for this calculation'',/,'' p='',e14.8,/,
     +           '' e='',e14.8,/,'' e_drip='',e14.8,/,'' p_drip='',
     +           e14.8)') FOUT(N)*EC,FENG1(N)*EC,EDRIP,EQDRIP
      STOP ' HT1000: DRIP PRESSURE II'

         END IF


      IF (DPDRIP <= FOUT(N-1) .AND. DPDRIP > FOUT(N)) THEN
      IF (IM05 >= 0) GOTO 1771
      IDPREM=N-1
      FOUT(IDPREM)=DPDRIP
      POPMFM      =DPDRIP*EC
      PRINT*,' IDPREM=',IDPREM,' Note: P_IDPREM   > P_join=',PDRIP,
     +       ' MeV/fm^3'
      PRINT*,'                      P_IDPREM+1 < P_join'
      IF (EJGCM3 /= 0.)
     +PRINT*,' Note: e_join < e_drip (hadronic crust below n-drip)'
      IF (EJGCM3 == 0.)
     +PRINT*,' Note: e_join = e_drip (hadronic crust at n-drip)'

      WRITE (6,'(//,''  -----------------  P < P_join encountered'',
     +              ''  ----------------'',/)')
      ISY1=ISY1+1
      IF (ISY1 /= 1) THEN
      WRITE (6,'(///,'' FATAL ERROR IN  HT1000:'',/,'' FAILED '',
     +           ''TO DETERMINE THE DRIP PRESSURE'',//,'' STOP '',
     +           ''CALCULATION;  ISY1='',I3,///)') ISY1
      STOP ' HT1000: DRIP PRESSURE'
      END IF
             END IF


 1771 CONTINUE
      IF (IDPREM < 0 .AND. PN <= 0. .AND. EDRIP /= (-1.)) THEN
      WRITE (6,'(///,'' FATAL ERROR IN  HT1000: '',/,'' DRIP '',
     +           ''PRESSURE OUT OF RANGE'',/,'' P(HAT,DRIP)='',E12.6,
     +           /,'' EOS: P(1)='',E12.6,/,'' EOS: P(C)='',E12.6,
     +           /,'' IDPREM='',I4,/,'' EDRIP='',E12.6,///)')
     +           DPDRIP,AUX4(1),AUX4(IEP),IDPREM,EDRIP
      STOP
      END IF


 1991 CONTINUE

      PNLLIM=DPCUTS

      IF (PN <= PNLLIM) THEN
          THOV1(MLF,12)   = PN*EC
          IF(PN > 0.) PN = - PNLLIM
          FOUT(N)         = PN
          GO TO 2000
      END IF


      FST(N-1)=EIP(IKP1)*FST251(N-1)**2

      IDIFF1= 1
      IF (IDIFF1 < 0) THEN

      IF (IREM < 4) THEN
      WT=DEC*FST251(N-1)**3/3.
      ELSE IF (IREM >= 4) THEN

!..INTEGRATE ENERGY DENSITY TO GET STAR'S MASS
!
!                        R             2
!                   WT = I [d R']  (R')   EPSILON(R'*R )
!                        0                            0

      EASS='HT1000: CALL SELSPL JJJ       '
!                                              **
      CALL IGRN(FST251,FST,IREM,X251,Y251,IREM,WT,C,NHU,10,1,AUX51,EASS)
!     RETURN:  WT
      END IF

!__ M A S S E S   (INTEGRATE MASS EQUATION):
!     GRAVITATIONAL MASS M(N-1)
      F251(N-1)  =VPI*A*WT

!     PROPER MASS: WILL NOT BE CALCULATED FOR  IDIFF1 < 0
      F251PM(N-1)=0.

!     CRUST MASS:  WILL NOT BE CALCULATED FOR  IDIFF1 < 0
      FDRIP(N-1) =0.

      ELSE       ! idiff1 > 0

!__ MASSES  (VIA SOLVING THE DIFFERENTIAL EQUATION):
!     GRAVITATIONAL MASS
      F251(N)=F251(N-1)+HRK*VPI*A*DRN2*DEN

!     PROPER STAR MASS
      F251PM(N)=F251PM(N-1)+HRK*VPI*A*DRN2*DEN /
     +                    SQRT(1.-2.*DMN*DGRAV/DRN)

!     MASS OF THE STAR'S CRUST (DENSITY < DRIP DENSITY)
      IF (IDPREM < 0) THEN
      FDRIP(N)=0.
      ELSE
      FDRIP(N)=FDRIP(N-1)+HRK*VPI*A*DRN2*DEN
      END IF

      
! For neutron stars with quark-hybrid composition only

      if (ioutse /= 0) then 

! Compute mass in pure quark matter core
      if( (feng1(n-1)*ec) >= thresh(nthr,2) ) then
      rm_qmc=F251(N-1)
      RMquark(N-1)=rm_qmc
      k_qmc =n-1

! Compute mass in mixed phase 
      else if( ((feng1(n-1)*ec) >= thresh(nthr,1)) .and. 
     +         ((feng1(n-1)*ec) < thresh(nthr,2)) ) then 
      rm_mix=F251(N-1)-rm_qmc
      RMMixed(N-1)=rm_mix
      k_mix =n-1

! Compute mass in nuclear liquid and crust 
      else if( (feng1(n-1)*ec) < thresh(nthr,1) ) then
      rm_nlc=F251(N-1)-rm_qmc-rm_mix
      RMNucl(N-1)=rm_nlc
      end if

      end if   ! end ioutse /= 0 case

      
      END IF 


 1000 CONTINUE


! ******************************************************************
!                         END OF INTEGRATION
! ******************************************************************


! NO ACCESS WRITE (6,'(5X,'' ** HT1000_ ERROR DETECTED  **'',////)')
! NO ACCESS WRITE (6,NAME1)
      IFAIL =-20
      RSPH  =-1.
      RQOZSP=-1.
      RETURN



 2000 CONTINUE
      write (6,'(''  M_qmc='',f7.4,'' M_mix='',f7.4,'' k_qmc='',i6,
     +           '' k_mix='',i6,'' M_nlc='',f7.4)') rm_qmc,rm_mix,
     +       k_qmc,k_mix,rm_nlc
      WRITE (6,'(''  surface of white dwarf/hyperon/hybrid/quark s'',
     +           ''tar encountered: irem='',I6,/)') irem


!__ DETERMINE RADIUS (in km) AND MASS OF THE SPHERICAL PART OF THE
!   ROTATING NEUTRON/HYBRID/QUARK STAR, OR DWARF

      K44  =5
      KHIGH=IREM+1
      KLOW =KHIGH-K44

      DO 581      I=KLOW,KHIGH
      IRUN         =I-KLOW+1
      AUX5(IRUN)   =FST251(I)
      AUX51(1,IRUN)=FOUT(I)
  581 CONTINUE

      IDIM=KHIGH-KLOW+1
      IOUT=40
      IYY =1
      WRITE (6,'(''  --> Computing stars radius'',//)')
!                                                        ****
      CALL ZERO22(AUX5,AUX51,IYY,IDIM,AUX4,AUX41,ST,IOUT,YOUT)
      FST251(IREM)=YOUT(IYY)

      RSPH  =FST251(IREM)*RMETER
      RQOZSP=F251(IREM)
      RMNS  =RQOZSP*RMSUN*ufe30/5.61012E29
      
!     AMOUNT OF THE STAR'S MASS WHICH IS BELOW THE DRIP DENSITY
      RMDRIS      =FDRIP(IREM)*RMSUN*EE55
      THT1(MLF,5) =RSPH/EE03
      THT1(MLF,10)=RQOZSP
      THT2(MLF,4) =F251PM(IREM)
      THT2(MLF,15)=RMDRIS
      THT2(MLF,18)=FDRIP(IREM)
      THT2(MLF,21)=0.
      IF (IDPREM > 0) THT2(MLF,21)=FST251(IDPREM)*P16
      IF (IDPREM == 1) THT2(MLF,21)=0.

      DO 481 I=0,NI
  481 RTH(I,1)=RSPH
      FOUT(IREM)=(FOUT(IREM-1)+FOUT(IREM-2))/10.


!>>>> PRINT CALCULATED FUNCTIONS
      IF (IOUT20 /= 0) THEN
      WRITE (6,'(//)')
      DO 67 I=1,IREM
   67 PRINT*,' I=',I,'  P=',FOUT(I),'  E=',FENG1(I),'  M=',F251(I)
      WRITE (6,'(//)')
      DO 68 I=1,IREM
   68 PRINT*,' I=',I,'  DP/DR=',DPHDRH(I),'  DNUE/DR=',DNUEDR(I),
     +       ' DM/DR=',DMSDR(I)
      WRITE (6,'(//)')
      END IF



      IF (IOUTPE /= 0) THEN

      OPEN (UNIT=64, FILE='pressure.dat', STATUS='unknown')

      WRITE (64,'(//)')
                       MM=2
      IF (IOUTPE == 1) MM=1
      IR=IOUTPE
      AUX5(1)=FST251(1)*P16
      AUX8(1)=FOUT(1)*EC
      AUX9(1)=FENG1(1)*EC
      DO 7871 I=1,IREM-1
      IF (I == IR) THEN
      AUX5(MM)=FST251(I)*P16
      AUX8(MM)=FOUT(I)*EC
      AUX9(MM)=FENG1(I)*EC
      IR      =IR+IOUTPE
      MM      =MM+1
      END IF
 7871 CONTINUE
      AUX5(MM)=FST251(IREM)*P16
      AUX8(MM)=FOUT(IREM)*EC
      AUX9(MM)=FENG1(IREM)*EC

      DO 7872 I=1,MM
      WRITE (64,686) AUX8(I),AUX9(I)
      IF (I == IDPREM .AND. IOUTPE == 1)
     +   WRITE (64,'('' ( ---- NEUTRON DRIP ---- )'')')
 7872 CONTINUE
      WRITE (64,'(//,''****************************'',//)')
      DO 788 I=1,MM
      WRITE (64,686) AUX5(I),AUX8(I)
      IF (I == IDPREM .AND. IOUTPE == 1)
     +   WRITE (64,'('' (---- NEUTRON DRIP ---- )'')')
  788 CONTINUE
      WRITE (64,'(//,''****************************'',//)')
      DO 687 I=1,MM
      WRITE (64,686) AUX5(I),AUX9(I)
      IF (I == IDPREM .AND. IOUTPE == 1)
     +   WRITE (64,'('' ( ---- NEUTRON DRIP ---- )'')')
  687 CONTINUE
  686 FORMAT (3X,E14.6,3X,E14.6)
      WRITE (64,'(//)')

               CLOSE ( 64 )
      END IF
!>>>> END OF WRITING TO OUTPUT FILE


!  added on 07/27/1996
      do 7010 i=1,IREM !-1
      aux5(i)=fst251(i)*P16
      aux9(i)=feng1(i)*EC*uf6
 7010 continue

!                              ****** ******
      call copy(aux5,aux9,irem,xstrob,ystrob)
!                              ^      ^
!                              r [km] e(r) (g/cm^3)

!__ Determine radial distance where mass density is eaual to e=2.4E14 g/cm**3
!   Don't do this for rotating white/strange dwarf calculations. Otherwise a
!     SIG 11 error will be produced when calling selspl.
! 05/09/2015
      if (iwhiteDwarf == 0) then ! for compact star calculations only
                                 ! skip this for white dwarf calculations
      
!                                          ****** 
       call selspl(ystrob,xstrob,irem,c,10,rinter,einter,10,eass)
!                                            ^      ^ input: crust density
!                                         r(e_nm) (km)
       end if
! crust thickness Delta_c
      rsurf        =tht1(mlf,5)
      tstrob(6,mlf)=rsurf-rinter
! end of supplement
 


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   INTEGRATE DIFFERENTIAL EQUATION OF THE METRIC FUNCTION  PHI(R)
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


      DO 500 N=2,IREM
      DMN=F251(N-1)
      DRN=FST251(N-1)
      DRN2=DRN*DRN
      DRN3=DRN2*DRN
      DPN=FOUT(N-1)

!.........................................
!>>>> PHI(R)
      HRK=FST251(N)-FST251(N-1)
      RNUED(N)=RNUED(N-1) + HRK * FN(DPN)
!.........................................

  500 CONTINUE


!   RENORMALIZATION OF THE METRIC FUNCTION
!   (ADDITIV CONSTANT)
      EM=1.-2.*F251(IREM)*DGRAV/FST251(IREM)
      IF (EM < 0.) THEN
      WRITE (6,'(//,4X,'' **  EM < 0 **'',///)')
      WRITE (6,NAME1)
      RENK=0.
      GO TO 800
      ELSE
      RENK=ALOG(EM)-RNUED(IREM)
      END IF

!__ RENORMALIZE  PHI(R)
      DO 700 N=1,IREM
  700 RNUED(N)=RNUED(N)+RENK



! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!              MOMENT OF INERTIA OF THE STAR (I)
!
! SIMPLIFIED TREATMENT: NEGLECTION OF THE DRAGGING OF THE INERTIAL
! FRAME. THE EXACT AND THE SIMPLIFIED TREATMENT ARE INVESTIGATED!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      WI=0.
  800 DO 600 J=1,IREM
      DRN=FST251(J)
      DRN2=DRN*DRN
      DRN4=DRN2*DRN2
      DMN=F251(J)
      DEN=FENG1(J)
      DPN=FOUT(J)
      W=RNUED(J)/2.
      IF (W > EXMAX) THEN
      W=EXMAX
      ELSE IF (W < EXMIN) THEN
      W=EXMIN
      END IF
      EXPOF=EXP(W)
      RAD=1.-2.*DGRAV*DMN/DRN
      IF (RAD <= 0.) RAD=1.E-08
  600 FST(J)=DRN4*(DEN+DPN)/(EXPOF*SQRT(RAD))


!__   PERFORM INTEGRATION
!
!                        R             4
!                   TM = I [D R']  (R')  FUNCTION(R'*R )
!                        0                            0

      EASS='HT1000: CALL SELSPL KKK       '
!                                              **
      CALL IGRN(FST251,FST,IREM,X251,Y251,IREM,TM,C,NHU,10,1,AUX51,EASS)
!     RETURN: TM
!__ COMPUTE  LOG_10 (I (g cm^2) )
      WI= ALOG10(2.*VPI*EC/3.)+5.*ALOG10(R0R0)+ALOG10(TM)-ALOG10(FMEVFM)
     +   -52.
      THT1(MLF,18)=WI



!__ COMPUTE THE CRITICAL ANGULAR VELOCITY  OMEGA(TILDE) = (M*G/R**3)**(1/2),
!   WHICH IS ROUGHLY THE ANGULAR VELOCITY AT WHICH EQUATORIAL SHEDDING (i.e.,
!   MASS LOSS) OCCURS (VERY RAPID ROTATION!), FROM THE PARAMETERS OF THE
!   N O N - ROTATING STAR (OV RESULTS):
      OGLUP=SQRT(F251(IREM)*DGRAV/FST251(IREM)) / FST251(IREM)
!__ Friedman, Ipser & Parker:
!     SET   OGLFR=0.659 * OGLTD,  WHICH IS EQUIVALENT TO SETTING
      OGLFR=24.*SQRT(F251(IREM)/(FST251(IREM)*R0R0*EEM18)**3) * 1.0E04
      OGLFR=OGLFR*R0R0*UF10
!__ RAY, DATTA:
      VRD=36.42*SQRT(F251(IREM)/(FST251(IREM)*R0R0*EEM18)**3) * 1.0E04
      OGLRD=SQRT(.27) * VRD
      OGLRD=OGLRD*R0R0*UF10

!>> CHOOSE  O N E  OF THE ABOVE COMPUTED FREQUENCIES: OGLUP,OGLFR,OGLRD
      OGLTD=OGLFR
      OGLT =OGLTD/R0R0
      
!   OMEGA(TILDE) IN  1/SEC
      OGLT =OGLT/UF10
      PRINT*,' OMEGA(TILDE)=',OGLT,' rad/s (--> begin solving HT)'
      THT1(MLF,6)=OGLT





! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! BEGIN OF THE RELATIVISITC STAR CALCULATIONS ACCORDING TO HARTLE AND
! THORNE (i.e., RELATIVISTIC ROTATING STAR TREATED IN THE THEORY OF
! GENERAL RELATIVITY)
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



! ******************************************************************
!           SOLVE DIFFERENTIAL EQUATION FOR OMEGA_BAR(R)
! ******************************************************************

!  STARTING VALUES
!  - UNCHANGED FOR ALL ROTATING STAR CALCULATIONS!
!  - CHOOSE  OGC=OMEGA_BAR[CENTER];  [OGC]=1/SEC
      OGC=6.0E03
      OGC=1.82342
      OGCFM=OGC*UF10
      OGCD=OGCFM*R0R0

      OGDR(1)=OGCD
      DRN=FST251(1)
      OGDR(2)=OGCD * (1.+16.*PI*A*DGRAV*(FENG1(1)+FOUT(1))*DRN*DRN)
      DOGDR(1)=0.
      DRN=FST251(2)
      DOGDR(2)=OGCD * 32.*PI*A*DGRAV*(FENG1(1)+FOUT(1))*DRN



      DO 905 I=3,IREM
      DRN=FST251(I-1)
      DMN=F251(I-1)
      RMGR=DMN*DGRAV/DRN
      W1=1.-2.*RMGR
      IF (W1 < 0.) THEN
         WRITE(6,'(///,'' HT1000: W1 < 0, SQRT COMPLEX!!'')')
         WRITE(6,'(//,'' W1='',E14.6,///)') W1
         STOP
      END IF
      DLNJDR=DGRAV * (DMN/DRN-DMSDR(I-1))/(W1*DRN) - DNUEDR(I-1)/2.
      C1R=FC1R(DRN,DLNJDR)
      C2R=FC2R(DRN,DLNJDR)
      HRK=FST251(I)-FST251(I-1)

!..............................................................
!   ETA(HAT;R) = D [OMEGA_BAR(R)] / D [R]
      DETDX    = - C1R*DOGDR(I-1) - C2R*OGDR(I-1)

!   SAVE SECOND DERIVATIVE OF  OMEGA_BAR(R)  FOR LATER USE
      D2OGDR(I-1) = DETDX

!   FIRST DERIVATIVE OF  OMEGA_BAR(R)
      DOGDR(I) = DOGDR(I-1) + HRK * DETDX
!   OMEGA_BAR(R)
      DOGDX    = DOGDR(I-1)
      OGDR(I)  = OGDR(I-1) + HRK * DOGDR(I-1)
!..............................................................

  905 CONTINUE


      N9111=6
      N9M1 =N9111-1
      DO 9107  I=1,N9111
      AUX5(I)   =FST251(I+N9M1)
 9107 AUX51(1,I)=DOGDR(I+N9M1)
      DO 9108  I=1,N9111
      E=FST251(I)
!                                       ***
      CALL INTER(1,AUX5,AUX51,1,1,N9111,DDO,E)
!     COPY EXTRAPOLATED VALUE OF  D[OMEGA]/D[R]
      DOGDR(I)=DDO(1)
 9108 CONTINUE

      N9111=10
      N9M1 =N9111-1
      DO 9111  I=1,N9111
      AUX5(I)   =FST251(I+N9M1)
 9111 AUX51(1,I)=D2OGDR(I+N9M1)
      DO 9112  I=1,N9111
      E=FST251(I)
!                                       ***
      CALL INTER(1,AUX5,AUX51,1,1,N9111,DDO,E)
!     COPY EXTRAPOLATED VALUE OF  D2[OMEGA]/D[R]2
      D2OGDR(I)=DDO(1)
 9112 CONTINUE


      N9111=5
      DO 8107  I=1,N9111
      AUX5(I)   =FST251(IREM-N9111+I)
 8107 AUX51(1,I)=D2OGDR(IREM-N9111+I)
      DO 8108  I=N9111,N9111
      E=FST251(IREM-N9111+I)
!                                       ***
      CALL INTER(1,AUX5,AUX51,1,1,N9111,DDO,E)
!     COPY EXTRAPOLATED VALUE OF  D2[OMEGA]/D[R]2
      D2OGDR(IREM-N9111+I)=DDO(1)
 8108 CONTINUE



!   COPY (SAVE CALCULATED VALUES FOR RESCALING PROCEDURE)
      DO 232 I=1,IREM
      OGAUX(I) =OGDR(I)
      DOGAUX(I)=DOGDR(I)
      D2OGAX(I)=D2OGDR(I)
  232 CONTINUE







! --------------------------------------------------------------------
!             BEGIN OF  S E L F - CONSISTENT CALCULATION
! NOTE: SELF-CONSISTENCY ONLY IF STAR MODELS ROTATING AT THEIR RESPEC-
!       TIVE GENERAL RELATIVISTIC KEPLER FREQUENCIES ARE CONSTRUCTED
      DO 5661 ISELFC=1,ISCMX7
! --------------------------------------------------------------------

!  THE EQUATIONS THAT ARE TO BE SOLVED SIMULTANEOUSLY ARE THOSE
!  DETERMINING:
!
!     (1)  OMEGA(STAR), as obtained from the dragging eqaution
!     (2)  OMEGA(STAR) = O_K_GR, self-consistency condition
!     (3)  J(DEFORMED STAR)



      IF (ISELFC == 1) THEN
!__ COMPUTE THE STAR'S ANGULAR MOMENTUM   J [OMEGA_C]
!   NOTE: THE FOLLOWING CALCULATIONS OF THE MOI AS WELL AS THE
!         STAR'S ANGULAR MOMENTUM REFER TO A SPHERICALLY SYMMETRIC
!         NEUTRON STAR. THE CASE OF A ROTATIONALLY DEFORMED STAR
!         IS TREATED BELOW.
      ZJD   =FST251(IREM)**4*DOGDR(IREM)/(6.*AGD)
      ZJD0  =ZJD
      ZJSAVE=ZJD


!__ COMPUTE THE STAR'S ANGULAR VELOCITY AS MEASURED BY A DISTANT
!   OBSERVER, AS GIVEN BY THE N O N - RESCALED FUNCTION omega_bar(r)
      OGLD=OGDR(IREM) + 2.*AGD*ZJD/FST251(IREM)**3
      OGL =OGLD/R0R0
!   OMEGA IN  1/SEC
      OGL =OGL/UF10

      OGLDSV=OGLD
      OGLSV =OGL

!   COMPUTE MOMENT OF INERTIA
!   NOTE: I  DOES  N O T  DEPEND ON OMEGA - AS LONG AS THE METRIC CORRECT
!         UP TO  O(OMEGA**3) IS CONSIDERED - AND CAN THEREFORE BE DETERMINED
!         FOR AN ARBITRARILY CHOSEN VALUE OF OMEGA (CHOOSE STARTING VALUE)
      ZISV =FST251(IREM)**4 * DOGAUX(IREM) / (6.*AGD*OGLDSV)
      ZISV0=ZISV
      END IF


!   INITIALIZE THE STAR'S ROTATIONAL FREQUENCY, i.e., OMEGA_STAR(NEW),
!   AND CALCULATE ITS PROPERTIES FOR THAT FREQUENCY
         K4A=1
         K4B=INPOG1
         IF (I5651 /= 1 .AND. ISCMX7 == 1) THEN
         K4A=MLF
         K4B=MLF
         END IF
         IF (I5651 == 1 .AND. ISCMX7 == 1) THEN
         K4A=MLF
         K4B=MLF
         END IF



! --------------------------------------------------------------------
! choose between different options: (1) self-consistent det. of Om_K
!                                   (2) M_rot and Omega  ==>  e_c
!                                   (3) e_c and Omega    ==>  M_rot
      DO 4006 K4=K4A,K4B
! --------------------------------------------------------------------

         K4OUT=K4
         OGV=RLOMG(K4)
         OGLT=OGV
         OGCOMS=OGLT
         TFREQ(K4,12)=OGLT
         OGLTD=OGLT*R0R0*UF10
         IF (ISELFC == 1) OGDIT(1)=OGLTD


      write (6,'(/,'' ----------------------------------------------'',
     +             ''-------------'')')
      print*,"Start solving Einstein's field equations selfconsistently"
      print*,'(--> find Kepler frequency of rotating stellar model)'
      write (6,'(  '' ----------------------------------------------'',
     +             ''-------------'')')
      print*,' ISELFC=',ISELFC,'   K4=',K4,'   e_c=',ec
      print*,' Number of radial grid points: IREM=',irem
      print*,' Omega=',oglt,' 1/s'
      print*,' Start rescaling scheme:'





!__ RESCALING OF THE FUNCTIONS
!                                 - OMEGA_BAR(R)
!                                 - D  [OMEGA_BAR(R)] / D [R]
!                                 - D2 [OMEGA_BAR(R)] / D [R]2
!                                 - J  [OMEGA_BAR(CENTER)]
!..........................
      SCF177=OGLTD/OGLDSV
!..........................
      DO 177 I=1,IREM
      IF(I == 1)
     +PRINT*,'                   rescale            omega_bar(r)'
      OGDR(I)  =OGAUX(I)  * SCF177
      IF(I == 1)
     +PRINT*,'                   rescale   d / d r  omega_bar(r)'
      DOGDR(I) =DOGAUX(I) * SCF177
      IF(I == 1)
     +PRINT*,'                   rescale  d2 / dr2  omega_bar(r)'
      D2OGDR(I)=D2OGAX(I) * SCF177
  177 CONTINUE



      IF (IOTOMG /= 0) THEN
      K908=1
      DO 908 I=1,IREM,K908
      IF (I == 1)
     +WRITE (6,'(///,'' HT1000: OMEGA_BAR(R), D[OG]/D[R], '',
     +           ''D2[OG]/D[R]2'',/)')
      PRTOG=OGDR(I)/(UF10*R0R0)
      WRITE (6,'('' I='',I4,''  OG='',E12.6,'' 1/S; D[OG]/D[R]='',
     +      E12.6,'' D2[OG]='',E12.6)') I,PRTOG,DOGDR(I),D2OGDR(I)
 908  CONTINUE
      END IF



!   RESCALING OF THE STAR'S ANGULAR MOMENTUM
      IF (ISELFC == 1) THEN
      PRINT*,'                   rescale            J'
      print*,' '
      ZJD=ZJSAVE * SCF177
      END IF




      IF (ISELFC == 1) THEN
!__ TEST:
!   RECOMPUTE THE STAR'S ANGULAR VELOCITY (AS MEASURED BY A DISTANT
!   OBSERVER) WHICH CORRESPONDS TO THE   R E S C A L E D  FUNCTION
!   OMEGA_BAR(R); USE  OMEGA_BAR( N E W ), J( N E W )
      OGLD=OGDR(IREM) + 2.*AGD*ZJD/FST251(IREM)**3
      OGL=OGLD/R0R0
!   OMEGA IN  1/SEC
      OGL=OGL/UF10
      WRITE (6,'(''  OGL [OMEGA_C]='',F8.2,'' 1/s'')') OGL

      DIFRC=ABS(OGL-OGLT)
      IF ( (DIFRC/OGL)  >=  1.E-05) THEN
      WRITE (6,'(//,'' OMEGA (RESCALED) DEVIATES FROM ITS NOMINAL'',
     +         '' VALUE!!'',/,''  DIFF='',E14.8,/,'' OGL='',E14.8,
     +         '' 1/SEC'',/,'' OGLT='',E14.8,'' 1/SEC'',/)') DIFRC,
     +       OGL,OGLT
      END IF

      ELSE
      OGLD=OGLTD
      OGL=OGLD/(R0R0*UF10)
      END IF



!__ ROTATIONAL ENERGY
      THT2(MLF,2)=ZJD*OGLD*A/2.


!__ COMPUTE ANGULAR VELOCITY OF THE LOCAL INERTIAL FRAMES
!   - OMEGA_CENTER
      OGSDC=OGLD-OGDR(1)
!   OMEGA IN  1/SEC
      OGSC =OGSDC/(R0R0*UF10)
!   OMEGA IN UNITS OF THE STAR'S ANGULAR VELOCITY
      OGSCUU=OGSDC/OGLD

!   - OMEGA_SURFACE
      OGSDS=OGLD-OGDR(IREM)
!   OMEGA IN  1/SEC
      OGSS=OGSDS/(R0R0*UF10)
!   OMEGA IN UNITS OF THE STAR'S ANGULAR VELOCITY
      OGSSUU     =OGSDS/OGLD
      THT1(MLF,7)=OGSSUU
      THT1(MLF,8)=OGSCUU

      IF (IOUT33 /= 0) THEN
      DO 333 I=1,IREM
      WQ=OGDR(I)/OGLD
      RQ=FST251(I)/FST251(IREM)
  333 PRINT*,' I=',I,' OG/OGL=',WQ,' r/R=',RQ
      END IF


!__ COMPUTE RADIUS OF GYRATION
      RGYD=SQRT(ZJD*A/(F251(IREM)*OGLD))
!     R[G] IN  cm
      RGY =RGYD*R0R0/EE13
!   R[G] / R
      RGDIVR=RGYD/FST251(IREM)

      THT1(MLF,4)=RGDIVR


!__ COMPUTE THE STAR'S MOMENT OF INERTIA
!   NOW: SPHERICAL STAR, DRAGGING TAKEN INTO ACCOUNT
!   - LOG ( I / (g cm**2) )
!>>>> WI= 2.*(ALOG10(RGYD)+ALOG10(R0R0/EE13)) + 33. + ALOG10(1.989) +
!>>>>+    ALOG10(F251(IREM))
      WI=ALOG10(ZJD)+ALOG10(EC*UF6)+5.*DLOG10(R0R0/EE13)-ALOG10(OGLD)
      THT1(MLF,19)=WI
      PRINT*,' log (I_old/g cm2)=', WI




! ******************************************************************
!     SOLVE SET OF COUPLED MONOPOLE EQUATIONS (L=0 CONTRIBUTION)
!     - MASS PERTURBATION FUNCTION  M_0(R)
!     - PRESSURE PERTURBATION FUNCTION  P*_0(R)
! ******************************************************************


!__ COMPUTE  D [EPS] / D [P]
!      NAX4=1
!      IF (IQM > 0) NAX4=45
!      XBA=AUX4(NAX4)
!      XBB=AUX4(IEP)
!      ARGINT='HT1000: PRESSURE GRID - X251  '
!                             ****    PRESSURE GRID  P(J)
!      CALL INTV(XBA,XBB,IREM,X251,0,ARGINT)
!      DO 25 J=1,IREM
!                                     ***         EPS(P=X251)
!      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,X251(J),1)
!      FST(J)=EIP(IKP1)
!   25 CONTINUE

!
!   COPY  EPS(P)  (RELABELING OF THE ELEMENTS!)
!      DO 967 I=1,IREM
!      KK      =IREM - (I-1)
!      X251(I) =FOUT(KK)
!  967 FST(I)  =FENG1(KK)
!  COMPUTE DERIVATIVE  D [E(int;P)] / D [P]


!    DERIVATIVE  D [EPS(R)] / D [R]
      IREMM1=IREM-1
      ISCD  =0
      DO 35 J=1,IREMM1
      Y251(J)=(FENG1(J+1)-FENG1(J)) / (FST251(J+1)-FST251(J))
      XNODE  =X251(J)
!    SAVE  D[E]/D[R]  FOR LATER USE IN SUBROUTINE MOIRDS
      DEHDRH(J)=Y251(J)
!                                         *****   RETURN: D [EPS_R] / D [R]
!      EASS='HT1000: CALL SELSPL LLL       '
!      CALL SELSPL(FST251,FENG1,IREM,C,10,YNODE,XNODE,20,EASS)
!      Y251(J)=YNODE

!>>> DERIVATIVE  D [EPS(P)] / D [P]
      Y251(J)=Y251(J)/DPHDRH(J)
      IF (Y251(J) < 0.) ISCD=ISCD+1


!TEST: COMPUTE DERIVATIVE ALTERNATIVELY BY CALCULATING THE DIFFERENTIAL
!      QUOTIENT OF THE ENERGY DENSITY IN THE STAR VS. PRESSURE
!      YG =(FENG1(J+1)-FENG1(J)) / (FOUT(J+1)-FOUT(J))
!      DYG=YG-Y251(J)
!      DYG=ABS(DYG)
!      PRINT*,' J=',J,' DE/DP=',Y251(J),' P=',FOUT(J),' DIFF=',DYG
!END OF TEST

!   SAVE  D[E]/D[P]  FOR LATER USE IN AUX51
      AUX51(1,J)=Y251(J)
   35 CONTINUE
      IF (IDPREM > 0) THEN
      IDRM1 = IDPREM - 1
      IDRP1 = IDPREM + 1
        DV2 = ( DEHDRH(IDRP1)-DEHDRH(IDPREM) ) /
     +        ( FST251(IDRP1)-FST251(IDPREM) )
      Y251(IDRM1) = DEHDRH(IDPREM) + DV2 *
     +              (FST251(IDRM1)-FST251(IDPREM))
      DEHDRH(IDRM1)  = Y251(IDRM1)
      Y251(IDRM1)    = Y251(IDRM1)/DPHDRH(IDRM1)
      AUX51(1,IDRM1) = Y251(IDRM1)
      END IF


      N9111=10
      N9M1 =N9111-1

      DO 7101  I=1,N9111
      AUX6(I)   =FOUT(I+N9M1)
 7101 AUX61(1,I)=AUX51(1,I+N9M1)
      DO 7102  I=1,N9111
      E=FOUT(I)
!                                       ***
      CALL INTER(1,AUX6,AUX61,1,1,N9111,DDO,E)
!     COPY EXTRAPOLATED VALUE OF  D[E]/D[P]
      AUX51(1,I)=DDO(1)
      Y251(I)   =AUX51(1,I)
 7102 CONTINUE
      N9111=5
      DO 7103  I=1,N9111
      AUX6(I)   =FOUT(IREMM1-N9111+I)
 7103 AUX61(1,I)=AUX51(1,IREMM1-N9111+I)
      DO 7104  I=N9111,N9111
      E=FOUT(IREM-N9111+I)
!                                       ***
      CALL INTER(1,AUX6,AUX61,1,1,N9111,DDO,E)
!     COPY EXTRAPOLATED VALUE OF  D[E]/D[P]
      AUX51(1,IREM-N9111+I)=DDO(1)
      Y251(IREM-N9111+I)   =AUX51(1,IREM-N9111+I)
 7104 CONTINUE


      N9111=10
      N9M1 =N9111-1
      DO 8101  I=1,N9111
      AUX6(I)   =FST251(I+N9M1)
 8101 AUX61(1,I)=DEHDRH(I+N9M1)
      DO 8102  I=1,N9111
      E=FST251(I)
!                                       ***
      CALL INTER(1,AUX6,AUX61,1,1,N9111,DDO,E)
!     COPY EXTRAPOLATED VALUE OF  D[E]/D[R]
      DEHDRH(I)=DDO(1)
 8102 CONTINUE
      N9111=5
      DO 8103  I=1,N9111
      AUX6(I)   =FST251(IREMM1-N9111+I)
 8103 AUX61(1,I)=DEHDRH(IREMM1-N9111+I)
      DO 8104  I=N9111,N9111
      E=FST251(IREM-N9111+I)
!                                       ***
      CALL INTER(1,AUX6,AUX61,1,1,N9111,DDO,E)
!     COPY EXTRAPOLATED VALUE OF  D[E]/D[R]
      DEHDRH(IREM-N9111+I)=DDO(1)
 8104 CONTINUE



!>>> CHECK VALUES OF CALCULATED DERIVATIVES
      IF (ISCD > IREM/5) THEN
         WRITE (6,'(///,'' D [EPS(P)] / D [P] < 0!!!  CHECK'',
     +              '' RESULTS!!!'',//)')
         DO 1717 I=1,IREM
 1717    WRITE (6,'('' I='',I4,'' DE/DP='',E12.6,'' P='',E12.6)')
     +         I,Y251(I),FOUT(I)
      ELSE
         DO 1718 I=2,IREM  ! 06/05/2014: used to be I=1,IREM
 1718    IF (Y251(I) < 0.) Y251(I)=Y251(I-1)
        END IF

!   COPY  D E / D P
      DO 3636 I=1,IREM
      DEDP(I)=Y251(I)
 3636 CONTINUE


      IF (IOUTES /= 0) THEN
      WRITE (6,'(///,'' HT1000: D[E]/D[P] VS. P, D[E]/D[R]'',/)')
      DO 36 I=1,IREM
      WRITE (6,'('' I='',I4,''  D[E]/D[P]='',E12.6,''  P='',E12.6,2X,
     +           ''D[E]/D[R]='',E12.6)') I,DEDP(I),FOUT(I),DEHDRH(I)
  36  CONTINUE
      END IF



!>>> COMPUTE DERIVATIVE  D [E(int;P)] / D [P]
!     D [E(int;P)]/D [P] = (D [E(int;P)]/D [R]) / (D [P]/D [R]))

!     COMPUTE DERIVATIVE  D [E(int;P)] / D [R]
      DO 1231 I=1,IREM
      XXX=FST251(I)
      EASS='HT1000: CALL SELSPL MMM       '
!                                         ***   RETURN: YYY=D[E_int]/D [R]
      CALL SELSPL(FST251,EINTHX,IREM,C,10,YYY,XXX,100,EASS)
      DEINDP(I)=YYY/DPHDRH(I)
 1231 CONTINUE




      DO 1808 I=1,IREM


!__ PRELIMINARY CALCULATIONS: FUNCTIONS   j(R),
!                                         D [ln (j(R))] / D [R],
!                                         D j(R)**2 / D [R]
      DRN=FST251(I)
      DRN2=DRN*DRN
      DRN3=DRN2*DRN
      DRN4=DRN3*DRN
      DEN=FENG1(I)
      DPN=FOUT(I)
      DMN=F251(I)

      RMGR=DMN*DGRAV/DRN
      W1=1.-2.*RMGR
      SQW1=SQRT(W1)

!   FUNCTION  j(R)
      RLJ(I)=SQW1 / EXP(RNUED(I)/2.)
!   FUNCTION  D [ln (j(R))] / D [R]
      DLNJ(I)=DGRAV * ( DMN/FST251(I)-DMSDR(I) ) / (W1*FST251(I))
     +        - DNUEDR(I)/2.
!   FUNCTION  D [j**2] / D [R]
      FDJ2=RLJ(I)*RLJ(I)*DNUEDR(I) +
     +     2.*DGRAV*(DMSDR(I) - DMN/FST251(I)) /
     +     (FST251(I) * EXP(RNUED(I)))
      DJ2DR(I)=-FDJ2
 1808 CONTINUE


!__ SOLVE THE COUPLED SET OF L=0 EQUATIONS
      DO 1908 I=1,IREM

      IF (I == 1) THEN

      RH    =FST251(I)
      RH2   =RH*RH
      RH3   =RH2*RH
      RH4   =RH2*RH2
      RH5   =RH4*RH
      DEN   =FENG1(I)
      DPN   =FOUT(I)
      DE1DP1=Y251(I)

      XRM00=VPI*(DEN+DPN)*(DE1DP1+2.) * (RLJ(I)*OGDR(I))**2 * A / 15.

!__ LIM  R -> 0
      RM0D(I)  =XRM00 * RH5
      DM0DR(I) =5. * XRM00 * RH4
      P0SD(I)  =(RLJ(I)*OGDR(I))**2 * RH2 / 3.
      DP0SDR(I)=2. * (RLJ(I)*OGDR(I))**2 * FST251(I) / 3.


      ELSE


      DRN =FST251(I-1)
      DRN2=DRN*DRN
      DRN3=DRN2*DRN
      DRN4=DRN3*DRN
      DEN =FENG1(I-1)
      DPN =FOUT(I-1)
      DMN =F251(I-1)

      PIP(IKP1)=DPN
      IF (PIP(IKP1) <= X251(1)) THEN
      EIP(IKP1)=Y251(1)
      ELSE IF (PIP(IKP1) >= X251(IREM)) THEN
      EIP(IKP1)=Y251(IREM)
      ELSE
!     INTERPOLATE                     ***         D [EPS] / D [P]
      CALL SELINT(X251,AUX51,1,1,IREM,EIP,DPN,1)
      END IF
      DE1DP1=EIP(IKP1)
!RRRRRR
      DE1DP1=Y251(I-1)
!RRRRRR

!__ COEFFICIENTS OF THE DIFFERENTIAL EQUATIONS
      BM1=VPI*A*DRN2*DE1DP1*(DEN+DPN)
      BM2=RLJ(I-1)*RLJ(I-1)*DRN4*DOGDR(I-1)*DOGDR(I-1) / (12.*DGRAV)
      BM3=-2.*DRN3*RLJ(I-1)*RLJ(I-1)*DLNJ(I-1)*OGDR(I-1)*OGDR(I-1)
     +    / (3.*DGRAV)

      BP01=-DGRAV*(1.+2.*VPI*A*DRN2*DGRAV*DPN) / (DRN2*W1*W1)
      BP02=-VPI*A*DGRAV*(DEN+DPN)*DRN / W1
      BP03=DRN3*RLJ(I-1)*RLJ(I-1)*DOGDR(I-1)*DOGDR(I-1) / (12.*W1)
      T1H=2.*DRN*RLJ(I-1)*OGDR(I-1)*(RLJ(I-1)*OGDR(I-1)+DRN*OGDR(I-1)
     +    *RLJ(I-1)*DLNJ(I-1)+DRN*RLJ(I-1)*DOGDR(I-1)) / W1
      T2H=2.*DGRAV*DRN*RLJ(I-1)*RLJ(I-1)*OGDR(I-1)*OGDR(I-1)
     +    *(DMSDR(I-1)-DMN/DRN) / (W1*W1)
      BP04=(T1H+T2H) / 3.

      HRK=FST251(I)-FST251(I-1)
!.................................................................
! MUE(0;R) AND P(0;R)
      DM0DR(I-1)=BM1*P0SD(I-1) + BM2 + BM3
      RM0D(I)=RM0D(I-1) + HRK*DM0DR(I-1)
      DP0SDR(I-1)=BP01*RM0D(I-1) + BP02*P0SD(I-1) + BP03 + BP04
      P0SD(I)=P0SD(I-1) + HRK*DP0SDR(I-1)
!.................................................................


      IF (I == IREM) THEN
      DM0DR(I) =BM1*P0SD(I)  + BM2          + BM3
      DP0SDR(I)=BP01*RM0D(I) + BP02*P0SD(I) + BP03 + BP04
      END IF

      END IF

      IF (IOTMP0 /= 0) PRINT*,' I=',I, ' MUE0=',RM0D(I),' P0=',P0SD(I)
 1908 CONTINUE



!__ TOTAL MASS OF ROTATING STAR WITH CENTRAL ENERGY DENSITY
!   EPSILON_C AND ANGULAR VELOCITY  OMEGA (=OGL)
      DELMH=RM0D(IREM) + ZJSAVE**2 * AGD*A / FST251(IREM)**3
!cc Begin of CHANGE
      DELMH=RM0D(IREM) + ZJD**2 * AGD*A / FST251(IREM)**3
!cc End of CHANGE
!   M + DELTA M
      RMTOTH=F251(IREM) + DELMH

! 06/05/2015
! For neutron stars with quark-hybrid composition only
      if (ioutse /= 0) then

         
! Mass contained in quark matter core
      delqmc=rm0d(k_qmc) + fst251(k_qmc)**5 * dogdr(k_qmc)**2 /
     +                     (36.*dgrav)
      rmtotq=rm_qmc + delqmc
      tabmm(1,mlf)=rmtotq

! Mass contained in mixed phase
      delmix=rm0d(k_mix) + fst251(k_mix)**5 * dogdr(k_mix)**2 /
     +                     (36.*dgrav)
      rmtotm=rm_mix + delmix - delqmc
      tabmm(2,mlf)=rmtotm

! Mass contained in nuclear liquid and crust
      rmtotc=rm_nlc + (delmh - delqmc - delmix)
      tabmm(3,mlf)=rmtotc

! Test
      tabmm(4,mlf)=rmtoth - (rmtotq + rmtotm + rmtotc)
      end if 
       
! STAR'S ANGULAR MOMENTUM  (IN UNITS OF  G * M**2)
      THT2(MLF,3)=ZJD * A / (DGRAV * RMTOTH**2)

! DELTA M / M(SPH)
      DMDIVM=DELMH/F251(IREM)
      THT1(MLF,11)=RMTOTH
      THT1(MLF,12)=DMDIVM

! ROTATING PROPER NS MASS
      THT2(MLF,14)=THT2(MLF,4) + DELMH



!__ COMPUTE FUNCTION  H(0;R)
!   - INSIDE THE STAR
      DO 951 I=1,IREM
      H0IN(I)=-P0SD(I)+(FST251(I)*OGDR(I))**2 / (3.*EXP(RNUED(I)))
  951 CONTINUE
!   - AT THE STAR'S SURFACE
      W1=1.-2.*F251(IREM)*DGRAV/FST251(IREM)
      H0ISF=-DELMH*DGRAV/(FST251(IREM)*W1)
     +      +(A*ZJD*DGRAV/FST251(IREM)**2)**2 / W1
      H0CON=H0ISF-H0IN(IREM)


!   - RENORMALIZE FUNCTION  H(0;R), AND COMPUTE  D [H(0;R)] / D [R]
!     FOR LATER USE IN  SUBROUTINE  MOIRDS
      DO 941 I=1,IREM
      H0IN(I)  = H0IN(I)+H0CON
      DRN=FST251(I)
      DRN2=DRN*DRN
      DREXP=3.*EXP(RNUED(I))
      DH0IN(I) = - DP0SDR(I) + DRN*DOGDR(I)**2 * (2.-DRN*DNUEDR(I))
     +             /DREXP + 2.*DRN2*OGDR(I)*DOGDR(I)/DREXP
  941 CONTINUE


!   CHECK CONTINUITY OF  H(0;R) AT THE STAR'S SURFACE
      EPSH0=1.E-05
      IF ( ABS(H0IN(IREM)-H0ISF)  >=  EPSH0) THEN
      WRITE (6,'(////,'' HT1000: FUNCTION  H(0;R) NOT CONTINUOUS'',
     +       '' AT THE STARS SURFACE:'',/,'' H(0,IN)='',E14.6,'' # '',
     +       '' H(0;OUTSIDE)='',E14.6,/,'' CONTINUE CALCULATION!!'',//
     +       )') H0IN(IREM), H0ISF
      END IF


! ************************************************************************
!     SOLVE THE SET OF COUPLED QUADRUPOLE EQUATIONS (L=2 CONTRIBUTION)
!     - FUNCTION  V(2:R)
!     - FUNCTION  H(2;R)
! ************************************************************************

!   INITIALIZE THE INTEGRATION CONSTANTS AT THE ORIGIN BY ASSIGNING THEM
!   AN ARBITRARY VALUE
      B2ARB=1.
      YZ=-B2ARB - VPI*(FENG1(1)+FOUT(1)) * (RLJ(1)*OGDR(1))**2*AGD/3.
      YN=ZPI * (FOUT(1)+FENG1(1)/3.) * AGD
      A2ARB=YZ/YN
      A2ARBH=-B2ARB / (ZPI*(FOUT(1)+FENG1(1)/3.) * AGD)


      DO 2908 I=1,IREM

      IF (I == 1) THEN

      RH =FST251(I)
      RH2=RH*RH
      RH3=RH2*RH
      RH4=RH2*RH2

      H200 =A2ARB*RH2
      V200 =B2ARB*RH4
      H200H=A2ARBH*RH2
      V200H=B2ARB*RH4


!__ LIM  R --> 0
      H2RIN(I)=H200
      V2RIN(I)=V200
      DH2DR(I)=2. * A2ARB * RH
      DV2DR(I)=4. * B2ARB * RH3

      H2HOM(I) =H200H
      V2HOM(I) =V200H
      DH2HOM(I)=2. * A2ARBH * RH
      DV2HOM(I)=4. * B2ARB  * RH3


      ELSE


      DRN =FST251(I-1)
      DRN2=DRN*DRN
      DRN3=DRN2*DRN
      DRN4=DRN3*DRN
      DEN =FENG1(I-1)
      DPN =FOUT(I-1)
      DMN =F251(I-1)

      RMGR=DMN*DGRAV/DRN
      W1  =1.-2.*RMGR

!__ COEFFICIENTS OF THE L=2 DIFFERENTIAL EQUATIONS
      BV21=-DNUEDR(I-1)
      BV22=1./DRN + DNUEDR(I-1)/2.
      BV23=-DRN3*DJ2DR(I-1)*OGDR(I-1)*OGDR(I-1)/3.
     +     + RLJ(I-1)*RLJ(I-1)*DRN4*DOGDR(I-1)*DOGDR(I-1)/6.

      BH21=-DNUEDR(I-1)
      BH22=DGRAV*(2.*VPI*A*(DEN+DPN)-4.*DMN/DRN3)
     +     / (W1*DNUEDR(I-1))
      BH23=-4./(DRN2*W1*DNUEDR(I-1))
      BH24=DRN3*RLJ(I-1)*RLJ(I-1)*DOGDR(I-1)*DOGDR(I-1)
     +     * (DNUEDR(I-1)*DRN/2. - 1. / (DRN*W1*DNUEDR(I-1))) / 6.
      BH25=-DRN2*DJ2DR(I-1)*OGDR(I-1)*OGDR(I-1)
     +     * (DNUEDR(I-1)*DRN/2. + 1. / (DRN*W1*DNUEDR(I-1))) / 3.



      HRK=FST251(I)-FST251(I-1)

!................  FUNCTIONS  V(2;R)  AND  H(2;R)  ...................

! COMPUTE  P A R T I C U L A R  SOLUTIONS!

      DV2DR(I-1)=BV21*H2RIN(I-1) + BV22*BV23
      V2RIN(I)  =V2RIN(I-1) + HRK*DV2DR(I-1)
      DH2DR(I-1)=(BH21+BH22)*H2RIN(I-1) + BH23*V2RIN(I-1) + BH24+BH25
      H2RIN(I)  =H2RIN(I-1) + HRK*DH2DR(I-1)


! COMPUTE  H O M O G E N E O U S  SOLUTIONS!

      DV2HOM(I-1)=BV21*H2HOM(I-1)
      V2HOM(I)   =V2HOM(I-1) + HRK*DV2HOM(I-1)
      DH2HOM(I-1)=(BH21+BH22)*H2HOM(I-1) + BH23*V2HOM(I-1)
      H2HOM(I)   =H2HOM(I-1) + HRK*DH2HOM(I-1)
!.....................................................................



      IF (I == IREM) THEN
      DV2DR(I) =BV21*H2RIN(I) + BV22*BV23
      DH2DR(I) =(BH21+BH22)*H2RIN(I) + BH23*V2RIN(I) + BH24+BH25
      DV2HOM(I)=BV21*H2HOM(I)
      DH2HOM(I)=(BH21+BH22)*H2HOM(I) + BH23*V2HOM(I)
      END IF

      END IF


      IF (IOTHV2 /= 0) THEN
      PRINT*,' I=',I,'  H2_PTC=',H2RIN(I),' V2_PTC=',V2RIN(I)
      PRINT*,'          H2_HOM=',H2HOM(I),' V2_HOM=',V2HOM(I)
      PRINT*,' '
      END IF
 2908 CONTINUE



!__ COMPUTE  G E N E R A L  SOLUTIONS OF THE  L=2  EQUATIONS
!   - H(2;R) = A' * H(2;HOMOG) + H(2;PARTIC)
!   - V(2;R) = A' * V(2;HOMOG) + V(2;PARTIC)
!   DETERMINE THE CONSTANTS OF INTEGRATION BY MATCHING THE SOLUTIONS
!   FOR THE INNER AND OUTER PART OF THE  L=2  EQUATIONS TO EACH OTHER!
      II=IREM
      SX0=F251(II)*DGRAV/FST251(II)
      SX1=1.+SX0

      T1=((AGD*ZJD/FST251(II)**2)**2/SX1 - H2RIN(II)) * DH2HOM(II)
     +   / H2HOM(II)
      T2=-(AGD*ZJD)**2 * (SX0/SX1-4.) / (FST251(II)**5*SX1)

      Y=1./SX0-1.
      IF (Y <= 1.) THEN
      WRITE(6,'(////,''  HT1000: ARGUMENT ERROR DETECTED! ASSOCIATED'',
     +      '' LEGENDRE POLYNOMIALS CANNOT BE'',/,'' CALCULATED'',//   ,
     +  '' DRN='',E14.6,'' DGRAV='',E14.6,'' SX0='',E14.6,'' Y='',E14.6,
     +          ///,'' ****** CONTINUE CALCULATION ******'',//,'' NOTE''
     +      ,'': WRONG RESULTS FOR  SPHERICAL SHAPE'',/,21X,''ECCENTRI''
     +      ,''CITY'',/,21X,''QUADRUPOLE MOMENT  OF THE RATATING STAR'')
     +      ') FST251(II),DGRAV,SX0,Y
      A2ARB =0.
      A2ARBP=0.
      GO TO 1441
      END IF

      Y2=Y*Y
      Z =1./(F251(II)*DGRAV)
      DIV=DQ22(Y,Y2,Z) - Q22(Y,Y2)*DH2HOM(II)/H2HOM(II)
      A2ARB =(T1+DH2DR(II)+T2) / DIV
      T1 =(AGD*ZJD/FST251(II)**2)**2 / SX1
      A2ARBP=(T1+A2ARB*Q22(Y,Y2) - H2RIN(II)) / H2HOM(II)


! THE  G E N E R A L  SOLUTION IS THEN GIVEN BY:
      DO 799 I=1,IREM
      H2RIN(I)= A2ARBP * H2HOM(I) + H2RIN(I)
      V2RIN(I)= A2ARBP * V2HOM(I) + V2RIN(I)
  799 CONTINUE



!__ COMPUTE NON-RADIAL MASS AND PRESSURE PERTURBATION FACTORS:
!   - FUNCTION  M(2;R)
!   - FUNCTION  P(2,*;R)

 1441 DO 971 I=1,IREM
      P2SD(I)=-H2RIN(I) - FST251(I)**2 * OGDR(I)**2
     +        / (3. * EXP(RNUED(I)) )
      RH=FST251(I)
      W1=1.-2.*F251(I)*DGRAV/RH
      RM2H(I)=RH*W1 * (-H2RIN(I) - RH**3*DJ2DR(I)*OGDR(I)**2/3.
     +                 +RH**4 * (RLJ(I)*DOGDR(I))**2 / 6.)
      IF (IOTMP0 /= 0) PRINT*,' I=',I,' P2=',P2SD(I),' M2=',RM2H(I)
  971 CONTINUE



!__ COMPUTE DERIVATIVES
!                           - D [M2;R)] / D [R]
!                           - D2 [j] / D [R]2
!                           - D2 [f2] / D [R]2
!   FOR LATER USE IN  SUBROUTINE   MOIRDS
      DO 9172 I=1,IREM
      DRN=FST251(I)
      DRN2=DRN*DRN
      DRN3=DRN2*DRN
      DRN4=DRN2*DRN2
      DMN=F251(I)
      RMGR=DMN*DGRAV/DRN
      W1=1.-2.*RMGR
      EK1=W1+2.*DGRAV*(DMN/DRN - DMSDR(I))
      EK2=-H2RIN(I) - DRN3*DJ2DR(I) * OGDR(I)**2 / 3.
     +              + DRN4 * (RLJ(I)*DOGDR(I))**2 / 6.
      DDJ=RLJ(I)*DLNJ(I) * ( (DEHDRH(I)+DPHDRH(I))/(FENG1(I)+FOUT(I))
     +                      -DLNJ(I) - DNUEDR(I) + 1./DRN )
      D2JDR2(I)=DDJ
      D2J22(I)=2. * ( DJ2DR(I)*DLNJ(I) + RLJ(I)*(D2JDR2(I)-RLJ(I)*
     +                                           DLNJ(I)*DLNJ(I)) )
      EK3=-DH2DR(I) - DJ2DR(I)*DRN2*OGDR(I) * (OGDR(I)+2.*DRN*
     +                                         DOGDR(I)/3.)
     +    +2.*DRN3 * (RLJ(I)*DOGDR(I))**2 * (1.+DRN*DLNJ(I)/2.) / 3.
     +    -DRN3 * ( OGDR(I)**2 * D2J22(I) - DRN*RLJ(I)**2 * DOGDR(I)*
     +                                 D2OGDR(I)/3. )
      DRM2H(I)=EK1*EK2 + DRN*W1*EK3
 9172 CONTINUE




!__ ROTATIONAL DEFORMATION OF THE STAR

!   - FUNCTIONS ZETA(0;R), ZETA(2;R)
      DO 847 I=1,IREM
      Z0AR(I)=-P0SD(I)*(FENG1(I)+FOUT(I)) / DPHDRH(I)
  847 Z2AR(I)=-P2SD(I)*(FENG1(I)+FOUT(I)) / DPHDRH(I)


!   ZETA0 AND ZETA2  AT  STAR'S SURFACE (P=0)
      IISF=IREM
      Z0H   =-P0SD(IISF)*(FENG1(IISF)+FOUT(IISF)) / DPHDRH(IISF)
      Z2H   =-P2SD(IISF)*(FENG1(IISF)+FOUT(IISF)) / DPHDRH(IISF)
      R1SURF=FST251(IISF)+Z0H
      R2SURF=Z2H+FST251(IISF)*(V2RIN(IISF)-H2RIN(IISF))

!   DELTA R / R
      DRDIVR=Z0H/FST251(IISF)
      THT1(MLF,9)=DRDIVR


!   ZETA0 AND ZETA2  AT STAR'S INNER SURFACE AT DENSITY  RHO_drip (P=P_drip)
      Z0HD=0.
      Z2HD=0.
      R1DRIP=0.
      R2DRIP=0.
      IF (IDPREM > 0) THEN
      IIIN=IDPREM
      Z0HD=-P0SD(IIIN)*(FENG1(IIIN)+FOUT(IIIN)) / DPHDRH(IIIN)
      Z2HD=-P2SD(IIIN)*(FENG1(IIIN)+FOUT(IIIN)) / DPHDRH(IIIN)
      R1DRIP=FST251(IIIN)+Z0HD
      R2DRIP=Z2HD+FST251(IIIN)*(V2RIN(IIIN)-H2RIN(IIIN))
      END IF

      DO 961 I=0,NI
      X    =THT(I)
      XDEG =X*(180./PI)
      POLY2=(3.*COS(2.*X)+1.)/4.

!  I)  CONSIDER SURFACE OF CONSTANT GIVEN DENSITY IN A PARTICULAR COORDINATE
!      SYSTEM; THIS INVARIANT PARAMETRIZATION OF THE SURFACE LEADS FOR THE
!      RADIUS OF THE ROTATING STAR TO:
!   STAR'S SURFACE AT ZERO DENSITY
      RTH(I,2)=(R1SURF+R2SURF*POLY2)*P16
!   STAR'S SURFACE AT NEUTRON DRIP DENSITY
      RTH(I,4)=(R1DRIP+R2DRIP*POLY2)*P16

!  II) CONSIDER SURFACE OF CONSTANT GIVEN DENSITY IN A PARTICULAR COORDINATE
!      SYSTEM; THIS LEADS FOR THE RADIUS OF THE ROTATING STAR TO:
      RTH(I,3)=(R1SURF+Z2H*POLY2)*P16
  961 CONTINUE



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!  compute R_polar vs. R_equator; save outcome in file 'radius.dat'
      if (iradot /= 0) then

      open (unit=70, file='radius.dat', status='unknown')

      i847=1
      if (idprem > 0) i847=2
      do 3113 jj=1,i847
      do 3119 i=0,ni
      if (jj == 1) then
      rath=rth(i,2)
      else
      rath=rth(i,4)
      end if
      xv=rath*sin(tht(i))
      yv=rath*cos(tht(i))
 3119 write (70,3118) xv,yv
 3118 format (3x,e14.6,3x,e14.6)
      write (70,3117)
 3117 format (' ------------------------------------------')
 3113 continue
      close ( 70 )
      end if
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv



      

!__ COMPUTE THE STAR'S ECCENTRICITY DUE TO ROTATION
      ECCENT= -3.*(V2RIN(IISF)-H2RIN(IISF)+Z2H/FST251(IISF))
      IF (ECCENT < 0.) THEN
      WRITE (6,'(////,'' HT1000: ECCENTRICITY DEFINED ACCORDING'',
     +       '' TO HARTLE AND THORNE 1968 CANNOT BE'',/,12X,''CALC'',
     +       ''ULATED:'',/ ,12X,''V2='',E14.6,'' H2='',E14.6,'' ZE'',
     +       ''TA2(R)='',E14.6,/,12X,''R='',E14.6,'' IIIN='',I5,
     +       '' IISF='',I5)') V2RIN(IISF),H2RIN(IISF),Z2H,FST251(IISF
     +       ),IIIN,IISF
      ECCENT=0.
      ECHT68=SQRT(ECCENT)
      ELSE
!     ECCENTRICITY ACCORDING TO  [HT68]
      ECHT68=SQRT(ECCENT)
                            END IF
      THT1(MLF,13)=ECHT68



!   INVARIANT PARAMETRIZATION
      RDR=1. - (RTH(0,2)/RTH(NI,2))**2
      IF (RDR < 0.) RDR=0.
      ECFIP1=SQRT(RDR)
      THT1(MLF,14)=ECFIP1
!     SURFACE
      THT1(MLF,15)=RTH(0,2)
      THT1(MLF,16)=RTH(NI,2)
!     DRIP
      THT2(MLF,19)=RTH(0,4)
      THT2(MLF,20)=RTH(NI,4)

!   PARAMETRIZATION IN A PARTICULAR COORDINATE SYSTEM
      RDR               =1. - (RTH(0,3)/RTH(NI,3))**2
      IF (RDR < 0.) RDR=0.
      ECFIP2            =SQRT(RDR)
      THT1(MLF,20)      =ECFIP2
      THT1(MLF,21)      =RTH(0,3)
      THT1(MLF,22)      =RTH(NI,3)


!__ COMPUTE THE STAR'S QUADRUPOLE MOMENT
!     DEFINED AS THE COEFFICIENT OF  R**(-3) * P (COS(THETA)) TERM IN THE
!     NEWTONAIN POTENTIAL                       2

      QH=8. * A2ARB * ( F251(IISF)*DGRAV / FST251(IISF) )**2 / 5.
     +   +    AGD   * ( ZJD              / FST251(IISF) )**2
      THT1(MLF,17) = QH
      tquad(MLF,1) = QH
      tquad(MLF,2) = QH * tht1(MLF,5)**2 * tht1(MLF,10)*rsunkm  ! Q in km^3


      IF (ICMOI /= 0) THEN
! *********************************************************************
!     COMPUTE THE MOMENT OF INERTIA OF A ROTATIONALLY  D E F 0 R M E D
!     STAR. THE MOMENT OF INERTIA IS CORRECT UP TO OMEGA(STAR)**3
! *********************************************************************
      PRINT*,' '
      PRINT*,' --> COMPUTE MOMENT OF INERTIA (GENERALIZED H-METHOD)'


      RAD1=0.
      RAD2=FST251(IREM)
!   RETURN                                      **** **** **** ****
      CALL MOIRDS(RAD1,RAD2,EC,QH,AUX3,AUX4,NHU,WI0H,WI2H,WJ0H,WJ2H,
!                 *** *** **
     +            WIH,WJH,WI,C,IREM,OGLTD,RSPH,SUMW1,ZJD)
!   SAVE  log(I_def_star/g cm2)
      THT1(MLF,27)=WI
      PRINT*,' log(I_def_star/g cm2)=',WI
      WI0GCM= ALOG10(WI0H)+5.*DLOG10(R0R0/EE13)+ALOG10(EC*UF6)
      PRINT*,' log(I_new/g cm2)     =',WI0GCM
      
!cc     NOTE: SET  ZISV0 = WI0H
      ZISV0=WI0H
      ZISV =WIH
      ZID0 =WI0H
      ZJD  =WJH
      ZJD0 =WJ0H
!cc     NOTE: SET ZISV = ZISV0 AND ZJD=WJ0H
      zisv =zisv0
      zjd  =wj0h
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WI= ALOG10(ZJD0)+ALOG10(EC*UF6)+5.*DLOG10(R0R0/EE13)
     +   -ALOG10(OGLD)
      DIFFI0=ABS(THT1(MLF,19)-WI)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      print*,' Moment of inertia: I_0(HAT)    =',WI0H
      print*,'                    I_2(HAT)    =',WI2H
      print*,'                    I_TOTAL(HAT)=',WIH
      END IF





!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!      COMPUTE GENERAL RELATIVISTIC KEPLER FREQUENCY OF THE STAR MODEL
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      RDCBEP=THT1(MLF,15)*EE18/R0R0
      RDCBEE=THT1(MLF,16)*EE18/R0R0
      W1PO  =1.-2.*RMTOTH*DGRAV/RDCBEP
      W1EQ  =1.-2.*RMTOTH*DGRAV/RDCBEE

      IF (W1PO < 0. .OR. W1EQ < 0.) THEN
      WRITE (6,'(/,'' METRIC FUNCTION CANNOT BE CALCULATED'',/,
     +             '' W1PO='',E14.6,'' W1EQ='',E14.6,/)') W1PO,W1EQ
      GO TO 4007
      END IF

!>>> COMPUTE  h_0(r=R_pole)  and   h_2(r=R_pole)
!    NOTE:  R_pole < R_sph
      RD=RDCBEP
      DO 1022 LH0H2=1,IREM
      LREM=LH0H2
      RS  =FST251(LH0H2)
      IF (RS >= RD) GO TO 1023
 1022 CONTINUE
 1023          CONTINUE

!>>> FUNCTION  H_0(R)  INSIDE THE STAR
      H0RCBP=H0IN(LREM)

!>>> QUADRUPOLE FUNCTION  H_2(R)  INSIDE THE STAR
      H2RCBP=H2RIN(LREM)

!>>> QUADRUPOLE FUNCTION  V_2(R)  INSIDE THE STAR
      V2RCBP=V2RIN(LREM)

!>>> METRIC FUNCTIONS  PHI_POLE  AND  PHI_EQUATOR
      PHIRCP=ALOG(W1PO)
      PHIRCE=ALOG(W1EQ)

!>>> omega_bar  OUTSIDE THE STAR; D [omega] / D [R]
!RRR NOTE: SET R=R_equator
      ODCBEE=OGLD - 2.*AGD*ZJD / RDCBEE**3
      DOGDCB=6.*AGD*ZJD / RDCBEE**4


!>>> FUNCTION  H_0(R)  OUTSIDE THE STAR
      H0RCBE=-DELMH*DGRAV/(RDCBEE*W1EQ)
     +       +(A*ZJD*DGRAV/RDCBEE**2)**2 / W1EQ

!>>> QUADRUPOLE FUNCTION  H_2(R)  OUTSIDE THE STAR
      COMT1 =(AGD*ZJD/RDCBEE**2)**2
      AMGR  =F251(IREM)*DGRAV/RDCBEE
      AQ    =1./AMGR-1.
      AQ2   =AQ*AQ
      H2RCBE=COMT1*(1.+1./AMGR) + A2ARB*Q22(AQ,AQ2)
      QH    =A2ARB*Q22(AQ,AQ2)
      H2NEGL=COMT1*(1.+1./AMGR)
      QHDH2 =QH/H2NEGL


!>>> QUADRUPOLE FUNCTION  V_2(R)  OUTSIDE THE STAR
      V2RCBE=-COMT1 + 2.*A2ARB*AMGR*Q21(AQ,AQ2)/SQRT(W1EQ)
      QV    =2.*A2ARB*AMGR*Q21(AQ,AQ2)/SQRT(W1EQ)
      V2NEGL=-COMT1
      QVDV2 =QV/V2NEGL
      IPQQ  = 0
      IF (IPQQ /= 0) THEN
      PRINT*,' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx',
     +       'xxxxxxxxxxxxxxx'
      PRINT*,' h_2=',h2rcbe,' h0_2=',h2negl,' Q_h/h0_2=',QHDH2
      PRINT*,' v_2=',v2rcbe,' v0_2=',v2negl,' Q_v/v0_2=',QVDV2
      PRINT*,' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx',    
     +       'xxxxxxxxxxxxxxx'
         END IF



      VSREQT=-1.
      OGDKGR=-1.
      RKV   = 1. + 2.*(H0RCBE-H2RCBE/2.)
      RKV2H2=(1. + 2.*(H0RCBE-H2RCBE/2.)) / (1. - (V2RCBE-H2RCBE))
      IF (RKV <= 0.) THEN
      WRITE (6,'(////,''  HT1000: KEPLERIAN VELOCITY  CANNONT BE '',
     +           ''CALCULATED;  SQRT < 0 !!'',/,'' RKV='',E14.6,///)')
     +      RKV
      GO TO 5787
      ELSE IF (RKV2H2 <= 0.) THEN
      WRITE (6,'(////,''  HT1000: KEPLERIAN FREQUENCY CANNONT BE '',
     +        ''CALCULATED;  SQRT < 0 !!'',/,'' RKV2H2='',E14.6,///)')
     +      RKV2H2
      GO TO 5787
      END IF


      DMETFC=  2.*F251(IREM)*DGRAV / (RDCBEE**2 * W1EQ)
      AWZ   =  EXP(PHIRCE)*DMETFC / (2.*RDCBEE)
      BWZ   =  ( RDCBEE*DOGDCB/2. )**2
      IF ((AWZ+BWZ) < 0.) THEN
      AWZ=0.
      BWZ=0.
      WRITE (6,'(//,'' NOTE: APPROXIMATE CALCULATION OF O_K_GR;'',
     +              ''       SQRT CANNOT BE PERFORMED'',/,
     +              ''       ARGUMENT < 0'',/,
     +              ''       SET ARGUMENT=0 AND CONTINUE CALCULATION'',
     +           //)')
      END IF

!//////////////////// COMPUTE Omeg_K  (Route I) \\\\\\\\\\\\\\\\\\\\\\\
!  O_K: MONOPOLE AND QUADRUPOLE FUNCTIONS ARE NEGLECTED
!     OGDKGR=  -RDCBEE*DOGDCB/6. + SQRT(AWZ+BWZ)
      OGDKGR=  -RDCBEE*DOGDCB/2. + SQRT(AWZ+BWZ) + 2.*AGD*ZISV*OGLTD
     +         /RDCBEE**3

!  second solution of O_K, V_eq < 0
      OGMINS=  -RDCBEE*DOGDCB/2. - SQRT(AWZ+BWZ) + 2.*AGD*ZISV*OGLTD
     +         /RDCBEE**3
      THT2(MLF,26)=OGMINS/(UF10*R0R0)


!  THE STAR'S VELOCITY AT THE EQUATOR (IN UNITS OF C) CALCULATED FOR
!  ROTATION AT THE GENERAL RELATIVISTIC KEPLER FREQUENCY
      VSREQT=  (-RDCBEE**2*DOGDCB/2. + RDCBEE*SQRT(AWZ+BWZ))
     +         / EXP(PHIRCE/2.)

!  second solution of quadratic equation for V_eq, V_eq <0
      VSMINS=  (-RDCBEE**2*DOGDCB/2. - RDCBEE*SQRT(AWZ+BWZ))
     +         / EXP(PHIRCE/2.)
      THT2(MLF,25)=VSMINS
!/////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


!RRR OUTPUT:
      WRITE(6,'(/)')
      PRINT*,'=======================================================',
     +       '======'
      PRINT*,'       BEGIN CALCULATION OF O_K_GR: ROUTE I     '
      PRINT*,'=======================================================',
     +       '======'
      PRINT*,' C U R V A T U R E  OF SPACE-TIME AT STARS EQUATOR:'
      EPMN=RDCBEE/EXP(PHIRCE/2.)
      PRINT*,' exp(psi-nu)/r_0=',EPMN
      DVN=DMETFC/2.
      DVP=1./RDCBEE
      PRINT*,' dnue/dr=',DVN
      OEQ =2.*AGD*ZJD / RDCBEE**3
      DOEQ=-6.*AGD*ZJD / RDCBEE**4
      PRINT*,' omega(R_eq)=',OEQ,'  d omega/d r=',DOEQ
      T1 =DVN/DVP
      T22=(DOGDCB*EPMN / (2.*DVP))**2
      WZ =T1+T22
      PRINT*,' T1=',T1,'  T22=',T22
      PRINT*,' sqrt(t1+t22)=',SQRT(WZ)
      P11=DVN/DVP
      P22=-SQRT(T22)
      PRINT*,' (dnue/dr)/(dpsi/dr) =',P11
      PRINT*,' (domega/dr exp(psi-nu)/(2 dpsi/dr)=',P22
      WRITE(6,'(/)')


      WRITE(6,'(/,'' --------------------------------------------'')')
      WRITE(6,'(  '' |                                          |'')')
      WRITE(6,'(  '' |  V(eq)/c='',F4.2,''  O(K;GR)='',F7.1,'' rad/s'',
     +            ''     |'')') VSREQT,OGDKGR/(UF10*R0R0)
      WRITE(6,'(  '' |                                          |'')')
      WRITE(6,'('' --------------------------------------------'',/)')

!__ COMPUTE KEPLERIAN FREQUENCY  W I T H O U T  NEGLECTING THE MONOPOLE
!   AND QUADRUPOLE FUNCTIONS
      EPH   =EXP(PHIRCE/2.)
      EPHINV=1./EPH
      TZ    =1.-V2RCBE+H2RCBE
      TN    =1.+2.*(H0RCBE-H2RCBE/2.)
      TZDTN =TZ/TN
      TWZ   =SQRT(TZDTN)

!   METRIC FUNCTIONS  nue(H)_eq  AND  psi(H)_eq  (Hartle's conventions)
      REQKM  = THT1(MLF,16)  ! equatorial radius in km
      ANUEH  = PHIRCE + ALOG(TN)
      APSIH  = 2.*ALOG(REQKM) + ALOG(TZ)


!>>> EXP [PSI - NUE]
      EPMN=R0R0 * RDCBEE * TWZ / EPH

!>>> EXP [NUE - PSI]
      ENMP=1. / EPMN

      VV=(AGD*ZJD)**2/RDCBEE**5
      BK1=DELMH*DGRAV/(RDCBEE**2 * W1EQ) - VV/W1EQ
      BK2=1.+2.*AMGR/W1EQ

!>>> D [H_0] / D [R]
      DH0OSE=BK1 * BK2 - 3. * VV / W1EQ

      BR1=-VV * (4.-1./AMGR) * (1.+1./AMGR)
      DZDR=1./(AMGR*RDCBEE)

!>>> D [H_2] / D [R]
      DH2OSE=BR1 + A2ARB * DQ22(AQ,AQ2,DZDR)

      BK1=1./RDCBEE + 1./(2.*W1EQ) - DQ21(AQ,AQ2,DZDR)

!>>> D [V_2] / D [R]
      DV2OSE=4.*VV - 2.*A2ARB*AMGR*BK1 / SQRT(W1EQ)

!>>> D [NUE] / D [R]
      DVN=DMETFC/2. + (DH0OSE - DH2OSE/2.) / RKV

!>>> D [PSI] / D [R]
      DVP=(1.- (2.* (V2RCBE-H2RCBE) + RDCBEE * (DV2OSE-DH2OSE))/2.)
     +    / (RDCBEE * (1.-V2RCBE+H2RCBE))


!   RADIAL DERIVATIVES OF METRIC FUNCTIONS
!   note: capital  H  in the following equations refers to Hartle's
!         notation
!                                           d nue(H)_eq / d r
!                                           d psi(H)_eq / d r
!rrrrrrrrrrrrrrrrrrrrrrrrrrr replace  ..,10)  -->  ..,11)
      ADPHIH = 2.*THT1(MLF,10)*RSUNKM / (REQKM**2 * EPH**2)
!rrrrrrrrrrrrrrrrrrrrrrrrrr  end of replace
      ADNUEH = ADPHIH + 2.*(DH0OSE - DH2OSE/2.) / (RKV*R0R0/EE18)

      AXDPSI= (1.- V2RCBE + H2RCBE - REQKM*(DV2OSE-DH2OSE)/
     +                                     (2.*R0R0/EE18))
     +        / (REQKM * (1.-V2RCBE+H2RCBE))
      ADPSIH = 2.*AXDPSI

!   SAVE metric functions as defined by [FIP86]
      ANUE         = ANUEH  / 2.
      ANUEFIP      = ANUE
      APSI         = APSIH  / 2.
      ADNUE        = ADNUEH / 2.
      ADPSI        = ADPSIH / 2.
      TAPPX(MLF,1) = ANUE
      TAPPX(MLF,2) = APSI
      TAPPX(MLF,3) = ADNUE
      TAPPX(MLF,4) = ADPSI

!   COMPUTE FUNCTIONS  f  and  g   (i.e., prefactors occuring in the
!   empirical formula for the general relativistic Kepler frequency)
      REQKM3 = REQKM**3
      REQKM4 = REQKM3*REQKM
      RM     = THT1(MLF,11)
      FFFAUX = (ADNUE/ADPSI)*EXP(2.*(ANUE-APSI)) * REQKM3 /
     +         (RM*RSUNKM)
      IF (FFFAUX < 0.) FFFAUX=0.
      FFF    = SQRT(FFFAUX)

      inot = 10
      if (im05 >= 0) inot = 0
      if (inot /= 0) then
      PRINT*,' _IDEFST_  -->  computing I_def (in km**3)'
 
! 06/29/2006 idefopt 
!   idefopt=1 >< compute I_total of the rotationally deformed star
!   idefopt=2 >< compute both I_crust and I_total 
      idefopt = 2

      if (iwhiteDwarf == 0) then ! for neutron and quark stars only (05/22/2015)
         
! 06/29/2006 ------------------------------------  IDPREM determined above. Set
!                                               |  IDPREM=0 if MOI of crust is
!                                              \ / not to be computed
      CALL IDEFST(RJI,NI,NJ,NIP1,NJP1,THT,IREM,idprem,NHU,C,AX1,AX2,
     +            RTH,MLF,EC,OGLTD,idefopt,R_90,e_90,R_00,e_00)
!     RETURN: COMMON/TABHT1/....,THT2(.,.)
      
      RAX   = THT2(MLF,24) + DLOG10(GCM2KM) - alog10(ufe30)
      RMOI  = 10**RAX
      RAX   = THT1(MLF,19) + DLOG10(GCM2KM) - alog10(ufe30)

      else                      ! for rotating white dwarfs only

      RAX   = THT1(MLF,19) + DLOG10(GCM2KM) - alog10(ufe30)
      RMOI  = 10**RAX  ! set I_rotating = I_spherical 
      end if
      
! I_sph  in units of  km**3
      RMOISS = 10**RAX



!rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr boc
!       rmoi=rmoiss
!rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr eoc

      PRINT*,' I_def=',RMOI,' km3','  I_sph=',RMOISS,' km3'
      ZIAUX  = THT2(MLF,24) - ALOG10(EC*UF6) - 5.*DLOG10(R0R0/EE13)
      ZISVDD = 10**ZIAUX
      ZJDDD  = ZISVDD*OGLTD
      end if

      GG1  = 1. + 3.*RMOI/(REQKM4*ADPSI) - 2.*RMOI/REQKM3
      GG2  = 3.*RMOI / (REQKM4*ADPSI)

      REQ  = THT1(MLF,5)
!rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr boc
      REQ  = REQKM
!rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr eoc
      REQ3 = REQ**3
      REQ4 = REQ**4
      GG1  = 1. + 3.*RMOI/(REQ4*ADPSI) - 2.*RMOI/REQ3
      GG2  = 3.*RMOI / (REQ4*ADPSI)

! 08/18/1993: Begin of changes 
      tt1 = ogsds*tht1(mlf,16)*ee18/r0r0
      tt2 = (tht1(mlf,16)*ee18/r0r0) / (rm*dgrav)
      ttz = 1. + tt1*tt1 * (1.+tt2/2.) / 2.
      ttn = 1. - tt1*tt1 * (1.+tt2/4.) / 2.
      tx3 = 1.5 * (ogsds/ogdkgr) * ttz / ttn
      gg1= 1. + tx3 - ogsds/ogld
      gg2= tx3
! 08/18/1993: End of changes


      GDIF = GG1*GG1 - GG2*GG2
      IF (GDIF > 0.) THEN
      GGG  = 1. / SQRT(GDIF)
      ELSE
      GGG  = 0.
      END IF

      FTG          = FFF * GGG
      TAPPX(MLF,5) = FFF
      TAPPX(MLF,6) = GGG
      TAPPX(MLF,7) = FTG
      TAPPX(MLF,9) = RMOI
      TAPPX(MLF,8) = 3.642339E+05 * FTG * SQRT(RMTOTH/REQKM3)



      WZ=DVN/DVP + (DOGDCB*EPMN / (2.*DVP*R0R0))**2
      IF (WZ < 0.) THEN
      WZ=0.
      WRITE (6,'(//,'' NOTE: E X A C T  CALCULATION OF O_K_GR;'',
     +              ''       SQRT CANNOT BE CALCULATED'',/,
     +              ''       ARGUMENT < 0'',/,
     +              ''       SET ARGUMENT=0 AND CONTINUE CALCULATION'',
     +           //)')
      END IF



!/////////////////// COMPUTE Omeg_K  (Route II) \\\\\\\\\\\\\\\\\\\\\\\
!  THE STAR'S VELOCITY AT THE EQUATOR (IN UNITS OF C) CALCULATED FOR
!  ROTATION AT THE  KEPLER FREQUENCY
      VSDM  =-DOGDCB * EPMN / (2.*DVP*R0R0) + SQRT(WZ)
      VSMODI=VSDM
      VSDMH =VSDM * R0R0
!  O_K:  MONOPOLE AND QUADRUPOLE FUNCTIONS  N O T  NEGLECTED
      OGDKDM=ENMP * VSDMH + 2. * OGLD * AGD * ZISV / RDCBEE**3
!>>>> ALTERNATIVELY CAN BE USED:
!>>>> OGDKDM=ENMP * VSDMH + OGLD * (1.-ODCBEE/OGLD)

      AXFWBW=(2.*OGLD*AGD*ZISV/RDCBEE**3) * EPMN/R0R0
      FRP   =EPH

!  FRACTIONAL REDSHIFT IN  B A C K W A R D  DIRECTION (DEFINITION OF FIP)
      X1=(1.+VSMODI) / (1.-VSMODI)
      IF (X1 < 0.) X1=0.
      X1=SQRT(X1)
      X2=1.+AXFWBW
      ZBW=X1 / (X2 * FRP) - 1.
      THT2(MLF,11)=ZBW

!  FRACTIONAL REDSHIFT IN  F O R W A R D  DIRECTION (DEFINITION OF FIP)
      X1=(1.-VSMODI) / (1.+VSMODI)
      IF (X1 < 0.) X1=0.
      X1=SQRT(X1)
      X2=1.-AXFWBW
      ZFW=X1 / (X2 *FRP) - 1.
      THT2(MLF,12)=ZFW
!//////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



!RRR OUTPUT:
      WRITE(6,'(/)')
      PRINT*,'=======================================================',
     +       '======'
      PRINT*,'       BEGIN CALCULATION OF O_K_GR: ROUTE II    '
      PRINT*,'=======================================================',
     +       '======'
      PRINT*,' C U R V A T U R E  OF SPACE-TIME AT STARS EQUATOR:'
      PRINT*,' exp(psi-nu)/r_0=',EPMN/R0R0
      PRINT*,' dh_0/dr (>)=',DH0OSE
      PRINT*,' dh_2/dr (>)=',DH2OSE,'  dv_2/dr (>)=',DV2OSE
      PRINT*,' dnue/dr=',DVN,'   dpsi/dr=',DVP
      PRINT*,' exp(-Phi/2)=',EPHINV,'  exp(-nu/2)=',EPHINV/SQRT(TN)
      OEQ =2.*AGD*ZJD / RDCBEE**3
      DOEQ=-6.*AGD*ZJD / RDCBEE**4
      PRINT*,' omega(R_eq)=',OEQ,'  d omega/d r=',DOEQ
      T1 =DVN/DVP
      T22=(DOGDCB*EPMN / (2.*DVP*R0R0))**2
      PRINT*,' T1=',T1,'  T22=',T22
      PRINT*,' sqrt(t1+t22)=',SQRT(WZ)
      P11=DVN/DVP
      P22=-SQRT(T22)
      PRINT*,' (dnue/dr)/(dpsi/dr) =',P11
      PRINT*,' (domega/dr exp(psi-nu)/(2 dpsi/dr)=',P22
      PRINT*,' v(eq)/c=',VSDM,' OG(K;GR)=',OGDKDM/(UF10*R0R0),' rad/s'
      WRITE(6,'(/)')


 5787 CONTINUE
      THT2(MLF,8) =VSREQT
      TFREQ(K4,11)=VSREQT
      THT2(MLF,10)=VSMODI
      THT2(MLF,7) =OGDKGR/(UF10*R0R0)
      TFREQ(K4,10)=OGDKGR/(UF10*R0R0)
      THT2(MLF,9) =OGDKDM/(UF10*R0R0)


!.................................................................
!   copy computed values of V_eq(i)  and  Om_K(i)  for further use
!   in the self-consistent iteration loop determining Om_K
      VSDMH =VSREQT
      OGDKDM=OGDKGR
!.................................................................


!__ SAVE
      TFREQ(K4,2)=THT1(MLF,15)
      TFREQ(K4,3)=EC
      TFREQ(K4,4)=EC*UF6
!   /5/M_TOTAL, /6/R_EQUATOR, /8/ECCENTRICITY
      TFREQ(K4,5)=THT1(MLF,11)
      TFREQ(K4,6)=THT1(MLF,16)
      TFREQ(K4,8)=THT1(MLF,14)
!    POLAR REDSHIFT (DEFINITON OF KD = FIP)
      A44=EXP(PHIRCP) * (1.+2.*(H0RCBP+H2RCBP))
      IF (A44 <= 0.) THEN
      WRITE (6,'(///,'' HT1000: POLAR REDSHIFT CANNONT BE '',
     +           ''CALCULATED;  SQRT < 0 !!'',/,'' A44='',E14.6,///)')
     +      A44
      A44=1.
      END IF
      TFREQ(K4,9)=1./SQRT(A44) - 1.
!    EQUATORIAL REDSHIFT (DEFINITION OF KD)
      A33=EXP(PHIRCE) * (1.+2.*(H0RCBE-H2RCBE/2.))
     +    - (RDCBEE*ODCBEE)**2 * (1.-(V2RCBE-H2RCBE))
      IF (A33 <= 0.) THEN
      WRITE (6,'(///,'' HT1000: EQUATORIAL REDSHIFT CANNONT BE '',
     +           ''CALCULATED;  SQRT < 0 !!'',/,'' A33='',E14.6,///)')
     +      A33
      A33=1.
      END IF
      TFREQ(K4,13)=1./SQRT(A33) - 1.
!    INJECTION ENERGY
      TFREQ(K4,1)=1. / (TFREQ(K4,9)+1.)**2
!    RATIO OF ROTATIONAL TO GRAVITATIONAL ENERGY: T/W
      TKIN=ZJD*OGLD*A/2.
      TFREQ(K4,7)=1. / (1.+(THT2(MLF,4)-THT1(MLF,10))/TKIN)


      IF (OGDKGR > OGLD .AND. OGDKDM > OGLD) GO TO 4007

! ----------------------------------------------------------------------
! end of option loop: (1) e_c              ==>  self-cons. det. of Om_K,
!                     (2) M_rot and Omega  ==>  e_c,
!                     (3) e_c and Omega    ==>  M_rot
 4006 CONTINUE
! ----------------------------------------------------------------------


 4007 CONTINUE



!  COMPUTE NEUTRON STAR FREQUENCY  OMEGA ( I+1 )
      IF (ISELFC == 1) THEN
      OGDIT(2)=OGDKDM
      ELSE
      OGDIT(2)=OGDIT(1) + ALN * (OGDKDM-OGDIT(1))
      END IF

!  TEST OF CONVERGENCY
      DIFFSC=(OGDIT(2)-OGDIT(1)) / OGDIT(1)
      DIFFSC=ABS(DIFFSC)
      S1=OGDIT(2)
      S2=OGDR(IREM) + 2.*AGD*ZJD/FST251(IREM)**3
      DIFFOG=ABS(S1-S2)
      IF (I5651 == 1 .AND. ISCMX7 /= 1) THEN
      PRINT*,' '
      PRINT*,' ---------  check convergency ---------'
      PRINT*,' DIFFSC=',DIFFSC
      PRINT*,' DIFFI0=',DIFFI0
      PRINT*,' ISELFC=',ISELFC
      PRINT*,' DIFFOG=',DIFFOG
      PRINT*,' --------------------------------------'
      PRINT*,' '
      END IF

      IF (DIFFSC <= ETASC7) THEN
!     SELF-CONSISTENT CALCULATION SUCCESSFUL
      THT2(MLF,13)=DIFFSC
      GO TO 5761
      ELSE IF (ISELFC == ISCMX7) THEN
!     ITERATION LIMIT ENCOUNTERED
      THT2(MLF,13)=DIFFSC
      END IF

!__ COPY
      OGDIT(1) = OGDIT(2)
      IF (I5651 == 1  .AND.  ISCMX7 /= 1)
     +   RLOMG(K4OUT)=OGDIT(1)/(R0R0*UF10)
      IF (OGDIT(1) <= 0.) THEN
      WRITE (6,'(///,'' OMEGA(STAR,SELF-CONSISTENT) < 0'',/,
     +           '' STOP CALCULATION FOR THIS E(C) VALUE'',/,
     +           '' OMEGA(STAR,SELF-CONSISTENT)='',E14.7,/,
     +           '' ISELFC='',I3,///)') OGDIT(1),ISELFC
      GO TO 8962
      END IF


! -------------------------------------------------------------------
 5661 CONTINUE
!                 END OF  S E L F - CONSISTENT CALCULATION
! -------------------------------------------------------------------




 5761 CONTINUE
!__ WRITE RESULTS TO OUTPUT FILE

      OPEN (UNIT=10, FILE='frequency.dat', STATUS='unknown')

      DO 419 I=1,K4OUT
  419 WRITE (10,519) I,TFREQ(I,1),TFREQ(I,2),TFREQ(I,3),TFREQ(I,4),
     +TFREQ(I,5),TFREQ(I,6)
  519 FORMAT (1X,I3,2X,F6.3,2X,F5.2,2X,F7.2,2X,E12.6,2X,F6.2,2X,F5.2)
      WRITE (10,'(//)')
      DO 420 I=1,K4OUT
  420 WRITE (10,520) I,TFREQ(I,7),TFREQ(I,8),TFREQ(I,9),TFREQ(I,10),
     +TFREQ(I,11),TFREQ(I,12),TFREQ(I,13),THT1(MLF,6)
  520 FORMAT (1X,I3,2X,F6.3,2X,F5.2,2X,F7.4,2X,F8.2,2X,F4.2,2X,F8.2,
     +        2X,F7.4,2X,F8.2)

      CLOSE ( 10 )

      IF (IOP30 /= 0) THEN

      OPEN (UNIT=30,FILE='plotting.top',STATUS='unknown')

      DO 2030 I=1,IREM,10
      RQ=FST251(I)/FST251(IREM)
      WQ=OGDR(I)/OGLTD
 2030 WRITE (30,2031) RQ,WQ
      WRITE (30,2037)
 2037 FORMAT (' --------------------------------------------')
      DO 2032 I=1,IREM,10
      RQ=FST251(I)/FST251(IREM)
      EOO=OGDR(I)/(EXP(RNUED(I)/2.)*OGLTD)
 2032 WRITE (30,2033) RQ,EOO
 2031 FORMAT (E14.6,3X,E14.6)
 2033 FORMAT (E14.6,3X,E14.6)
      WRITE (30,2037)
      DO 2034 I=1,IREM,10
      RQ=FST251(I)/FST251(IREM)
 2034 WRITE (30,2035) RQ,Z0AR(I)
 2035 FORMAT (E14.6,3X,E14.6)
      WRITE (30,2037)
      DO 2038 I=1,IREM,10
      RQ=FST251(I)/FST251(IREM)
 2038 WRITE (30,2035) RQ,Z2AR(I)
      WRITE (30,2037)
      DO 2039 I=1,IREM,10
      RQ=FST251(I)/FST251(IREM)
      ECCENT= -3.*(V2RIN(I)-H2RIN(I)+Z2H/FST251(I))
      IF (ECCENT < 0.) ECCENT=0.
!     ECCENTRICITY ACCORDING TO  [HT68]
      IF (ECHT68 > 0.) THEN
      QECC=SQRT(ECCENT)/ECHT68
      ELSE
      QECC=-1.
      END IF
 2039 WRITE (30,2035) RQ,QECC

      CLOSE ( 30 )

      END IF

!   PLOT ENERGY DENSITY VERSUS PRESSURE (DIMENSIONLESS)
      IOP33=0
      IF (IOP33 /= 0) THEN

      OPEN (UNIT=31,FILE='PLOTN.EP',STATUS='unknown')

      DO 3037 I=1,IREM
      Y=FENG1(I)
      X=FOUT(I)
      IF (X <= 0. .OR. Y <= 0.) GO TO 3037
      Y=ALOG10(Y)
      X=ALOG10(X)
      WRITE (31,3038) X,Y
 3037 CONTINUE
 3038 FORMAT (E14.6,3X,E14.6)

      CLOSE ( 31 )

      END IF



      if (inot /= 0) then
!__ COMPUTE TOTAL BARYON NUMBER OF THE NON-ROTATING AS WELL AS THE ROTATING
!   NEUTRON STAR
      PRINT*,' '
      PRINT*,' --> COMPUTING THE STARS TOTAL BARYON NUMBER'

!                                                       ****** ******
      CALL BYNUMB(A,AR0,DGRAV,ZJD,RDSP,C,NHU,IREM,DELMH,RLATOT,RLASPH,
!                 *****
     +            DELEB,EDDY,PRESS,EINTH,PDEN,IEP,MLF,IATOTY)
!     RETURN: RLATOT=LOG [A_TOTAL]
!             RLASPH=LOG [A_SPH]
!             DELEB= DELTA E_BINDING/M_SUN
      end if



!__ COMPUTE THE STAR'S CRUST MASS (i.e., MASS AT DENSITIES BELOW THE
!   NEUTRON DRIP DENSITY)
!   NOTE: ACCESS ONLY IF A ROTATING QUARK STAR CALCULATION IS
!         BEING PERFORMED
      IF (IDPREM > 0) THEN
      PRINT*,' '
      PRINT*,' --> COMPUTING THE STARS CRUST MASS (RHO < RHO_drip)  AND'
      PRINT*,'     THE CRUSTS MOMENT OF INERTIA (perform  r* - O*',
     +       '  integration)'
      PRINT*,'     Omega_input =', OGLTD/(R0R0*UF10),' 1/sec'

      CALL CRUST(RJI,NI,NJ,NIP1,NJP1,THT,IREM,IDPREM,NHU,C,AX1,AX2,RTH,
!                ** ****
     +           WT,WT28,MLF,EC)
!      RETURN: WT  =STAR'S CRUST MASS IN UNITS OF THE SOLAR MASS
!              WT28=SAME, BUT IN UNITS OF 10**28 gramm
      THT2(MLF,16) = WT
      THT2(MLF,17) = WT28


      PRINT*,' --> COMPUTE MOMENT OF INERTIA OF DEFORMED STAR (call',
     +       'ing  _IDEFST_)'
      CALL IDEFST(RJI,NI,NJ,NIP1,NJP1,THT,IREM,IDPREM,NHU,C,AX1,AX2,
     +            RTH,MLF,EC,OGLTD,2,R_90,e_90,R_00,e_00)
!      RETURN: COMMON/TABHT1/....,THT2(.,.)
      END IF




      IF (IVONR /= 0) THEN
!  choose IVONR=10
!  compute  I(r), save results for later use in output routine
      PRINT*,' '
      PRINT*,' --> COMPUTE I_s(r) and A_s(r); results are being tabula-'
      PRINT*,'     ted // written on output file in subroutine OUTRHT!'
      DO 1739 I=1,IREM
      RIHAT = FST251(I)**4 * RLJ(I) * DOGDR(I) / OGLD
      RIHAT = RIHAT / (6. * DGRAV * A)
! log_10 (I(r)/g cm2)
      YOUT(I) = ALOG10(RIHAT) + ALOG10(EC*UF6) + 5.*DLOG10(R0R0/EE13)
! radial distance from origin in meters
      ST(I)   = FST251(I) * R0R0 / EE15
 1739 CONTINUE
! compute number of particles contained in sphere with radius r, LOG A_s(r)
      IADD=11
      DO 2978 IA=5,IREM
      RLASR=0.
      IF (IA == IADD .OR. IA == IREM) THEN
      IADD=IADD+10
      IGAS=IA
      EASS='HT1000: A_s(r)                '
! note: AXBN1 was being computed in subroutine BYNUMB!
!                                                ***
      CALL IGRN(FST251,AXBN1,IGAS,X251,Y251,IGAS,ASR,C,
     +          NHU,10,1,AUX51,EASS)
      RLASR=3.*DLOG10(AR0) + ALOG10(VPI*ASR)
      END IF
      PDENT(IA)=RLASR
 2978 CONTINUE
      END IF


      ITST68=10
      IF (ITST68 /= 0) THEN
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                CHECK SOLUTION  OMEGA_BAR(R)
      SUM124=0
      DO 124 I=1,IREM
      X251(I)=FST251(I)**4 * RLJ(I) * DOGDR(I)
  124 Y251(I)=4. * FST251(I)**3 * OGDR(I)
      I124M1=IREM-1
      DO 125 I=1,I124M1
      FST(I)=(X251(I+1)-X251(I)) / (FST251(I+1)-FST251(I))
      AUX5(I)=(RLJ(I+1)-RLJ(I))  / (FST251(I+1)-FST251(I))
!      DJDR=(-RLJ(I)/2.) * ( DNUEDR(I) + (2.*DGRAV/FST251(I)) *
!     +     (DMSDR(I)-F251(I)/FST251(I)) /
!     +     (1. - 2.*F251(I)*DGRAV/FST251(I)) )
!      D00=AUX5(I)-DJDR
      AUX5(I)=AUX5(I)*Y251(I)
      SUM124=SUM124 + FST(I)+AUX5(I)
  125 CONTINUE
      SUM124=ABS(SUM124)/FLOAT(IREM)
      END IF
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


! 06/05/05: 
!            write frame draggine frequency omega(r) as a function of r
!            to output file omega_vs_r.dat  
      if (ioutse /= 0) then

      open(unit=92, file='omega_vs_r.dat', status='unknown')

      do 199 i=1,irem
      if (i == 1) write(92,'(/,'' Frame dragging frequency  omega(r)/'',
     +                         ''Omega_K '',/)')
!      stellar_freq=ogld/(UF10*r0r0)        ! Stellar frequency in 1/s
      stellar_freq=ogld                     ! Stellar frequency 
  199 write (92,6210) fst251(i)/fst251(irem), ogdr(i)/stellar_freq
 6210 format (3x,e12.6,3x,e12.6)

      close ( 92 )
                    end if
! end of module


 8962 RETURN

      END

!MOIRDS
!                                                     **** **** ****
      SUBROUTINE MOIRDS(RAD1,RAD2,EC,QH,AUX3,AUX4,NHU,WI0H,WI2H,WJ0H,
!                       **** *** *** **
     +                  WJ2H,WIH,WJH,WI,C,IREM,OMGLH,RSPH,SUMW1,ZJD)


!-----------------------------------------------------------------------
!
! PURPOSE:
!            COMPUTE THE MOMENT OF INERTIA OF A ROTATIONALLY DISTORTED
!            STAR
!
!
! INPUT:
!        RADII: RAD1, RAD2, []=1,
!        OMEGA(STAR)=OMGLH,
!        OMEGA_BAR(R),
!        QUADRUPOLE MOMENT OF THE STAR=QUADM
!        MONOPOLE AND QUADRUPOLE DISTURBATION FUNCTIONS,
!        EQUATION OF STATE  P(E),
!        ENERGY DENSITY  E(R), PRESSURE  P(R),
!        DERIVATIVES: D[P]/D[E], D[E]/D[R], D[P]/D[R]
!
!
!  TO BE CALCULATED:
!                    W_1(HAT;R): FUNCTION DESCRIBING THE PERTURBATIONS OF
!                    ORDER  OMEGA(NEUTRON STAR)**3
!
! RETURN:
!         MOMENT OF INERTIA OF THE ROTATIONALLY DISTORTED STAR (i.e., 
!         MOI OF THE MATERIAL BETWEEN TWO SURFACES OF CONSTANT DENSITY
!         E(R1=RAD1)  AND  E(R2=RAD2))
!         WI0H=I_0(HAT), WI2H=I_2(HAT)
!         WJ0H=J_0(HAT), WJ2H=J_2(HAT)
!         WIH=I_HAT,     WJH=J_HAT   (TOTAL CONTRIBUTIONS)
!         WI =LOG_10 ( I / g cm**2)
!
!-----------------------------------------------------------------------


      REAL * 8 DGRAV,A,AR0,EE18,GCM2KM,EEM18,EE14,EE34,EE44
      REAL * 8 RMSUN,GRAV,GRAVK,GRAVKM 
      real * 8 ee03,ee13,ee15
      
      DIMENSION AUX3(NHU),AUX4(NHU),C(4,NHU),DMX(3,4),DCOEFF(3,3)

      COMMON/ExpoMaxMin/EXMAX,EXMIN
      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      common/ConversionFactor30/ufe30
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/GRAVITATIONAL/GRAVK,GRAV,GRAVKM,GCM2KM
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/SOLARPROPERTIES/RMSUN,RSUN,RSUNKM
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/MISCINT2/J651,J650,IGIVMX,IVONR,KC20
      COMMON/IPLOT1/IRADOT,IOUT20,IOTOMG,IOUT33,IOUTES,IOTMP0,IOP30,IOTH
     +              V2,IOUTPE,IOUTSE
      COMMON/RADIALSCALINGFACTOR/R0R0
      COMMON/REDSEP/H0RCBE,H0RCBP,PHIRCP,PHIRCE,RDCBEP,RDCBEE,ODCBEE,
     +              H2RCBE,H2RCBP,V2RCBE,V2RCBP

      PARAMETER (m651=160501)
      COMMON/CNSTA1/F251(M651),FOUT(M651),FST(M651),FST251(M651),
     +              FENG1(M651),X251(M651),Y251(M651),F251PM(M651),
     +              FDRIP(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)
      COMMON/BN01/AXBN1(M651),AXBN2(M651),PDENEX(M651),EINTHX(M651),
     +            DEINDP(M651),DEHDRH(M651),PDENT(M651)
      COMMON/HT68A/DNUEDR(M651),RLJ(M651),DOGDR(M651),DJ2DR(M651)
      COMMON/HT68B/DEDP(M651),RNUED(M651),OGDR(M651),DLNJ(M651)
      COMMON/HT68C/RM0D(M651),P0SD(M651),DM0DR(M651),DP0SDR(M651)
      COMMON/HT68D/DMSDR(M651),DPHDRH(M651),P2SD(M651),RM2H(M651)
      COMMON/HT68E/H0OUT(M651),H0IN(M651),DH0IN(M651)
      COMMON/HT68F/DV2DR(M651),DH2DR(M651),V2RIN(M651),H2RIN(M651)
      COMMON/HT68G/Z0AR(M651),Z2AR(M651)
      COMMON/HT68H/DV2HOM(M651),DH2HOM(M651),V2HOM(M651),H2HOM(M651)
      COMMON/HT68I/OGAUX(M651),DOGAUX(M651),D2OGAX(M651)
      COMMON/MRDS1/W1H(M651),D1W1H(M651),D2W1H(M651),DH0R(M651),
     +             DH2R(M651),ARRIH(M651)
      COMMON/AXX17/AUXW11(M651),AUXW12(M651)
      COMMON/MRDS2/DRM2H(M651),D2OGDR(M651),D2JDR2(M651),D2J22(M651)


      DATA RMC/0./
      DATA FMKM/1.E18/,FMEVFM/5.61012/

      NAMELIST/NAME1/RMSUN,FOUT,FST251,J650,J651 




! ...................... AUXILIARY  FUNCTIONS ..........................


!   FUNCTIONS  C(1,R)  AND  C(2,R)
      FC1R(X,Y) = (4. + X*Y) / X
      FC2R(X,Y) = 4. * Y / X

!   FUNCTION  B(Y,Z)
      FUNB(Y,Z) =  (3./8.)
     +          * ( - 25./18. + 1./Y + 1./(2.*Y*Y) - 1./(9.*Y*Y*Y)
     +              + ALOG(1./Y) * (-1./Y-1./(Y*Y)+1./(3.*Y*Y*Y)) )
     +          + 3./(2.*Z*Z) + 1./(Z*Z*Z)-1./(2.*Z*Z*Z*Z)

!   D [B] / D [R]
      DFUNB(Y,DY,Z,DZ) =  (3./8.)
     +           * (DY*ALOG(1./Y)*(1.+2./Y-1./(Y*Y))/(Y*Y))
     +           - 3.*DZ*(1.+1./Z-2./(3.*Z*Z)) / (Z*Z*Z)

!   D2 [B] / D [R]2
      D2FUNB(Y,DY,DDY,Z,DZ,DDZ) = (3./8.)
     +           * (  (ALOG(1./Y)*(Y*DDY-2.*DY*DY)-DY*DY) *
     +                (1.+2./Y-1./(Y*Y)) / (Y*Y*Y)
     +              - 2.*DY*DY*ALOG(1./Y)*(Y-1.)/Y**5  )
     +           - 3.*(Z*DDZ-3.*DZ*DZ)*(1.+1./Z-2./(3.*Z*Z))/Z**4
     +           + 3.*DZ*DZ*(Z-4./3.)/Z**6


! ....................... END OF AUX FUNCTIONS .........................





!__ CONSTANTS
      AR0   =R0R0
      A     =EC*AR0**3/(RMSUN*ufe30)
      DGRAV =(GRAV/ufe30)*(RMSUN*ufe30)/(UF1*AR0)


!__ COMPUTE STARTING VALUES FOR
!                                -  W(1,HAT;R=0)
!                                -  D [W(1,HAT;R=0)] / D[R]
!                                -  D2 [W(1,HAT;R=0)] / D[R]2
      DO 99 IORIG=1,3
      DRNS=FST251(IORIG)
      DRNS2=DRNS*DRNS
      DMNS=F251(IORIG)
      RMGRS=DMNS*DGRAV/DRNS
      W1S=1.-2.*RMGRS

!   FUNCTIONS  C(1,R)  AND  C(2,R)  AT  R=R_s
      DLNJDR=DLNJ(IORIG)
      C1R=FC1R(DRNS,DLNJDR)
      C2R=FC2R(DRNS,DLNJDR)

!   FUNCTION  D(0,HAT)  AT  R=R_s
      RUCK=1. + 2.*DGRAV*(DMNS/DRNS - DMSDR(IORIG)) / W1S
      ECKK=DM0DR(IORIG) - RM0D(IORIG) * RUCK / DRNS
      XRB1=DGRAV*ECKK / (DRNS*W1S) + DH0IN(IORIG)
      XRB2=  2.*RM0D(IORIG)*DGRAV / (DRNS*W1S)
     +     + P0SD(IORIG) * (1. + DEDP(IORIG))
     +     + 2.*(OGDR(IORIG)*DRNS)**2 / (3.*EXP(RNUED(IORIG)))
      D0= - DOGDR(IORIG)*RLJ(IORIG)*XRB1
     +    + 4.*RLJ(IORIG)*DLNJ(IORIG)*OGDR(IORIG)*XRB2/DRNS
      DH0R(IORIG)=D0
!   FUNCTION  D(2,HAT)  AT  R=R_s
      Z1=  4.*(DV2DR(IORIG)-DH2DR(IORIG))
     +    +(RM2H(IORIG)/DRNS-DRM2H(IORIG)) / (DRNS*W1S)
      Z2=  2.*RM2H(IORIG)*DGRAV*(DMNS/DRNS - DMSDR(IORIG))
     +     / (DRNS2*W1S*W1S)
      Z3= -DH2DR(IORIG)
      XRB1=Z1+Z2+Z3
      XRB2=  2.*RM2H(IORIG) / (DRNS*W1S)
     +      + P2SD(IORIG) * (1. + DEDP(IORIG))
     +      - 2.*(OGDR(IORIG)*DRNS)**2 / (3.*EXP(RNUED(IORIG)))
      D2=   DOGDR(IORIG)*RLJ(IORIG)*XRB1
     +    + 4.*RLJ(IORIG)*DLNJ(IORIG)*OGDR(IORIG)*XRB2/DRNS
      DH2R(IORIG)=D2
      DMX(IORIG,1) = 1.
      DMX(IORIG,2) = C1R
      DMX(IORIG,3) = C2R
      DMX(IORIG,4) = (D0 - D2/5.) / RLJ(IORIG)

   99 CONTINUE

!   COMPUTE DETERMIANTS OF THE COEFFICIENT METRICES
      DO 77 I=1,4
         DO 88 J=1,3
            DO 88 K=1,3
            DCOEFF(J,K) = DMX(J,K)
   88    CONTINUE
      IF (I /= 1) THEN
      DO 66 L=1,3
   66 DCOEFF(L,I-1) = DMX(L,4)
      END IF
!  COMPUTE DETERMINANT OF MATRIX  DCOEFF
!                            ******
      CALL DETMAX(DCOEFF,3,3,DETERM)

      IF (I == 1) DETD =DETERM
      IF (I == 2) DETD1=DETERM
      IF (I == 3) DETD2=DETERM
      IF (I == 4) DETD3=DETERM

      IF (DETD == 0.) THEN
      WRITE (6,'(//,'' DETERMIANT < 0,  MATRIX INVERSION FAILED'',
     +           '' IN  _MOIRDS_'',/,'' DETD='',E14.8,/,''DETD1='',
     +           E14.8,/,'' DETD2='',E14.8,/,'' DETD3='',E14.8)')
     +      DETD,DETD1,DETD2,DETD3
      STOP 'DET < 0, _MOIRDS_'
      END IF

   77 CONTINUE

!................................
!   COMPUTE STARTING VALUES
      D2W1H(1)=DETD1/DETD
      D1W1H(1)=DETD2/DETD
      W1H(1)  =DETD3/DETD
!................................



! *********************************************************************
! START SOLVING THE SECOND-ORDER DIFFERENTIAL EQUATION OF   W(1,HAT;R)
! IN ORDER TO GET THE  I N T E R I O R  SOLUTION
! *********************************************************************
      DO 905 I=2,IREM
      IIM1=I-1
      DRN=FST251(IIM1)
      DRN2=DRN*DRN
      DMN=F251(IIM1)
      RMGR=DMN*DGRAV/DRN
      W1=1.-2.*RMGR

      IF (W1 < 0.) THEN
         WRITE(6,'(///,'' _MOIRDS_ W1 < 0, SQRT COMPLEX!!'')')
         WRITE(6,'(//,'' W1='',E14.6,///)') W1
         STOP
      END IF

      DLNJDR=DLNJ(IIM1)
      C1R=FC1R(DRN,DLNJDR)
      C2R=FC2R(DRN,DLNJDR)

!   FUNCTION  D(0,HAT)
      RUCK=1. + 2.*DGRAV*(DMN/DRN - DMSDR(IIM1)) / W1
      ECKK=DM0DR(IIM1) - RM0D(IIM1) * RUCK / DRN
      XRB1=DGRAV*ECKK / (DRN*W1) + DH0IN(IIM1)
      XRB2=  2.*RM0D(IIM1)*DGRAV / (DRN*W1)
     +     + P0SD(IIM1) * (1. + DEDP(IIM1))
     +     + 2.*(OGDR(IIM1)*DRN)**2 / (3.*EXP(RNUED(IIM1)))
      D0= - DOGDR(IIM1)*RLJ(IIM1)*XRB1
     +    + 4.*RLJ(IIM1)*DLNJ(IIM1)*OGDR(IIM1)*XRB2/DRN
      DH0R(IIM1)=D0

!   FUNCTION  D(2,HAT)
      Z1=  4.*(DV2DR(IIM1)-DH2DR(IIM1))
     +    +(RM2H(IIM1)/DRN-DRM2H(IIM1)) / (DRN*W1)
      Z2=2.*RM2H(IIM1)*DGRAV*(DMN/DRN - DMSDR(IIM1))
     +   / (DRN2*W1*W1)
      Z3=-DH2DR(IIM1)
      XRB1=Z1+Z2+Z3
      XRB2=  2.*RM2H(IIM1) / (DRN*W1)
     +     + P2SD(IIM1) * (1. + DEDP(IIM1))
     +     - 2.*(OGDR(IIM1)*DRN)**2 / (3.*EXP(RNUED(IIM1)))
      D2=   DOGDR(IIM1)*RLJ(IIM1)*XRB1
     +    + 4.*RLJ(IIM1)*DLNJ(IIM1)*OGDR(IIM1)*XRB2/DRN
      DH2R(IIM1)=D2

      HRK=FST251(I)-FST251(IIM1)
!...............................................................
!   ETA(HAT;R)  = D [W(1,HAT;R)] / D [R]
      SOURCE    = (D0 - D2/5.) / RLJ(IIM1)
      DETDX     = - C1R*D1W1H(IIM1) - C2R*W1H(IIM1) + SOURCE
!-- D [W(1,HAT;R + H)] / D [R]
      D1W1H(I ) = D1W1H(IIM1) + HRK * DETDX

!-- W(1,HAT;R + H)
      DW1DX     = D1W1H(IIM1)
      W1H(I)    = W1H(IIM1) + HRK * DW1DX

!-- D2 W(1,HAT;R + H) / D[R]2
      D2W1H(I)  = - C1R*D1W1H(I) - C2R*W1H(I) + SOURCE
!...............................................................

  905 CONTINUE




      IF (IOTOMG /= 0) THEN
      DO 908 I=1,IREM
      IF (I == 1)
     +WRITE (6,'(///,'' _MOIRDS_ W(1,HAT;R),  D[W1]/D[R],  '',
     +                ''D2[W1]/D[R]2'',/)')
      KK908=I
      IF (I > 50) KK908=(-1)**I
      IF (KK908 < 0) GO TO 908
      PRTOG=W1H(I)/(UF10*R0R0)
      WRITE (6,'('' I='',I4,''  W1='',E12.6,'' 1/S; DW1='',E12.6,
     +           '' D2W1='',E12.6)') I,PRTOG,D1W1H(I),D2W1H(I)
 908  CONTINUE
      WRITE (6,'(/,'' IREM='',I6,/)') IREM
      PRTOG=W1H(IREM)/(UF10*R0R0)
      WRITE (6,'('' IREM='',I6,''  W1='',E12.6,'' 1/S; DW1='',E12.6,
     +       '' D2W1='',E12.6)') IREM,PRTOG,D1W1H(IREM),D2W1H(IREM)
      END IF



! *********************************************************************
! JOIN   I N T E R I O R   SOLUTION  W(<)=W(1,HAT;R)  SMOOTHLY WITH
! THE  E X T E R I O R   SOLUTION  W(>)=W(1,HAT;R)  AT THE STELLAR
! SURFACE, i.e., AT  R=R_s.
! *********************************************************************
      DRNS=FST251(IREM)
      DRNS2=DRNS*DRNS
      DMNS=F251(IREM)
      RMGRS=DMNS*DGRAV/DRNS
      W1S=1.-2.*RMGRS

!   COMPUTE  F(R)  AT  R=R_s
      X=DRNS/(DGRAV*DMNS)
      Y=(X+1.)/(X-1.)
      Z=X+1.
      QUADM=QH*DRNS**2*DMNS - (ZJD*A)**2/DMNS
      FCB=FUNB(Y,Z)
      FCF = - (4./5.) * (ZJD*A)**3 * DGRAV**2 / (DMNS*DRNS**6)
     +      + (48./35.) * (ZJD*DGRAV*A)**3 / DRNS**7
     +      - 3. * ZJD*A * QUADM * FCB /
     +        (DMNS**3 * DGRAV * DRNS**3)

!   COMPUTE  D [F] / D [R]  AT  R=R_s
      RAD=1.-1./(DMNS*DGRAV/DRNS)
      DY=-2. / (DMNS*DGRAV * RAD*RAD)
      DZ=1./(DMNS*DGRAV)
      DBDR=DFUNB(Y,DY,Z,DZ)
      DFCF =   (24./5.) * (ZJD*A)**3 * DGRAV**2 / (DMNS * DRNS**7)
     +       - (48./5.) * (ZJD*DGRAV*A)**3 / DRNS**8
     +       + 9. * ZJD*A * QUADM * FCB  / (DMNS**3*DGRAV*DRNS**4)
     +       - 3. * ZJD*A * QUADM * DBDR / (DMNS**3*DGRAV*DRNS**3)

!   COMPUTE  D2 [F] / D [R]2   AT  R=R_s
      DDY=- (2./(DMNS*DGRAV))**2 / RAD**3
      DDZ=0.
      D2BDR =D2FUNB(Y,DY,DDY,Z,DZ,DDZ)
      BDBD2B=12.*FCB/DRNS - 6.*DBDR + DRNS*D2BDR
      D2FCF =- (168./5.) * (ZJD*A)**3 * DGRAV**2 / (DMNS * DRNS**8)
     +       + (348./5.) * (ZJD*DGRAV*A)**3 / DRNS**9
     +       - 3. * ZJD*A * QUADM * BDBD2B / (DMNS**3*DGRAV*DRNS**4)




!>> COMPUTE  DELTA(J)  (CONSTANT OF INTEGRATION)  AT  R=R_s
!.............................................................
      DELJH = (FCF - W1H(IREM)) *DRNS**3 / (2 * A * DGRAV)
!.............................................................



!>> VALUE OF  D [W(1,HAT;R)] / D [R]  AT  R=R_s,
!   D [W(<)] / D [R] = D [W(>)] / D [R]:
      DW1EXT= 6.*A*DGRAV * DELJH / DRNS**4 + DFCF
      DIFF=ABS(D1W1H(IREM)-DW1EXT)


!>> VALUE OF  D2 [W(1,HAT;R)] / D [R]2  AT  R=R_s,
!   D2 [W(<)] / D [R]2 = D2 [W(>)] / D [R]2:
      D2W1EX= -24.*A*DGRAV * DELJH / DRNS**5 + D2FCF
      DIFF=ABS(D2W1H(IREM)-D2W1EX)




! *****************************************************************
! COMPUTE THE STAR'S MOMENT OF INERTIA:
!       I(R,OMEGA) = I_0(R) + I_2(R) * OMEGA(STAR)**2
! *****************************************************************
      DO 122 I=1,IREM
      DRN=FST251(I)
      DRN2=DRN*DRN
      DRN3=DRN2*DRN
      DRN4=DRN3*DRN

!   FIRST STEP: COMPUTE  I_0(R)
      ZID0=DRN4 * RLJ(I) * DOGDR(I) / (6.*A*DGRAV * OMGLH)


!   SECOND STEP: COMPUTE  I_2(R)
      DMN=F251(I)
      RMGR=DMN*DGRAV/DRN
      W1=1.-2.*RMGR
      RM0GR=RM0D(I)*DGRAV / (DRN*W1)
      RK2=V2RIN(I)-H2RIN(I)
      ZID2=( DRN4*RLJ(I) * (D1W1H(I) + DOGDR(I)*( H0IN(I)+RM0GR
     +      +(4.*RK2-H2RIN(I)-RM2H(I)/(DRN*W1)) / 5. ) )
     +      +4.*DRN3*RLJ(I)*DLNJ(I)*DOGDR(I)*(Z0AR(I)-Z2AR(I)/5.) )
     +      / (-6. * OMGLH**3)


!   SAVE  I_0(HAT)  AND  I_2(HAT)  IN AUXILIARY ARRAYS
      AUXW11(I)=ZID0
      AUXW12(I)=ZID2

!   ADD MONOPOLE AND QUADRUPOLE CONTRIBUTIONS OF THE MOMENT OF INERTIA
      ARRIH(I)=ZID0 + ZID2*OMGLH**2 / (A*DGRAV)

  122 CONTINUE



!>> COMPUTE THE MOMENT OF INERTIA OF A REGION BETWEEN SURFACES OF
!   CONSTANT DENSITY  E(R1)  AND  E(R2):
!   I(R1,R2;OMEGA)=I(R1,OMEGA)-I(R2,OMEGA)
      IF (RAD1 == 0. .AND. RAD2 == FST251(IREM)) THEN
      WI0H = AUXW11(IREM)
      WI2H = AUXW12(IREM)
      WJ0H = WI0H * OMGLH
      WJ2H = WI2H * OMGLH**3 / (A*DGRAV)
      WIH  = ARRIH(IREM)
      WJH  = WJ0H + WJ2H

      ELSE

      WRITE (6,'(///,'' R(1)#0 OR R(2)#R(SPH)  N O T  ALLOWED'')')
      STOP 'ERROR STOP IN _MOIRDS_'

      END IF



!   RETURN: LOG_10 ( I / (g cm**2) )
      WI= ALOG10(WJH)+ALOG10(EC*UF6)+5.*DLOG10(R0R0/EE13)-ALOG10(OMGLH)
      PRINT*,' --> SUB  MOIRDS:  log I_total(new)=',WI


      ITSTTA=0
      IF (ITSTTA /= 0) THEN
!**********************************************************************
      WRITE (6,'(//)')
      WRITE (6,'('' THE VALUES OF THIS TABLE MUST AGREE WITH THOSE'',
     +'' GIVEN IN TABLE I IN [Har73]'',/,'' (Astrophys. and Space'',
     +'' Sci. 24 (1973) 385, p. 396 there)'',/)')
      WRITE (6,50)
   50 FORMAT (' ---------------------------------------------------')
      YEC=EC*UF6
      YRS=RSPH/EE03
      YMS=F251(IREM)
      RI0=ALOG10(WI0H*YEC) + 5.*DLOG10(R0R0/EE13) - 44.
     +    + ALOG10(7.421943)
      YI0=10.**RI0
      RI2=ALOG10(WI2H)     + 5.*DLOG10(R0R0*EEM18)

      YI2=10.**RI2
      YYA=YI0**2 / (4.*YI2)
      WRITE (6,51)
   51 FORMAT ('  e_c      R_s      M_s      I_0       I_2       A')
      WRITE(6,52)
   52 FORMAT (' g/cm**3    km     M_sun     km**3     km**5     km')
      WRITE (6,50)
      WRITE (6,53) YEC,YRS,YMS,YI0,YI2,YYA
   53 FORMAT (1X,E10.3,2X,F5.2,2X,F6.4,2X,F5.2,2X,E10.3,2X,E10.3)
      WRITE (6,'(///)')
      END IF
!**********************************************************************




      ITST68=10
      IF (ITST68 /= 0) THEN
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                CHECK SOLUTION  W(1,HAT;R)
      SUMW1=0
      DO 124 I=1,IREM
      DRN=FST251(I)
      DRN2=DRN*DRN
      DMN=F251(I)
      RMGR=DMN*DGRAV/DRN
      W1=1.-2.*RMGR
      
!   FUNCTIONS  C(1) AND C(2)
      DLNJDR=DLNJ(I)
      C1R=FC1R(DRN,DLNJDR)
      C2R=FC2R(DRN,DLNJDR)
      AAA=D2W1H(I) + C1R*D1W1H(I) + C2R*W1H(I)
      BBB=(DH0R(I)-DH2R(I)/5.) / RLJ(I)
      DIZ=ABS(AAA-BBB)
      SUMW1=SUMW1 + DIZ
  124 CONTINUE
      SUMW1=SUMW1/FLOAT(IREM)
      END IF
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX





      RETURN

      END


!DETMAX                         ******
      SUBROUTINE DETMAX(D,ID,JD,DETERM)
!
! ----------------------------------------------------------------------
!
!__ PURPOSE:
!             TO CALCULATE THE DETERMIANT OF AN  ID X JD  MATRIX
!             WHOSE ELEMENTS ARE CONTAINED IN   D(ID,JD).
!
!   RETURN:
!             DETERM  CONTAINS THE VALUE OF THE DETERMINANT
!
! ----------------------------------------------------------------------


      DIMENSION D(ID,JD)

      IERR=0
      IF (ID /= JD .OR. ID /= 3) THEN
      WRITE (6,'(//,'' INPUT MATRIX  D  IN  _DETMAX_  IS NOT A'',
     +           '' SQUARE MATRIX'',/,''ID='',I3,'', JD='',I3)')
     +      ID,JD
      IERR=-1
      END IF
      IF (IERR < 0) STOP 'FATAL ERROR DETECTION IN  _DETMAX_'


!  COMPUTE DETERMINANT
      T1 = D(1,1) * ( D(2,2)*D(3,3) - D(2,3)*D(3,2) )
      T2 = D(2,1) * ( D(1,2)*D(3,3) - D(1,3)*D(3,2) )
      T3 = D(3,1) * ( D(1,2)*D(2,3) - D(1,3)*D(2,2) )
      DETERM = T1 - T2 + T3


      RETURN

      END


!heating_project
      subroutine heating_project(irem,ec,anuefip,k_qmc,k_mix,i4006a,RMoI
     +                                                                 )

      real n_i, M_i, Mrot_i, nueeq_i, nueeqFIP_i

      real * 8 rmsun,grav,gravk,gravkm,ee18,gcm2km,eem18
      real * 8 ee14,ee34,ee44,dgrav
      real * 8 ee03,ee13,ee15

      parameter (m651=160501)
      common/cnsta1/f251(m651),fout(m651),fst(m651),fst251(m651),feng1
     +              (m651),x251(m651),y251(m651),f251pm(m651),fdrip(m6
     +              51)
      common/iplot1/iradot,iout20,iotomg,iout33,ioutes,iotmp0,iop30,ioth
     +              v2,ioutpe,ioutse
      common/cnsta4/aux71(1,m651),aux7(m651),aux8(m651),aux9(m651)
      common/ht68a/DNUEDR(m651),RLJ(m651),DOGDR(m651),DJ2DR(m651)
      common/ht68b/dedp(m651),rnued(m651),ogdr(m651),dlnj(m651)
      common/ht68c/RM0D(m651),P0SD(m651),DM0DR(m651),DP0SDR(m651)
      common/ht68e/H0OUT(M651),H0IN(M651),DH0IN(M651)
      common/ht68f/DV2DR(m651),DH2DR(m651),V2RIN(m651),H2RIN(m651)
      common/ht68g/Z0AR(m651),Z2AR(m651)
      common/bn01/axbn1(m651),axbn2(m651),pdenex(m651),einthx(m651),
     +            deindp(m651),dehdrh(m651),pdent(m651)

      common/heating/Req_i(m651),Rpo_i(m651),Rsph_i(m651),Asph_i(m651),
     +               Arot_i(m651),RMTOTH_i(m651)
      COMMON/massfractions/RMquark(m651),RMMixed(m651),RMNucl(m651)

      parameter (n27=50,nptb3=10,nttby2=15)
      common/tabbbb/tabb(n27,nptb3)
      common/tbynb1/ttby(n27,nttby2)

      common/Gravitational/GRAVK,GRAV,GRAVKM,GCM2KM
      common/PowersTen/ee03,ee15,ee18,eem18,ee13,ee14,ee34,ee44
      common/RadialScalingFactor/r0r0
      common/ConversionFactor30/ufe30
      common/ConversionFactors/uf1,uf12,uf13,uf4,uf6,uf10
      common/SolarProperties/RMSUN,RSUN,RSUNKM
      common/ctype/typeos,teoscb

      character * 50 typeos
      character * 45 teoscb

      data spol/2.99792458e+05/  ! speed of light in km/s

      open (unit=19, file='heatingProject.dat', status='unknown')

      print *,' --> computing input for latent heat project'

      DGRAV =(GRAV/ufe30)*(RMSUN*ufe30)/(UF1*R0R0)

!       write (19,'('' jobname='',a50)') typeos
!       write (19,'('' Omega='',f8.2,'' 1/sec'')') tabb(i4006a,1)
      write (19,'(f8.2)') tabb(i4006a,1)
!       write (19,'('' log_10 A_rot='',f9.5)') ttby(i4006a,4)
!      write (19,4) k_qmc, k_mix
! 4    format (2x,i4,2x,i4)
      write (19,'(f8.4)') RMoI
      write (19,'(i5)') irem

!      write (19,5) anuefip  ! Metric nue(r=R_eq) as defined in [FIP86]
! 5    format (2x,'nu(R_eq)=',1x,f8.5,4x,'(ds^2 = - exp(2*nue)) dt^2 +',
!     + ' ... )')


      RMNucl(irem)   = RMNucl(irem-1)
      Arot_i(irem)   = Arot_i(irem-1)
      rmtoth_i(irem) = rmtoth_i(irem-1)

      p16            = r0r0/ee18

      do j=1,irem

      e_i    = feng1(j)*ec ! energy density(r) in MeV/fm^3
      p_i    = fout(j)*ec  ! pressure(r) in MeV/fm^3
      n_i    = pdenex(j)   ! baryon number density(r) in 1/fm^3
      M_i    = f251(j)     ! non-rotating stellar mass(r) in M_sol
      Mrot_i = rmtoth_i(j) ! rotating stellar mass(r) in M_sol
      Phi_i  = rnued(j)    ! metric Phi(r) according to Hartle's notation:
!                            ds^2 = - e^Phi dt^2 + ...]

!__  compute and tabulate star's deformation in equatorial and polar 
!    direction
      if(j == 1) print *,'     1. polar and equatorial deformation'
      req_i = FST251(j) + Z0AR(j) + ( Z2AR(j) + FST251(j) * (V2RIN(j)
     +                               - H2RIN(j)) ) * (-0.5)
      req_i = req_i*p16  ! equatorial direction in km
      rpo_i = FST251(j) + Z0AR(j) + ( Z2AR(j) + FST251(j) * (V2RIN(j)
     +                               - H2RIN(j)) ) 
      rpo_i = rpo_i*p16  ! polar direction in km

!__  compute and tabulate metric function nue(r), as defined in [FIP86],
!    in equatorial direction
!      W1EQ = 1. - 2.* RMTOTH_i(j)*DGRAV / (req_i(j)*EE18/r0r0) 
!      Phieq_i = ALOG(W1EQ)
      if(j == 1) print *,'     2. metric functions nu and Psi'
      Phieq_i     = Phi_i 
      TN          = 1. + 2. * (H0IN(j) - H2RIN(j)/2.)
      nueeq_i     = Phieq_i + ALOG(TN)   ! Hartle's convention
      nueeqFIP_i  = nueeq_i / 2.         ! Change to [FIP86], [Weber's book]
      PhiFIP_i    = Phi_i / 2.           ! Change to [FIP86], [Weber's book]

      term1  = fst251(j)*p16
      term12 = term1*term1
      term2  = 1. - v2rin(j) + h2rin(j)
      expPsi = term12 * term2      ! e^Psi  (Hartle's convention)
      expnu  = exp(nueeq_i)        ! e^nu   (Hartle's convention)
      term3  = ogdr(j)*ogdr(j) / (uf10*r0r0*spol)
      if(j == 1) print *,'     3. metric tensor g_00'
      g_00_rot = - expnu + expPsi * term3  ! g_00 = -e^nu + e^Psi omega^2
      g_00_sph = - exp(Phi_i)              ! g_00 = -e^Phi


!__ compute and tabulate masses of quark matter core(r), mixed phase(r), 
!   and hadronic liquid(r) [requires input value for ioutse /= 0]
      if(j == 1) print *,'     4. fraction of matter in various phases'
      delm       = rm0d(j) + fst251(j)**5 * dogdr(j)**2 / (36.*dgrav)
      if (j <= k_qmc) then
      RMquark(j) = RMquark(j) + delm
      else if (j <= k_mix) then
      RMMixed(j) = RMMixed(j) + delm
      else 
      RMNucl(j)  = RMNucl(j)  + delm
      end if

      
      if (j == 1) write (19,11) 
 11   format ('   # ',1x,'r_eq(A)',3x,'r_p(A)',2x, 'log_10 A',6x, 
     + 'e(A)',10x,'p(A)',9x,'n(A)',10x,'M(A)',5x,'M_quark(A)',4x,
     + 'M_mixed(A)',3x,'M_nuclear(A)',4x,'nu_eq',6x,'Phi',4x,'g_00_sph',
     + 2x,'g_00_rot') 
      write (19,10) j,req_i(j),rpo_i(j),Arot_i(j),e_i,p_i,n_i,Mrot_i,
     +              RMquark(j),RMMixed(j),RMNucl(j),nueeqFIP_i,PhiFIP_i,
     +              g_00_sph,g_00_rot
 10   format (i4,2x,f7.4,2x,f7.4,2x,f8.5,2x,e12.6,2x,e12.6,2x,e12.6,2x,
     +        e9.4,2x,e12.6,2x,e12.6,2x,e12.6,2x,f8.5,2x,f8.5,2x,f8.5,
     +        2x,f8.5)
      end do

      if (ioutse /= 0) total_mass = RMquark(k_qmc)+RMMixed(k_mix)
     +                              +RMNucl(irem)
!      write (19,200) total_mass
! 200  format (2x,'Check stars total mass:',1x,f7.5)

      close (19)

      return
              end


!BYNUMB                                                       ******
      SUBROUTINE BYNUMB(A,AR0,DGRAV,ZJD,RDSP,C,NHU,IREM,DELMH,RLATOT,
!                       ****** *****
     +                  RLASPH,DELEB,EDDY,PRESS,EINTH,PDEN,IEP,MLF,IA
     +                  TOTY)


! ----------------------------------------------------------------------
!
!__ PURPOSE:
!             1. TO CALCULATE THE TOTAL BARYON NUMBER OF THE DEFORMED
!             (ROTATING) NEUTRON STAR
!             2. Compute/save data needed for latent heat studies;
!                the changes are indicated by 'cc----+/-'; (06/2007)
!
!
!   RETURN: RLATOT = LOG ( TOTAL BARYON NUMBER OF ROTATING STAR)
!           RLASPH = LOG ( TOTAL BARYON NUMBER OF SPHERICAL STAR)
!           DELEB  = DELTA E(BINDING) / M_SUN
!
! ----------------------------------------------------------------------


      REAL * 8 DGRAV,A,AR0,RMSUN,RMNTH,DA,ATOT,ASPHQ,EE18,EEM18
!---+
      real * 8 ATOTrad,DArad,ASPHQrad
!---+
      REAL * 8 EE14,EE34,EE44
!---+
      REAL * 8 GRAV,GRAVK,GRAVKM,GCM2KM
!---+
      real * 8 ee03,ee13,ee15
!---+      
      REAL C(4,NHU)
      REAL EDDY(NHU),PRESS(NHU),EINTH(NHU),PDEN(NHU)


      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      common/ConversionFactor30/ufe30
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
!c---+
      COMMON/GRAVITATIONAL/GRAVK,GRAV,GRAVKM,GCM2KM
!c----
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/SOLARPROPERTIES/RMSUN,RSUN,RSUNKM
!c---+
      COMMON/RADIALSCALINGFACTOR/R0R0
!c----

      PARAMETER (m651=160501)
      COMMON/CNSTA1/F251(M651),FOUT(M651),FST(M651),FST251(M651),
     +              FENG1(M651),X251(M651),Y251(M651),F251PM(M651),
     +              FDRIP(M651)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA3/AUX61(1,M651),AUX6(M651),DPDRAX(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)
      COMMON/BN01/AXBN1(M651),AXBN2(M651),PDENEX(M651),EINTHX(M651),
     +            DEINDP(M651),DEHDRH(M651),PDENT(M651)
      COMMON/HT68A/DNUEDR(M651),RLJ(M651),DOGDR(M651),DJ2DR(M651)
      COMMON/HT68B/DEDP(M651),RNUED(M651),OGDR(M651),DLNJ(M651)
      COMMON/HT68C/RM0D(M651),P0SD(M651),DM0DR(M651),DP0SDR(M651)
!c---+
      common/heating/Req_i(m651),Rpo_i(m651),Rsph_i(m651),Asph_i(m651),
     +               Arot_i(m651),RMTOTH_i(m651)
!c----

      PARAMETER (N27=50,NTTBY2=15)
      COMMON/TBYNB1/TTBY(N27,NTTBY2)

      CHARACTER * 30 EASS



      IF (IATOTY == -10) THEN
      DO 389 I=1,N27
         DO 389 J=1,NTTBY2
            TTBY(I,J)=0.
  389    CONTINUE
      GO TO 9999
      END IF



!__ COMPUTE INTEGRANDS

!c---+
      tbns=0.
      tbnr=0.
!c----
      DO 10 I=1,IREM

      DRN=FST251(I)
      DRN2=DRN*DRN
      DRN3=DRN2*DRN
      DRN4=DRN3*DRN
      DMN=F251(I)
      RAD=1.-2.*DGRAV*DMN/DRN
      SRAD=SQRT(RAD)
      SRAD3=SRAD**3
      DEN=FENG1(I)
      DPN=FOUT(I)

      AXBN1(I)=DRN2 * PDENEX(I) / SRAD
!c---+  compute and store log A_s(r)
      if(i < irem) then
      tbns = tbns + axbn1(i) * (fst251(i+1)-fst251(i))
      RLASPHrad=3.*DLOG10(AR0) + ALOG10(VPI*tbns)
      Asph_i(i)=RLASPHrad
      r_meter=fst251(i)*(R0R0/EE18)*1.E03
      Rsph_i(i)=r_meter
!      print *, ' i=',i,' r=',r_meter,' log_10 A_s=',RLASPHrad
      end if
!c----
      B1=(DEN+DPN) * P0SD(I) * (  DEDP(I)*(1./SRAD-1.)
     +                          - DEINDP(I)/SRAD )
      B2=(DEN-EINTHX(I)) * ( RM0D(I)*DGRAV/DRN +
     +                       (RLJ(I)*DRN*OGDR(I))**2/3. )
     +                   / SRAD3
      B3=( (RLJ(I)*DRN2*DOGDR(I))**2/12. -
     +     (DJ2DR(I)*DRN3*OGDR(I)**2)/3. )
     +   / (VPI*A*DGRAV*DRN2)
      AXBN2(I)=(B1 + B2 - B3) * DRN2
!c---+  compute and store log A_rot(r) and M_rot(r) 
      if(i < irem) then
      tbnr = tbnr + axbn2(i) * (fst251(i+1)-fst251(i))
      DGRAV =(GRAV/ufe30)*(RMSUN*ufe30)/(UF1*AR0)
      AGD   =A*DGRAV
      ZJDrad = FST251(i)**4 * DOGDR(i) / (6.*AGD)
      DELEBrad = - (ZJDrad*A)**2 * DGRAV / FST251(i)**3 + VPI*A*tbnr
      RMNTH=AMNEUT/(RMSUN * ufe30)
      DELMHrad=RM0D(i) + ZJDrad**2 * AGD*A / FST251(i)**3
      RLGDArad=ALOG10(DELMHrad + DELEBrad) - DLOG10(RMNTH)
      ASPHQrad=10.**(RLASPHrad/4.)
      ASPHQrad=ASPHQrad**4
      if (RLGDArad /= 0.) then
         DArad  =10.**(RLGDArad/4.)
         DArad  =DArad**4
         else
         DArad=0.
      end if
      ATOTrad=ASPHQrad+DArad 
      RLATOTrad=DLOG10(ATOTrad)
!      r_meter=fst251(i)*(R0R0/EE18)*1.E03
      Arot_i(i)=RLATOTrad
!      print *, ' i=',i,' r=',r_meter,' lg A_s=',RLASPHrad, 'lg A_rot=',
!     +RLATOTrad
      DELMH_i=RM0D(i) + ZJDrad**2 * AGD*A / FST251(i)**3
      RMTOTH_i(i)=F251(i)+DELMH_i
      end if
!c----

   10 CONTINUE


!__ COMPUTE  LOG [A_s]
      EASS='BYNUMB: CALL SELSPL 10        '
!                                                ****
      CALL IGRN(FST251,AXBN1,IREM,X251,Y251,IREM,ASPH,C,
     +          NHU,10,1,AUX51,EASS)
!     RETURN:  ASPH
      IF (ASPH > 0.) THEN
      RLASPH=3.*DLOG10(AR0) + ALOG10(VPI*ASPH)

      ELSE
      RLASPH=0.
      END IF
      TTBY(MLF,2)=RLASPH


!__ COMPUTE CHANGE OF BINDING ENERGY  E(B)
      EASS='BYNUMB: CALL SELSPL 20        '
!                                                *****
      CALL IGRN(FST251,AXBN2,IREM,X251,Y251,IREM,DELEB,C,
     +          NHU,10,1,AUX51,EASS)
!     RETURN:  DELEB
      RDSPH=FST251(IREM)
      DELEB= - (ZJD*A)**2 * DGRAV / RDSPH**3 + VPI*A*DELEB
      TTBY(MLF,1)=DELEB


!__ COMPUTE CHANGE OF BARYON NUMBER (DUE TO ROTATION)
      RMNTH=AMNEUT/(RMSUN * ufe30)
      RLGDA=ALOG10(DELMH + DELEB) - DLOG10(RMNTH)

      TTBY(MLF,3)=RLGDA


!__ TOTAL NUMBER OF BARYONS: (A_TOT=A_s+DELTA_A)
      IF (RLASPH /= 0.) THEN
      ASPHQ=10.**(RLASPH/4.)
      ASPHQ=ASPHQ**4
      ELSE
      ASPHQ=0.
      END IF

      IF (RLGDA /= 0.) THEN
      DA  =10.**(RLGDA/4.)
      DA  =DA**4
      ELSE
      DA=0.
      END IF
      TTBY(MLF,5)=SNGL (DA/ASPHQ)

      ATOT=ASPHQ+DA 
      TTBY(MLF,6)=SNGL (ATOT / 1.e30 / 1.e27)


      IF (ATOT > 0.) THEN
      RLATOT=DLOG10(ATOT)

      ELSE
      RLATOT=0.
      END IF
      TTBY(MLF,4)=RLATOT


      ISTOP=10
      IF (ISTOP < 0) STOP


 9999 RETURN

      END




!IGRN                                *****
      SUBROUTINE IGRN(X,Y,N,XX,YY,NN,VALUE,C,IC,IREF,IOPT,AUX51,EASS)


! ----------------------------------------------------------------------
!
!__ PURPOSE:
!              TO INTEGRATE A FUNCTION WHICH IS SPECIFIED NUMERICALLY AT
!              N => 4 POINTS
!
!              IREF < 0: USE CUBIC SPLINE APPROXIMATION,
!              IREF > 0: USE POLYNOMIAL APPROXIMATION TO COMPUTE THE
!                        INTEGRAND YY AT NN GRID POINTS
!
!              IOPT : SEE  _INTER_
!
!
!   NOTE:   NN-1  ==  INTEGRAL NUMBER!
!
!   INPUT: GRID POINTS  X(N), ORDINATE VALUES  Y(N);
!          XX(NN), YY(NN) ARE AUXILIARY ARRAYS
!
!   RETURN: ON RETURN VALUE CONTAINS THE ESTIMATE OF THE INTEGRAL
!
! ----------------------------------------------------------------------


      REAL X(N),Y(N),XX(NN),YY(NN),RESLT(1),C(4,IC),AUX51(1,N)


      CHARACTER * 30 EASS


      VALUE = 0.
      NNNN = NN


! ***********************************************************************
!__ OPEN TEST SEGMENT
      IF ( (N-4) < 0 ) GO TO 100
      IF ( (IREF < 0) .AND. (IC /= N) ) GO TO 100
      GO TO 76
!__ CLOSE TEST SEGMENT
! ***********************************************************************



!__ HANDLE ERROR DETECTION
  100 WRITE (6,800) N,NN,DEL,IREF,IC
  800 FORMAT (///,2X,'N=',I5,2X,'NN=',I5,2X,'DEL=',E12.6,2X,'IREF=',
     +       I3,2X,'IC=',I4,////)
      DO 101 I=1,N
  101 WRITE (6,102) I,X(I),Y(I)
  102 FORMAT (2X,'I=',I4,2X,'X=',E13.6,2X,'Y=',E13.6)
      WRITE (6,'(//,'' _IGRN_ FATAL ERROR --> ERROR STOP'')')
      STOP 'ERROR STOP'


   76 CONTINUE
!__ COPY  Y VALUES
      DO 222 I=1,N
      AUX51(1,I) = Y(I)
  222 CONTINUE


! NNNN  ==  ODD NUMBER!
      IF ( ((-1)**NNNN) > 0 ) NNNN = NNNN-1

      DEL = (X(N)-X(1)) / FLOAT(NNNN-1)
      IF (DEL <= 0.) GO TO 100

      DO 200 I=1,NNNN
      XX(I) = X(1) + DEL*FLOAT(I-1)
  200 CONTINUE



      IF (IREF < 0) THEN

!__ COMPUTE INTEGRAND CONTAINED IN Y(N) AT GRID POINTS XX(N) BY A CUBIC
!   SPLINE APPROXIMATION
!                               **
      CALL SPLINE(X,Y,N,C,IC,XX,YY,NNNN,10,EASS)

      ELSE IF (IREF > 0) THEN

!__ USE A POLYNOMIAL APPROXIMATION
      DO 300 K=1,NNNN
      P = XX(K)
      N1 = 1
!                               *****
      CALL SELINT(X,AUX51,1,1,N,RESLT,P,IOPT)
      YY(K) = RESLT(N1)
  300 CONTINUE

      ELSE

      GO TO 100

      END IF



!                     X(N)
!__ COMPUTE INTEGRAL   I  D[X] Y(X)
!                     X(1)

      VALUE = SIMP(DEL,NNNN,YY)

      RETURN

      END
!SPLINE                                  ****
      SUBROUTINE SPLINE(X,Y,NX,C,IC,XOUT,YOUT,IOUT,IREF,EASS)


! ----------------------------------------------------------------------
!__ PURPOSE:
!              TO COMPUTE AN INTERPOLATORY APPROXIMATION OF A SET OF POINTS
!              BY A CUBIC SPLINE AND
!
!           A) IF IREF<10 COMPUTE ONLY SPLINE COEFFICIENTS
!           B) IF IREF=10 EVALUATE THE CALCULATED SPLINE AT IOUT POINTS
!           C) IF IREF>10 EVALUATE THE FIRST DERIVATIVE OF THE CUBIC SPLINE
!              AT IOUT POINTS
!
! SPLINE MATRIX:
!     C:  MATRIX C(4,IC), CONTAINS THE SPLINE COEFFICIENTS (OUTPUT)
!     IC: ROW DIMENSION OF MATRIX C (INPUT), IC=NX
!
!
! ERROR INFORMATION:
!         IER=129  ><  IC#NX
!         IER=130  ><  NX<2
!         IER=131  ><  ABSCISSA ARE NOT ORDERED SO THAT X(1)<X(2)< ... <X(NX)
!         IER=33   ><  XOUT(I) < X(1)
!         IER=34   ><  XOUT(I) > X(NX)
! ----------------------------------------------------------------------


      REAL X(NX),Y(NX),C(4,IC),XOUT(IOUT),YOUT(IOUT)

      CHARACTER * 30 EASS

      DATA RLIM/1.E-02/

      NAMELIST/NAME1/IER,NXNL,ICNL,IREFNL
      NXNL   = NX
      ICNL   = IC
      IREFNL = IREF
      IV4    = NX

!//////////////////////// BEGIN OF TEST SEGMENT \\\\\\\\\\\\\\\\\\\\\\
      IER = 0
      IF (IC /= NX) IER = 129
      IF (NX < 2)  IER = 130
      CALL SORTOR(X,NX,-10,IERSO)
      IF (IERSO /= 0) IER = 131
! 08/02/1993: Begin of changes
         IF (IER == 131 .AND. NX == 5) THEN
            IER = 0
            IV4 = NX-1
            DO 498 I=1,NX-1
               IF (X(I) == X(I+1)) THEN
               IF (I == NX-1) GO TO 598
               DO 398 J=I+1,NX
                  X(J) = X(J+1)
  398             Y(J) = Y(J+1)
                     END IF
  498                   CONTINUE
  598                     CALL SORTOR(X,IV4,-10,IERSO)
                             END IF
! 08/02/1993: End of changes 


      IF (XOUT(1) < X(1)) THEN
      CABSL=ABS((XOUT(1)-X(1))/X(1))
         IF (CABSL > RLIM) THEN
         IER = 33
         ELSE
         XOUT(1)=X(1)
            END IF
      ELSE IF (XOUT(IOUT) > X(IV4)) THEN
      CABSU=ABS((XOUT(IOUT)-X(NX))/X(IV4))
         IF (CABSU > RLIM) THEN
         IER = 34
         ELSE
         XOUT(IOUT)=X(IV4)
            END IF
                   END IF
!////////////////////////// END OF TEST SEGMENT \\\\\\\\\\\\\\\\\\\\\\


      IF (IER) 10,10,2000

 2000 WRITE (6,NAME1)
      DO 15 I=1,IV4
   15 PRINT*,' XIN=',X(I),' YIN=',Y(I),'  I=',I
      DO 16 I=1,IOUT
   16 PRINT*,' XOUT=',XOUT(I),' YOUT=',YOUT(I),'  I=',I
      PRINT*,' WARNING:  _SPLINE_ FATAL ERROR DETECTION  '
      WRITE (6,'(/,1X,''LOCATION IN PROGRAM: '',A30)') EASS
      STOP


!__ COMPUTE SPLINE COEFFICIENTS OF CUBIC SPLINE FUNCTIONS
   10 CONTINUE
!                  NX      *
      CALL CSPLINE(IV4,X,Y,C)
      IF (IREF-10) 30,40,50
   30 RETURN


!__ EVALUATE THE CALCULATED SPLINE AT IOUT POINTS (CONTAINED IN XOUT)
   40 CONTINUE
!                   NX            ****
      CALL SEVU(X,Y,IV4,C,IC,XOUT,YOUT,IOUT)
!     RETURN: YOUT_IOUT
      RETURN


!__ EVALUATE THE FIRST DERIVATIVE OF THE CUBIC SPLINE AT IOUT POINTS.
   50 CONTINUE
!                    NX             ****
      CALL DSEVU(X,Y,IV4,C,IC,XOUT,YOUT,IOUT)
!     RETURN: YOUT_IOUT = D Y_IOUT/D X


      RETURN

      END


!SPLINE8                                  ****
      SUBROUTINE SPLINE8(X,Y,NX,C,IC,XOUT,YOUT,IOUT,IREF,EASS)


! ----------------------------------------------------------------------
!__ PURPOSE: DOUBLE PRECISION VERSION OF  SPLINE
!
!              TO COMPUTE AN INTERPOLATORY APPROXIMATION OF A SET OF POINTS
!              BY A CUBIC SPLINE AND
!
!           A) IF IREF<10 COMPUTE ONLY SPLINE COEFFICIENTS
!           B) IF IREF=10 EVALUATE THE CALCULATED SPLINE AT IOUT POINTS
!           C) IF IREF>10 EVALUATE THE FIRST DERIVATIVE OF THE CUBIC SPLINE
!              AT IOUT POINTS
!
! SPLINE MATRIX:
!     C:  MATRIX C(4,IC), CONTAINS THE SPLINE COEFFICIENTS (OUTPUT)
!     IC: ROW DIMENSION OF MATRIX C (INPUT), IC=NX
!
!
! ERROR INFORMATION:
!         IER=129  ><  IC#NX
!         IER=130  ><  NX<2
!         IER=131  ><  ABSCISSA ARE NOT ORDERED SO THAT X(1)<X(2)< ... <X(NX)
!         IER=33   ><  XOUT(I) < X(1)
!         IER=34   ><  XOUT(I) > X(NX)
! ----------------------------------------------------------------------


      IMPLICIT REAL * 8 (A-H,O-Z)

      REAL * 8 C(4,IC)
      REAL * 8 X(NX),Y(NX),XOUT(IOUT),YOUT(IOUT)

      CHARACTER * 30 EASS

      NAMELIST/NAME1/IER,NXNL,ICNL,IREFNL
      NXNL   = NX
      ICNL   = IC
      IREFNL = IREF

!__ OPEN TEST SEGMENT
      IER = 0
      IF (IC /= NX) IER = 129
      IF (NX < 2)  IER =130
! no access  CALL SORTOR(X,NX,-10,IERSO)
! no access  IF (IERSO /= 0) IER = 131
      DO 555 I=1,IOUT
      IF (XOUT(I) < X(1))  IER = 33
  555 IF (XOUT(I) > X(NX)) IER = 34
!__ CLOSE TEST SEGMENT

      IF (IER) 10,10,2000

 2000 WRITE (6,NAME1)
      DO 15 I=1,NX
   15 PRINT*,' XIN=',X(I),' YIN=',Y(I),'  I=',I
      DO 16 I=1,IOUT
   16 PRINT*,' XOUT=',XOUT(I),' YOUT=',YOUT(I),'  I=',I
      PRINT*,'  _SPLINE8_  -->  FATAL ERROR DETECTION   '
      WRITE (6,'(/,1X,''LOCATION IN PROGRAM: '',A30)') EASS
      STOP


!__ COMPUTE SPLINE COEFFICIENTS OF CUBIC SPLINE FUNCTIONS
   10 CONTINUE
!                       *
      CALL CSPL8(NX,X,Y,C)
      IF (IREF-10) 30,40,50
   30 RETURN


!__ EVALUATE THE CALCULATED SPLINE AT IOUT POINTS (CONTAINED IN XOUT)
   40 CONTINUE
!                                 ****
      CALL SEVU8(X,Y,NX,C,IC,XOUT,YOUT,IOUT)
!     RETURN: YOUT_IOUT
      RETURN


!__ EVALUATE THE FIRST DERIVATIVE OF THE CUBIC SPLINE AT IOUT POINTS.
   50 CONTINUE
!                                  ****
      CALL DSEVU8(X,Y,NX,C,IC,XOUT,YOUT,IOUT)
!     RETURN: YOUT_IOUT = D Y_IOUT/D X


      RETURN

      END

!INTER                                   ***
      SUBROUTINE INTER(IOPT,X,Y,IY,N1,N2,ERG,P)

! --------------------------------------------------------------------
!
! PURPOSE:
!          INTER IS AN INTERPOLATION ROUTINE WHICH ALLOWS FOR
!          DIFFERENT CHOICES OF INTERPOLATING OR EXTRAPOLATING A
!          FUNCTION.
!          IOPT DETERMINES THE METHOD OF INTERPOLATION AS FOLLOWS:
!
!  IOPT:(<0)     LOGARITHMIC INTERPOLATION
!       (=1)     LINEAR INTERPOLATION
!       (=2)     POLYNOMIAL INTERPOLATION
!       (=3)     SPLINE INTERPOLATION
!       (=10+NZ) POLYNOMIAL INTERPOLATION AT THE NEXT (TO NODE P) NZ
!                GRID POINTS
!
!  X(N2):     GRID POINTS XI; I=1,N2
!  Y(IY,N2):  Y VALUES AT THE GRID POINTS
!  IY:        FIRST DIMENSION OF Y IN THE CALLING PROGRAM
!  N1:        NUMBER OF FUNCTIONS (CURVES);  N1 <= IY
!  N2:        NUMBER OF GRID POINTS
!  ERG(N1):   CONTAINS THE RESULT
!  P:         X VALUE AT WHICH THE INTERPOLATION IS TO BE PERFORMED
! --------------------------------------------------------------------



      PARAMETER (NGES=112)
      PARAMETER (NMMM=2*NGES)
      PARAMETER (LCK=NMMM+4,LWRK=6*NMMM+16)

      DIMENSION X(N2),Y(IY,N2),ERG(IY),W1(NMMM),Y1(NMMM)

      COMMON/MISCPARAMETERS/LIM007,ICRO22,IINE9,DZER9L,DZER9P

      

! ***********************************************************************
!__ OPEN TEST SEGMENT
      IF (P < X(1) .OR. P > X(N2)) THEN
         IF ( IINE9 <= LIM007 ) THEN
         WRITE (6,'(/,''  _INTER_  -->  P OUT OF RANGE:'')')
         WRITE (6,499) P,X(1),X(N2)
  499    FORMAT (' P=',E13.6,' X(1)=',E13.6,' X(N2)=',E13.6)
         END IF
         IINE9 = IINE9+1
      END IF
!__ CLOSE TEST SEGMENT
! ***********************************************************************

      NERR=999
      S=0.

      NTEST = 10

      IF (N1 > IY) THEN
         NERR = 5
         GO TO 140

      END IF

      IF (IOPT < 0) THEN

!  LOGARITHMIC INTERPOLATION

         IF (IY /= 1) THEN
            WRITE (6,'(///,'' IY /= 1 IS NOT ALLOWED: IY='',I4,//)') IY
            GO TO 140
         END IF

!  SEARCH FOR THE GRID POINTS  X(J)  AND  X(J+1), SUCH THAT
!  X(J) < P < X(J+1)
         N2M1=N2-1

         DO 370 J=1,N2M1
         JJ=J
         IF ( (P >= X(J)) .AND. (P <= X(J+1)) ) GO TO 371
         IF (J == N2M1) THEN
         WRITE (6,'(///,'' ERROR STOP IN SUBROUTINE INTER: P OUT OF'',
     +             '' RANGE'',/,'' LOGARITHMIC INTERPOLATION FAILED'',
     +             /,'' X(1)='',E12.6,'' X(N2)='',E12.6,'' P='',E12.6
     +             ,'' N2='',I4)') X(1),X(N2),P,N2
         GO TO 140
         END IF
  370    CONTINUE

  371    JJP1=JJ+1

!  INTERPOLATE
         TPZ=ALOG(Y(IY,JJP1)) - ALOG(Y(IY,JJ))
         TPN=ALOG(X(JJP1)) - ALOG(X(JJ))
         TPF=ALOG(P) - ALOG(X(JJ))

         RLGP=ALOG(Y(IY,JJ)) + TPZ * TPF / TPN
         ERG(IY)=EXP(RLGP)


      ELSE IF (IOPT == 1) THEN


!  LINEAR INTERPOLATION

         DO 30 N = 1,N1
!  SEARCH FOR THAT GRID POINT WHICH IS CLOSEST TO  P
            NMIN1 = 1
            ABW = ABS(X(1)-P)
            DO 10 K = 1,N2
               DIFF = ABS(X(K)-P)
               IF (DIFF < ABW) THEN
                  NMIN1 = K
                  ABW = DIFF
               END IF

   10       CONTINUE

!  SEARCH FOR THAT GRID POINT WHICH IS THE SECOND CLOSEST POINT TO  P

            IF (NMIN1 /= 1) THEN
               NMIN2 = 1
               ABW = ABS(X(1)-P)

            ELSE
               NMIN2 = 2
               ABW = ABS(X(2)-P)
            END IF

            DO 20 K = 1,N2
               IF (K == NMIN1) GO TO 20
               DIFF = ABS(X(K)-P)
               IF (DIFF < ABW) THEN
                  NMIN2 = K
                  ABW = DIFF
               END IF

   20       CONTINUE

            DXXXX = X(NMIN1)-X(NMIN2)
            IF (ABS(DXXXX) < DZER9L) DZER9L = ABS(DXXXX)
            ERG(N) = (Y(N,NMIN1) * (P-X(NMIN2)) -
     +               Y(N,NMIN2) * (P-X(NMIN1))) / DXXXX
   30    CONTINUE
         IF (NTEST == 0) RETURN
         IF (NMIN1 == NMIN2) THEN
            NERR = 1
            GO TO 140

         ELSE IF (ABS(X(NMIN1)-X(NMIN2))  <=  1.0E-26) THEN
            NERR = 2
            GO TO 140

         ELSE IF (N1 < 1) THEN
            NERR = 3
            GO TO 140

         ELSE IF (N2 < 2) THEN
            NERR = 4
            GO TO 140

         END IF

      ELSE IF (IOPT == 2) THEN

!  POLYNOMIAL INTERPOLATION ACCORDING TO AITKEN NEVILLE

         DO 60 N = 1,N1
            IF (N2 > NMMM) THEN
               NERR = 5
               GO TO 140

            END IF
!  EXTRAPOLATION TABLE

            DO 50 I = 1,N2
               W1(I) = Y(N,I)
               DO 40 K = I - 1,1,-1
                  DXXXX = X(I)-X(K)
                  IF (ABS(DXXXX) < DZER9P) DZER9P = ABS(DXXXX)
                  W1(K) = W1(K+1) + (W1(K+1)-W1(K)) * (P-X(I))/
     +                    DXXXX
   40          CONTINUE
   50       CONTINUE
            ERG(N) = W1(1)
   60    CONTINUE

!  TESTS

         IF (NTEST == 0) RETURN
         IF (N1 < 1) THEN
            NERR = 3
            GO TO 140

         ELSE IF (N2 < 2) THEN
            NERR = 4
            GO TO 140

         END IF

      ELSE IF (IOPT == 3) THEN
         IF (N2 > NMMM) THEN
            NERR = 5
            GO TO 140

         END IF


!  SPLINE INTERPOLATION USING NAG-ROUTINE E01BAF

         DO 80 N = 1,N1
!  LOADING OF Y1
            DO 70 M = 1,N2
               Y1(M) = Y(N,M)
   70       CONTINUE
            IFAIL = 0
! NO ACCESS CALL E01BAF(N2,X,Y1,RK,C,LCK,WRK,LWRK,IFAIL)
            IF (IFAIL /= 0) THEN
               NERR = 6
               GO TO 140

            END IF

! NO ACCESS CALL E02BBF(N2+4,RK,C,P,S,IFAIL)
            IF (IFAIL /= 0) THEN
               NERR = 7
               GO TO 140

            END IF


            ERG(N) = S
   80    CONTINUE

      ELSE IF (IOPT > 10) THEN

!  POLYNOMIAL INTERPOLATION AT THE NZ GRID POINTS WHICH ARE NEXT TO P

         NZ = IOPT - 10
         IF ((NZ < 2) .OR. (NZ > N2) .OR. (NZ > NMMM)) THEN
            NERR = 8
            GO TO 140

         END IF

         DO 130 N = 1,N1

!  SEARCH FOR THE NEXT GRID POINT, AND TEST IF GRID POINTS ARE ORDERE

            NMIN1 = 1
            ABW = ABS(X(1)-P)
            DO 90 K = 2,N2
               DIFF = ABS(X(K)-P)
               IF (DIFF < ABW) THEN
                  ABW = DIFF
                  NMIN1 = K
               END IF

               IF (X(K) <= X(K-1)) THEN
                  NERR = 9
                  GO TO 140

               END IF

   90       CONTINUE
            NMIN2 = NMIN1

!  SEARCH FOR THE NEXT NZ-1 GRID POINTS

            DO 100 K = 1,NZ - 1
               IF (NMIN1 == 1) THEN
                  NMIN2 = NMIN2 + 1

               ELSE IF (NMIN2 == N2) THEN
                  NMIN1 = NMIN1 - 1

               ELSE
                  ABW1 = ABS(X(NMIN1-1)-P)
                  ABW2 = ABS(X(NMIN2+1)-P)
                  IF (ABW1 < ABW2) THEN
                     NMIN1 = NMIN1 - 1

                  ELSE
                     NMIN2 = NMIN2 + 1
                  END IF

               END IF

  100       CONTINUE

!  AITKEN-NEVILLE

            DO 120 I = 1,NZ
               W1(I) = Y(N,NMIN1-1+I)
               DO 110 K = I - 1,1,-1
                  W1(K) = W1(K+1) + (W1(K+1)-W1(K))* (P-X(NMIN1-1+I))/
     +                    (X(NMIN1-1+I)-X(NMIN1-1+K))
  110          CONTINUE
  120       CONTINUE
            ERG(N) = W1(1)
  130    CONTINUE

      ELSE
         WRITE (6,'(''  _INTER_   --> FATAL ERROR '')')
         STOP
      END IF


      RETURN

!  05/23/2015: added 4 lines additional after 140 continue      
 140  continue
      if (x(nmin1) == x(nmin2)) then
         p = y(iy,nmin1)
         return
      end if

      WRITE (6,'(''  _INTER_   --> FATAL ERROR  '')')
      WRITE (6,9000) NERR, IOPT, X(NMIN1), X(NMIN2), Y(IY,NMIN1),
     +               Y(IY,NMIN2), NMIN1, NMIN2, N1, N2, IY, P
      STOP

 9000 FORMAT ('  NERR= ',I3,' IOPT= ',I3, /,'  X(NMIN1),  X(NMIN2),  Y('
     +        'IY,NMIN1),  Y(IY,NMIN2):', /, 4 (1PE15.6,2X),  /, '  NMI'
     +        'N1,  NMIN2:', 2 (I3,2X), / ,'  N1,  N2,  IY: ', 3 (I3,2X)
     +        , '  P: ', 1PE15.6) 
      END

!INTV                           ***
      SUBROUTINE INTV(XA,XB,NAB,AXK,IREF,ARGINT)


! ----------------------------------------------------------------------
!
! PURPOSE:
!          CALCULATION OF THE GRID POINTS OF INTERVAL XB-XA, (XB>XA)
!
!          IREF=0: AEQUIDISTANT GRID POINTS
!          IREF=1: WITH SIN(X) WEIGHTED GRID POINTS
!
! RETURN:
!         ARRAY AXK CONTAINS THE GRID POINTS
!         NOTE: AXK(1)=XA, AXK(NAB)=XB
!
! ----------------------------------------------------------------------


      DIMENSION AXK(NAB)

      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2

      CHARACTER * 30 ARGINT


!__ TEST BOUNDARIES
      IF (XA >= XB) GO TO 17


      IF (IREF == 0) THEN

!__ AEQUIDISTANT GRID POINTS
      DEL=(XB-XA)/FLOAT(NAB-1)
      DO 1000 I=1,NAB
      AXK(I)=XA+DEL*FLOAT(I-1)
 1000 CONTINUE
      IF ( ABS(AXK(NAB)-XB)/XB  >=  1.E-03 ) GO TO 17
      AXK(NAB)=XB
      RETURN


      ELSE IF (IREF == 1) THEN

!__ COMPUTE WEIGHTED GRID POINTS
      DO 700 I=1,NAB
      ZZ=(FLOAT(I-1)/FLOAT(NAB-1)) * PID2
  700 AXK(I)=XA+(XB-XA)*SIN(ZZ)
      IF ( ABS(AXK(NAB)-XB)  >=  1.E-03 ) GO TO 17
      RETURN


      ELSE

!__ ERROR RETURN
      PRINT*,'  _INTV_  -->  FATAL ERROR DETECTED;  IREF=', IREF
      GO TO 17
      END IF



   17 CONTINUE
      DO 2000 I=1,NAB
 2000 WRITE (6,20) I,NAB,XA,XB,DEL,AXK(I)
      WRITE (6,10) ARGINT
   10 FORMAT (//, 1X,'  _INTV_  --> ERROR STOP **',/,'  CALLING ROUT',
     +           'INE:', /, 18X, A30, ///)
   20 FORMAT (3X,I4,2X,I4,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5)
      STOP

      END


!SIMP
      FUNCTION SIMP(DEL,N,A)


! ----------------------------------------------------------------------
!__ PURPOSE:
!              APPLY SIMPSON'S RULE OF INTEGRATION.
! NOTE: ON ENTRY, A(N) MUST CONTAIN THE VALUES OF THE INTEGRAND;
!       UNCHANGED ON EXIT.
!
!     DEL = (B - A) / (N - 1);  N-1  ==  INTEGRAL NUMBER
! ----------------------------------------------------------------------


      REAL A(N)


      L = N
      SUM = A(1) - A(L)

!__ INTEGRATION LOOP
      DO 1 I=2,L,2
    1 SUM = SUM + 4.0*A(I) + 2.0*A(I+1)


      SIMP =  (DEL * SUM) /3.0


      RETURN

      END

!RNSNEW                                                       **** **
      SUBROUTINE RNSNEW(EDDY,PRESS,IEP,EC,AUX3,AUX4,NHU,IFAIL,RQOZ,WI,C,
!                                             ***                     **
     +                  AUX31,AUX41,OMEGA,P,R,RTH,T,E,NI,NIP1,NJ,NJP1,IO
!                       **
     +                  V1,MMM,IREM,BORD0,BORD1)


!-----------------------------------------------------------------------
!
! PURPOSE:
!          ROTATING QUARK STAR IN THE NEWTONIAN LIMIT:
!          SOLVE THE NON-RELATIVISTIC REDUCTION OF THE OPPENHEIMER-
!          VOLKOFF EQUATIONS OF STELLAR STRUCTURE FOR A  R O T A T I N G
!          STAR (i.e., AS A NEW FEATURE ONE ENCOUNTERS A NON-HOMOGENEOUS
!          MASS DISTRIBUTION); ROTATIONAL SYMMETRY OF THE STAR AROUND
!          THE ROTATIONAL Z-AXIS REMAINS VALID (PHI-DIRECTION).
!
!
! INPUT:
!        - EOS P(EPSILON)  (IN MeV/fm**3): EDDY,PRESS
!        - CENTRAL ENERGY DENSITY  EC=EPSILON_C  (IN MeV/fm**3)
!          (150 MeV/fm**3 < E_C < 5000 MeV/fm**3)
!        - ROTATIONAL FREQUENCY  OMEGA  OF THE STAR (IN  1/fm)
!
!
! CONSTANTS:
!        GRAVITATIONAL CONSTANT: GRAV (IN fm**2)
!        SCALING FACTOR: R0R0 (IN fm)
!        MASS OF THE SUN:  RMSUN (IN MeV)
!
!
! RETURN:
!         RTH=R(NS;SURFACE;0<=THETA<=#) (IN METER)
!         RQOZ=M_NS/M_SUN
!         WI= LOG (MOMENT OF INERTIA OF THE ROTAING STAR / (g cm**2) )
!
!         IFAIL<0  <=>  ERROR RETURN
!-----------------------------------------------------------------------



      REAL * 8 DGRAV,A,AR0,EE18,GCM2KM,EEM18,EE14,EE34,EE44
      REAL * 8 ee03,ee13,ee15      
      REAL * 8 RMSUN,RMNS,GRAV,GRAVK,GRAVKM

      DIMENSION EDDY(NHU),PRESS(NHU),AUX3(NHU),AUX4(NHU),C(4,NHU)
      DIMENSION AUX31(1,NHU),AUX41(1,NHU),R(0:NJ,0:NI),P(0:NJ,0:NI)
      DIMENSION T(0:NI),RTH(0:NI,4),E(0:NJ,0:NI)

      PARAMETER (KP1=1)
      DIMENSION EIP(KP1),PIP(KP1)

      COMMON/ExpoMaxMin/EXMAX,EXMIN
      COMMON/POWERSTEN/EE03,EE15,EE18,EEM18,EE13,EE14,EE34,EE44
      common/ConversionFactor30/ufe30
      COMMON/MULTIPI/PI,PIQ,ZPI,VPI,ZPI3,PID2
      COMMON/CONVERSIONFACTORS/UF1,UF12,UF13,UF4,UF6,UF10
      COMMON/NUCLEONMASSES/AMNUC,AMNEUT
      COMMON/GRAVITATIONAL/GRAVK,GRAV,GRAVKM,GCM2KM
      COMMON/SOLARPROPERTIES/RMSUN,RSUN,RSUNKM
      COMMON/MiscInt1/I100,I7,I27,IVAR,INPOG1,IKP1,ICRUST
      COMMON/MISCINT2/J651,J650,IGIVMX,IVONR,KC20
      COMMON/RADIALSCALINGFACTOR/R0R0

      PARAMETER (N27=50,N10=15)
      COMMON/TABN1/TN1(N27,N10)

      PARAMETER (m651=160501)
      COMMON/CNSTA1/F251(M651),FOUT(M651),FST(M651),FST251(M651),
     +              FENG1(M651),X251(M651),Y251(M651),F251PM(M651),
     +              FDRIP(M651)
      COMMON/CNSTA2/AUX51(1,M651),AUX5(M651)
      COMMON/CNSTA3/AUX61(1,M651),AUX6(M651),DPDRAX(M651)
      COMMON/CNSTA4/AUX71(1,M651),AUX7(M651),AUX8(M651),AUX9(M651)
      COMMON/SPL1/ST(M651),YOUT(M651),XSTROB(M651),YSTROB(M651)

      CHARACTER * 30 ARGINT
      CHARACTER * 30 EASS


      DATA RMC/0./
      DATA FMKM/1.E18/,FMEVFM/5.61012/

      NAMELIST/NAME1/RMNS,RMSUN,FOUT,FST251 



!.................... AUXILIARY FUNCTIONS ........................

!            ROTATING (PSEUDO) OPPENHEIMER-VOLKOFF STAR

!                           D [P] / D [R]

!     FOVRAD(DEN,DPN)=-DGRAV*(DEN+DPN)*(DMN+VPI*A*DRN3*DPN)
!    +                /(DRN2*(1.-F2EFF*DGRAV*DMN/DRN))
!    +                +DRN*(DEN+0. )*DOM2*SINT*SINT

      FOVRAD(DEN,DPN)=(-DGRAV*DEN*DMN/DRN2 + DRN*DEN*DOM2*SINT*SINT)
     +                * (1.+DPN/DEN) * (1.+VPI*A*DRN3*DPN/DMN)
     +                / (1.-2.*DGRAV*DMN/DRN)

!                           ROTATING NEWTONIAN STAR

!                             D [P] / D [R]

      FH(DEN)=-DGRAV*DEN*DMN/DRN2 + DRN*DEN*DOM2*SINT*SINT

!                            D [P] / D [THETA]

      FT(DEN)= DOM2*DRN2*SIN2T*DEN/2.

!                     2
!                    D  [M] / D [R] D [THETA]

      FM(DEN)=ZPI*A*DRN2*SINT*DEN

!.................. END OF AUXILIARY FUNCTIONS ...................



      WRITE (6,'('' _RNSNEW_ -->  ROTATING STAR CALCULATION'',
     +          '';  MMM='',I3)') MMM

      IOV1 =-10
      ISEC1=0
      ISEC2=0
      KT=1
      ISPH=10
      IDIFF1=0
      IMESS1=0
      IFAIL=10
      IF (IEP > NHU) THEN
      IFAIL=-10
      WRITE (6,'('' ** _RNSNEW_ FIELD LENGTH OUT OF RANGE **'')')
      WRITE (6,'(5X,'' IEP='',I4,5X,''NHU='',I4,///)') IEP,NHU
      RETURN
      END IF

!__ INITIALIZE
      DO 1 I=0,NI
      DO 2 J=0,NJ
      R(J,I)=0.
      P(J,I)=0.
    2 E(J,I)=0.
      T(I)=0.
    1 RTH(I,2)=0.


!__ CONSTANTS
      IF (IKP1 /= KP1) THEN
      WRITE (6,'(//,'' ****  _RNSNEW_  IKP1 # KP1  ****'',//,
     +              '' IKP1='',I4,3X,'' KP1='',I4)') IKP1,KP1
      IFAIL=10
      RETURN
      END IF
      AR0=R0R0
      A=EC*AR0**3/(RMSUN*ufe30)
      DGRAV =(GRAV/ufe30)*(RMSUN*ufe30)/(UF1*AR0)
      DRMC=RMC/(RMSUN*ufe30)
      RMETER=R0R0/EE15
      DOM=OMEGA*R0R0
      DOM2=DOM*DOM


!__ ENERGY DENSITY AND PRESSURE (DIMENSIONLESS)
      DO 10 I=1,IEP
      AUX3(I)=EDDY(I)/EC
      AUX4(I)=PRESS(I)/EC
      AUX31(1,I)=AUX3(I)
      AUX41(1,I)=AUX4(I)
   10 CONTINUE
      DEC=EC/EC


!__ COMPUTE GRID FOR RADIAL INTEGRATION,  R=[0,1]
      BORD00=BORD0
      BORD11=BORD1
      ARGINT='RNSNEW: RADIAL GRID - FST251  '
!                                  ******
      CALL INTV(BORD00,BORD11,J651,FST251,1,ARGINT)
!     RETURN: FST251
      FST251(1)=FST251(2)/2.


!__ COMPUTE GRID POINTS FOR THETA INTEGRATION
      ABORD=0.
      BBORD=PID2
      ARGINT='RNSNEW: THETA GRID  [0,#/2]   '
      CALL INTV(ABORD,BBORD,NIP1,X251,0,ARGINT)
      DO 530 I=0,NI
  530 T(I)=X251(I+1)




! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   SIMULTANEOUS INTEGRATION OF THE SET OF PARTIAL DIFFERENTIAL EQUATIONS
!   OF A ROTATING STAR (NEWTONIAN DESCRIPTION)
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


!__ STARTING VALUES
      F251(1)=DRMC
      EIP(1)=DEC
      FST(1)=0.
! CENTRAL PRESSURE  PIP=P(C)         ***
      CALL SELINT(AUX3,AUX41,1,1,IEP,PIP,DEC,1)
      DO 723 I=0,NI
      P(1,I)=PIP(IKP1)
  723 E(1,I)=DEC
      DPDRS=-DGRAV*VPI*A*DEC**2*FST251(1)


      LF=2
      YMS=0.

 1010 CONTINUE
      NN2=LF



! ******************************************************************
!                         START INTEGRATING
! ******************************************************************
      DO 1000 N=NN2,J650
! RADIAL LOOP


! ...................... SECURITY COUNTER CHECK ....................
      ISEC1=ISEC1+1
      IF (ISEC1 > J651) THEN
      WRITE (6,'(//'' _RNSNEW_ SECURITY COUNTER LIMIT REACHED!!'',
     +      /,'' ERROR RETURN TO CALLING SUBROUTINE'')')
      IFAIL=-10
      RETURN
      END IF
! ........................... END OF SCC ...........................


      IREM=N-1
      IF (ISPH > 0) THEN
      DMN=F251(N-1)
      ELSE IF (ISPH < 0) THEN
      DMN=DMN+YMS
      END IF

      DRN=FST251(N-1)
      DRN2=DRN*DRN
      DRN3=DRN2*DRN
      HRK=FST251(N)-FST251(N-1)


      KTLL=KT
      DO 1500 I=KTLL,NIP1
! THETA LOOP (RADIAL DISTANCE IS KEPT FIXED)

      IF (I /= NIP1) THEN
      DELTH=T(I)-T(I-1)
      ELSE
      DELTH=T(I-1)-T(I-2)
      END IF

      SINT=SIN(T(I-1))
      SIN2T=SIN(2.*T(I-1))



!__ P(N,I)
      IF (N == 2) THEN
         DPDR=DPDRS
         ELSE
         IF (IOV1 < 0) DPDR=FOVRAD(E(N-1,I-1),P(N-1,I-1))
         IF (IOV1 > 0) DPDR=FH(E(N-1,I-1))
      END IF
      DPDT=FT(E(N-1,I-1))

! ..................................................
      IF (I == KTLL) THEN
         P(N,I-1)=P(N-1,I-1) + HRK*DPDR
         ELSE
         P(N,I-1)=P(N  ,I-2) + DELTH*DPDT
      END IF
! ..................................................


!     COMPUTE ENERGY DENSITY  EPS(.,THETA)
      PIP(IKP1)=P(N,I-1)
      DPN      =PIP(IKP1)

      IF (PIP(IKP1) <= AUX4(1)) THEN
      EIP(IKP1)=AUX3(1)
      ELSE
!                                    ***
      CALL SELINT(AUX4,AUX31,1,1,IEP,EIP,DPN,1)
      END IF

! ..................................................
      E(N,I-1)=EIP(IKP1)
! ..................................................


 1500 CONTINUE


 1020 KTLL=KT

      PNP1=P(N,KTLL-1)
      PNP2=P(N-1,KTLL-1)


! ...................... CHECK  P(N) < P(N-1) ........................
      IF (PNP1 > PNP2) THEN

      IF (ISEC2 > 4) THEN
      WRITE (6,'(///,'' _RNSNEW_ P(N) > P(N-1)!!  ABORT CALCU'',
     +      ''LATION'',/,'' CONTINUE WITH NEXT VALUE FOR THE CE'',
     +      ''NTRAL ENERGY DENSITY'',///)')
      IFAIL=-10
      RETURN
      ELSE
      ISEC2=ISEC2+1
      END IF

      END IF
! ............................. E N D ................................


!>>>>>>
       IPT33=0
       IF (IPT33 /= 0) THEN
       PRINT*,' '
       PRINT*,' N=',N,' KT=',KT,' KTLL=',KTLL
       PRINT*,' PNP1=',PNP1,' PNP2=',PNP2
       END IF

      IF (PNP1 <= 0. .AND. PNP2 > 0.) THEN
!     - P(N-1)>0 .AND. P(N)<0  =>  R(THETA)

         IF (KTLL == 1) THEN
!        - END OF SPHERICAL CALCULATION
         ISPH=-10
         RMSPH=F251(IREM)
         RRSPH=FST251(IREM)*RMETER
         END IF

!     COMPUTE THE STAR'S RADIUS  R(THETA)  FOR A GIVEN VALUE OF THETA
      K44=5
      KHIGH=IREM+1
      KLOW=KHIGH-K44

      DO 5811 II=KLOW,KHIGH
      IRUN=II-KLOW+1
      AUX3(IRUN)=FST251(II)
      AUX41(1,IRUN)=P(II,KTLL-1)
 5811 CONTINUE

      IDIM=KHIGH-KLOW+1
      IOUT=20
      IYY=1
! COMPUTE ZEROE POINT ON FUNCTION AUX31 VS. AUX3         ****
      CALL ZERO22(AUX3,AUX41,IYY,IDIM,X251,AUX51,ST,IOUT,YOUT)
      RTH(KTLL-1,2)=YOUT(IYY)

!>>>>>
      IPT44=0
      IF (IPT44 /= 0) THEN
      PRINT*,' R(OLD)=',FST251(IREM)*RMETER
      AUSD=RTH(KTLL-1,2)*RMETER
      PRINT*,' R=',AUSD,' KTLL-1=',KTLL-1,' N=',N
      END IF

!     - CONTINUE CALCULATION WITH  THETA(I)->THETA(I+1)
      KT=KT+1
      LF=IREM
!>>>>>
      IPT55=0
      IF (IPT55 /= 0)
     +PRINT*,' P(I)>0+P(I+1)<0 R=',RTH(KTLL-1,2),' N=',N,' KT=',KT

      IF (KTLL-NIP1) 1020,2000,2000

      ELSE
!     - SIGN OF P(I) AND P(I+1) OTHERWISE!

      IF (ISPH > 0) GO TO 1515
!     - DEAD CODE AS LONG AS THE SPHERICAL PART OF THE STAR IS TO BE
!       CALCULATED

         IF (PNP1 < 0. .AND. PNP2 < 0.) THEN
!        - P(N)<0 .AND. P(N-1)<0!
         WRITE (6,'(//,'' _RNSNEW_ -NO ACCESS TO THIS NESTED  IF-'',
     +         //)')
         LPR1=-10

         IF (LPR1 > 0) THEN
         WRITE (6,'('' PNP1='',E14.6,'' PNP2='',E14.6,///)') PNP1,PNP2
         DO 7003 KF=0,NI
 7003    WRITE (6,7004) KF,RTH(KF,2)*RMETER,E(N,KF),E(N-1,KF)
 7004    FORMAT (' I=',I4,2X,'R(I)=',E13.6,2X,'E(N-1,I)=',E13.6,2X,
     +   'E(N,I)=',E13.6)
         STOP
         ELSE
         RTH(KTLL-1,2)=RTH(KTLL-2,2)
         KT=KT+1
         IF (KT-(NIP1+1)) 1020,2000,2000
         END IF

         ELSE IF (PNP1 > 0. .AND. PNP2 > 0.) THEN
!        - ADD MASS OF A SPHERICAL SHELL
!        - CONTINUE CALCULAITON WITH  R(N)->R(N+1)
         IDIFF1=1

!>>>>>>>>
         IPT66=0
         IF (IPT66 /= 0) PRINT*,' ADD SPHERICAL MASS SHELL'

         END IF

      END IF



 1515 IF (IDIFF1 < 0) THEN

      FST(N)=EIP(IKP1)*FST251(N)**2
      IF (IMESS1 == 0) THEN
      WRITE (6,'(///'' _RNSNEW_  ***  CALCULATION OF THE STARS'',
     +       '' MASS BY MEANS OF A DIRECT INTEGRATION'',/,'' OVER'',
     +       '' THE ENERGY DENSITY FUNCTION  ***'',///)')
      IMESS1=10
      END IF
      IF (IREM < 4) THEN
      WT=DEC*FST251(N)**3/3.
      ELSE IF (IREM >= 4) THEN

!..INTEGRATE ENERGY DENSITY TO GET STAR'S MASS
!
!          #                   R (THETA)        2
!     WT = I [D THETA] SIN(THETA) I [D R']  (R')   EPSILON(R'*R ;THETA)
!          0                      0                            0

      EASS='RNSNEW: CALL IGRN  10         '
!                                              **
      CALL IGRN(FST251,FST,IREM,X251,Y251,IREM,WT,C,NHU,10,1,AUX51,EASS)
!     RETURN:  WT
      END IF

!__  M(N)
      F251(N)=VPI*A*WT
      GO TO 1000

      ELSE IF (IDIFF1 == 0) THEN
!     - SPHERICAL PART OF THE DEFORMED STAR
!     - M(N) (VIA SOLVING THE DIFFERENTIAL EQUATION FOR THE STAR'S MASS)
      YMS=0.
      DO 1600 I=KTLL,NIP1
      IF (I /= NIP1) THEN
      DELTH=T(I)-T(I-1)
      ELSE
      DELTH=T(I-1)-T(I-2)
      END IF
      SINT=SIN(T(I-1))
      SIN2T=SIN(2.*T(I-1))
      D2DRDT=FM(E(N-1,I-1))
      YMS=YMS+2.*D2DRDT*HRK*DELTH
 1600 CONTINUE
      F251(N)=F251(N-1)+YMS
!>>>>>>>
!>>>>>>>  F251(N)=F251(N-1)+HRK*VPI*A*DRN2*E(N-1,0)
      GO TO 1000

      ELSE IF (IDIFF1 > 0) THEN
!     - SPHERICAL MASS SHELL
      YMS=0.
      DO 1700 I=KTLL,NIP1
      IF (I /= NIP1) THEN
      DELTH=T(I)-T(I-1)
      ELSE
      DELTH=T(I-1)-T(I-2)
      END IF
      SINT=SIN(T(I-1))
      YMS=YMS+VPI*A*DRN2*SINT*E(N-1,I-1)*HRK*DELTH
 1700 CONTINUE
!>>>>>
      IPT77=0
      IF (IPT77 /= 0) PRINT*,' YMS=',YMS,' N=',N
      END IF


      IF (PNP1 > 0. .AND. PNP2 > 0.) THEN
      LF=LF+1
      GO TO 1010
      END IF



! ******************************************************************
!                         END OF INTEGRATION
! ******************************************************************
 1000 CONTINUE



! NO ACCESS WRITE (6,'(5X,'' ** _RNSNEW_ ERROR DETECTION **'',////)')
! NO ACCESS WRITE (6,NAME1)
      IFAIL=-20
      RRNS=-1.
      RQOZ=-1.
      RETURN


!__ RADIUS and MASS OF THE ROTATING STAR

 2000 CONTINUE

      RQOZ=RMSPH
      RMNS=RQOZ*RMSUN*ufe30/5.61012E29
      TN1(MMM,8)=DMN
      TN1(MMM,9)=(DMN-RMSPH)/RMSPH

!__ CRITICAL ANGULAR MOMENTUM AT WHICH EQUATORIAL SHEDDING OCCURS  (1/SEC)
      CAVS=SQRT(RMSPH*DGRAV/RTH(0,2)) / RTH(0,2)
      TN1(MMM,6)=CAVS/(R0R0*UF10)





      WI=0.
      MOI=10
      IF (MOI /= 0) THEN
!............. MOMENT OF INERTIA OF THE ROTATING STAR ...............


      DO 606 I=0,NI
      DO 606 J=0,NJ
      IF (P(J,I) <= 0.) E(J,I)=0.
  606 CONTINUE

      TM=0.
      DO 600 J=1,NJP1
      DRN=FST251(J)
      DRN2=DRN*DRN
      DRN4=DRN2*DRN2
      HRK=FST251(J+1)-FST251(J)

      SUM=0.
      DO 700 I=1,NIP1
      DEN=E(J-1,I-1)
      IF (DEN == 0.) GO TO 700
      IF (I /= NIP1) THEN
      DELTH=T(I)-T(I-1)
      ELSE
      DELTH=T(I-1)-T(I-2)
      END IF
      SINT=SIN(T(I-1))
      SIN3T=SINT**3


!__ PERFORM INTEGRATION
!
!                        R             4
!                   TM = I [D R']  (R')  EPSILON(R'*R )
!                        0                           0
!
!      EASS='RNSNEW: CALL IGRN  20         '
!                                              **
!>>>>>>>  CALL IGRN(FST251,FST,IREM,X251,Y251,IREM,TM,C,NHU,10,1,AUX51,EASS)
      SUM=SUM + DELTH*SIN3T*DEN
  700 CONTINUE

      TM=TM + SUM * (HRK*DRN4)
  600 CONTINUE


!__ COMPUTE  LOG (I (g cm^2))
      WI= ALOG10(VPI*EC)+5.*ALOG10(R0R0)+ALOG10(TM)-ALOG10(FMEVFM)
     +   -52.

      END IF
!....................... END OF I-CALCULATION .........................




!__ CHANGE DIMENSIONS
      DO 2050 I=0,NI
      RTH(I,2)=RTH(I,2)*RMETER
 2050 CONTINUE


!__ WRITE RESULTS TO OUTPUT FILES
      IOP20=0
      IF (IOP20 /= 0) THEN

      OPEN (UNIT=20, FILE='DRQSN.dat', STATUS='NEW')
               OPEN (UNIT=30, FILE='plotting.top', STATUS='unknown')

      DO 2030 I=0,NI
      ANG=(180./PI)*T(I)
      WRITE (20,2031) I,ANG,RTH(I,2)
      RU=RTH(I,2)/EE03
      XV=RU*SIN(T(I))
      YV=RU*COS(T(I))
 2030 WRITE (30,2032) XV,YV
 2031 FORMAT ('  I=',I4,'  THETA=',F6.2,'  R(THETA)=',F9.2,' M')
 2032 FORMAT (E14.6,3X,E14.6)

      CLOSE ( 20 )
                   CLOSE ( 30 )
      END IF


      RETURN

      END

!SORTOR                 *
      SUBROUTINE SORTOR(X,N,IREF,IERSO)


! ----------------------------------------------------------------------
!
! PURPOSE:
!
!         IREF<0: CHEK ORDERING OF THE ELEMENTS OF ARRAY X (X(1)<X(2)<..
!                 ...<X(N)).
!         IREF=0: ORDER THE ELEMENTS OF ARRAY X (X(1)<X(2)<..<X(3)).
!         IREF>0: PERFORM BOTH OF THE ABOVE INTRUCTIONS.
!
!     INPUT:  X(N) (ELEMENTS OF THIS ARRAY ARE TO BE CHECKED FOR
!                   ORDERING, OR WILL BE ORDERED, DEPENDING ON  IREF)
!     OUTPUT: X(N) (UNCHANGED (IREF<0) RESPECTIVELY ORDERED (IREF=>0)
!
! ----------------------------------------------------------------------


      DIMENSION X(N)


      IERSO=0

      IF (IREF /= 0) THEN
      DO 50 K=1,N-1
      D=X(K+1)-X(K)

      IF (D > 0.) THEN
      GO TO 50
      ELSE
      IERSO=10
      GO TO 60
      END IF

   50 CONTINUE
      END IF

   60 IF (IREF < 0) RETURN

      IF (IREF > 0) THEN
      IF (IERSO == 0) GO TO 1000
      ELSE
      IERSO=0
      END IF


      I=1
  100 IF (I == N) GO TO 1000
      A=X(I)
      C=X(I+1)
      IF (C > A) THEN
      I=I+1
      ELSE IF (C < A) THEN
      X(I)=C
      X(I+1)=A
      I=1
      ELSE IF (C == A) THEN
      PROM=ABS(C)/1.E10
      IF (PROM == 0.) PROM=1.E-10
      C=C+PROM
      X(I+1)=C
      END IF
      GO TO 100


 1000 RETURN

      END

!fun
      function yfun (ttt)

! ----------------------------------------------------------------------
!
! purpose:
!          compute an approximate value of a function, which is
!          numerically specified at node ttt, via a cubic spline
!          approximation. the spline coefficients must be available
!          in matrix  cc.
!  note: interpolation  a n  d  extrapolation is permitted, i.e.,
!        ttt<x(1) and/or ttt>x(n).
!
! ----------------------------------------------------------------------

      parameter (kscc=10)
      common/spl/nn,cc(4,kscc),t(kscc)


      do 5 i=1,nn-1
      if (ttt  >=  t(i) .and. ttt  <=  t(i+1) ) go to 10
    5 continue

      if (ttt < t(1))then
      tx = t(1)
      sl = yfunp(tx)
      yfun = cc(1,1) + sl*(ttt-tx)
      return
      end if

      if (ttt > t(nn)) then
      tx = t(nn)
      sl = yfunp(tx)
      yfun = cc(1,nn) + sl*(ttt-tx)
      return
      end if

   10 x1 = t(i)
      c1 = cc(1,i)
      c2 = cc(2,i)
      c3 = cc(3,i)
      c4 = cc(4,i)

      dt = ttt-x1
      yfun = c1 + dt*(c2 + dt*(c3 + dt*c4))

      return

      end


!yfunp
      function yfunp(ttt)


! ----------------------------------------------------------------------
!
! purpose:
!          compute an approximate value of the first derivative of
!          a function, which is numerically specified at node ttt, by
!          a cubic spline approximation. the spline coefficients must be
!          available in matrix  cc.
!  note: interpolation  a n  d  extrapolation is permitted, i.e.,
!        ttt<x(1) and/or ttt>x(n).
!
! ----------------------------------------------------------------------

      parameter (kscc=10)
      common/spl/nn,cc(4,kscc),t(kscc)

      do 5 i=1,nn-1
      if (ttt  >=  t(i) .and. ttt  <=  t(i+1) ) go to 10
    5 continue

      if (ttt < t(1)) i = 1
      if (ttt > t(nn)) i = nn-1

   10 x1 = t(i)
      c2 = cc(2,i)
      c3 = cc(3,i)
      c4 = cc(4,i)

      dt = ttt-x1
      yfunp = c2 + dt*(c3*2 + dt*c4*3)

      return

      end

!sevu                                   ****
      subroutine sevu(x,y,nn,cc,ic,xout,yout,iout)


! ----------------------------------------------------------------------
!
! purpose:
!          evaluate the calculated spline at iout points (contained in
!          array xout)
!
! return: yout
!
! ----------------------------------------------------------------------


      real x(nn),y(nn),cc(4,ic),xout(iout),yout(iout)


      do 1 j=1,iout
      ttt = xout(j)



      do 5 i=1,nn-1

!__ search for the location of ttt (where the spline is to be calc.)
      ii = i
      if (ttt  >=  x(i) .and. ttt  <=  x(i+1) ) go to 10
    5 continue


      if ( ttt < x(1) .or. ttt > x(nn) ) then
         print*,' _sevu_ error call (argument out of range)'
         print*,' ttt=',ttt,' x(1)=',x(1),' x(nn)=',x(nn)
         stop
      end if


   10 x1 = x(ii)
      c1 = cc(1,ii)
      c2 = cc(2,ii)
      c3 = cc(3,ii)
      c4 = cc(4,ii)


      dt = ttt-x1
      yval = c1 + dt*(c2 + dt*(c3 + dt*c4))
      yout(j) = yval

    1 continue

      return

      end

!dsevu                                   ****
      subroutine dsevu(x,y,nn,cc,ic,xout,yout,iout)


! ----------------------------------------------------------------------
!
! purpose:
!          evaluate the first derivative of the cubic spline at iout
!          points
!
! return: yout
!
! ----------------------------------------------------------------------


      real x(nn),y(nn),cc(4,ic),xout(iout),yout(iout)


      do 1 j=1,iout
      ttt = xout(j)



      do 5 i=1,nn-1

!__ search for the location of ttt (where the spline is to be calc.)
      ii = i
      if (ttt  >=  x(i) .and. ttt  <=  x(i+1) ) go to 10
    5 continue


      if ( ttt < x(1) .or. ttt > x(nn) ) then
         print*,' _dsevu_ error call (argument out of range)'
         print*,' ttt=',ttt,' x(1)=',x(1),' x(nn)=',x(nn)
         stop
      end if


   10 x1 = x(ii)
      c2 = cc(2,ii)
      c3 = cc(3,ii)
      c4 = cc(4,ii)

      dt = ttt-x1
      ypval = c2 + dt*(c3*2 + dt*c4*3)
      yout(j) = ypval

    1 continue

      return

      end


!sevu8                                   ****
      subroutine sevu8(x,y,nn,cc,ic,xout,yout,iout)


! ----------------------------------------------------------------------
!
! purpose: double precision version of  sevu.
!          evaluate the calculated spline at iout points (contained in
!          array xout)
!
! return: yout
!
! ----------------------------------------------------------------------

      implicit real * 8 (a-h,o-z)

      real * 8 cc(4,ic)
      real * 8 x(nn),y(nn),xout(iout),yout(iout)


      do 1 j=1,iout
      ttt = xout(j)



      do 5 i=1,nn-1

!__ search for the location of ttt (where the spline is to be calc.)
      ii = i
      if (ttt  >=  x(i) .and. ttt  <=  x(i+1) ) go to 10
    5 continue


      if ( ttt < x(1) .or. ttt > x(nn) ) then
         print*,' _sevu8_ error call (argument out of range)'
         print*,' ttt=',ttt,' x(1)=',x(1),' x(nn)=',x(nn)
         stop
      end if


   10 x1 = x(ii)
      c1 = cc(1,ii)
      c2 = cc(2,ii)
      c3 = cc(3,ii)
      c4 = cc(4,ii)


      dt = ttt-x1
      yval = c1 + dt*(c2 + dt*(c3 + dt*c4))
      yout(j) = yval

    1 continue

      return

      end

!dsevu8                                   ****
      subroutine dsevu8(x,y,nn,cc,ic,xout,yout,iout)


! ----------------------------------------------------------------------
!
! purpose: double precision version of  dsevu.
!
!          evaluate the first derivative of the cubic spline at iout
!          points
!
! return: yout
!
! ----------------------------------------------------------------------

      implicit real * 8 (a-h,o-z)

      real * 8 cc(4,ic)
      real * 8 x(nn),y(nn),xout(iout),yout(iout)


      do 1 j=1,iout
      ttt = xout(j)



      do 5 i=1,nn-1

!__ search for the location of ttt (where the spline is to be calc.)
      ii = i
      if (ttt  >=  x(i) .and. ttt  <=  x(i+1) ) go to 10
    5 continue


      if ( ttt < x(1) .or. ttt > x(nn) ) then
         print*,' _dsevu8_ error call (argument out of range)'
         print*,' ttt=',ttt,' x(1)=',x(1),' x(nn)=',x(nn)
         stop
      end if


   10 x1 = x(ii)
      c2 = cc(2,ii)
      c3 = cc(3,ii)
      c4 = cc(4,ii)

      dt = ttt-x1
      ypval = c2 + dt*(c3*2 + dt*c4*3)
      yout(j) = ypval

    1 continue

      return

      end


!CSPLINE                       *
      SUBROUTINE CSPLINE(N,X,Y,C)


! ----------------------------------------------------------------------
!
! PURPOSE:
!          CSPLINE COMPUTES COEFFICIENTS C(4,N) ON CUBIC SPLINE FUNCTIONS.
!          CONDITIONS ON THE SPLINE FUNCTIONS ARE THAT THE SECOND
!          DERIVATIVE VANISHES ON THE BOUNDARIES.
!
! INPUT:  X(N), Y(N), N, SPLINE MATRIX  C(4,N)
!
! RETURN: COEFFICIENT-MATRIX C(J,I); I COUNTS THE INTERVALS WHERE THE
!         SPLINE FUNCTION IS COMPUTED. THE X VALUE IN EACH INTERVAL HAS
!         TO BE COUNTED FROM THE INTERVAL BORDER.
!
! ----------------------------------------------------------------------


      IMPLICIT REAL * 4 (A-H,O-Z)

      PARAMETER (NP100=100,NPVAR=3*NP100)

      REAL * 4 K(NPVAR), K1(NPVAR), K2(NPVAR), DX(NPVAR), DY(NPVAR)
      REAL * 4 C(4,N),X(N),Y(N)


       DO 10 I=1,N
       K (I) = 0.
       K1(I) = 0.
       K2(I) = 0.
       DO 10 J=1,4
       C(J,I) = 0.
   10  CONTINUE

       DO 11 I=1,N
   11  C(1,I) = Y(I)

      DX1 = X(2)-X(1)
      DY1 = Y(2)-Y(1)
      DX(1) = DX1
      DY(1) = DY1
      K2(2) = 1.
      K1(1) = 1.5*DY1/DX1
      K2(1) = -0.5

!  Since we have not enough conditions on one border
!  to determine K(1) and K(2), we have to go through all intervls
!  to get the value of K(2)
!  at the other border. Therefore we need two values (first) for each
!  coefficient, one the 'constant' K1, the other the factor in front
!  of K(2) (K2(i,j)).

      N1 = N-1

      DO 12 I=3,N
      DX2 = DX1
      DY2 = DY1
      DX1 = X(I) - X(I-1)
      DY1 = Y(I) - Y(I-1)
      DX(I-1) = DX1
      DY(I-1) = DY1
      DXX  = DX1/DX2
      K1(I) = 3.*(DY1/DX1+DY2/DX2*DXX) -
     +           2.*(1.+DXX)*K1(I-1) - DXX*K1(I-2)
      K2(I) = - 2.*(1.+DXX)*K2(I-1) - DXX*K2(I-2)
   12 CONTINUE



! Now we can determine K(2) = D0
      D0 = (3.*DY1/DX1-K1(N1)-2.*K1(N)) / (K2(N1)+2.*K2(N))
      DO 13 I=1,N
   13 K(I) = K1(I) + D0*K2(I)

!  Determination of the coefficients themselves
      DO 14 I=1,N1
      C(2,I) = K(I)
      C(3,I) = (3.*DY(I)/DX(I) - 2.*K(I) - K(I+1)) / DX(I)
      C(4,I) = (K(I)+K(I+1)-2.*DY(I)/DX(I)) / DX(I) / DX(I)
   14 CONTINUE


!     write (6,3000) d0
!3000 format (' D0 ',5(1x,e12.4)/5(1x,e12.4))
!     write (6,3001) (x(i),y(i),i=1,n)
!3001 format (' X und Y'/2(2x,e12.4))
!     do 99 j=1,n1
!  99 write (6,3002) (c(i,j),i=1,4)
!3002 format (4(1x,e12.4))

      RETURN

      END

!CSPL8                       *
      SUBROUTINE CSPL8(N,X,Y,C)


! ----------------------------------------------------------------------
!
! PURPOSE: DOUBLE PRECISION VERSION OF  CSPLINE
!
!          CSPL8 COMPUTES COEFFICIENTS C(4,N) ON CUBIC SPLINE FUNCTIONS.
!          CONDITIONS ON THE SPLINE FUNCTIONS ARE THAT THE SECOND
!          DERIVATIVE VANISHES ON THE BOUNDARIES.
!
! INPUT:  X(N), Y(N), N, SPLINE MATRIX  C(4,N)
!
! RETURN: COEFFICIENT-MATRIX C(J,I); I COUNTS THE INTERVALS WHERE THE
!         SPLINE FUNCTION IS COMPUTED. THE X VALUE IN EACH INTERVAL HAS
!         TO BE COUNTED FROM THE INTERVAL BORDER.
!
! ----------------------------------------------------------------------


      IMPLICIT REAL * 8 (A-H,O-Z)

      PARAMETER (NP100=100,NPVAR=3*NP100)

      REAL * 8 K(NPVAR), K1(NPVAR), K2(NPVAR), DX(NPVAR), DY(NPVAR)
      REAL * 8 C(4,N)
      REAL * 8 X(N),Y(N)


       DO 10 I=1,N
       K (I) = 0.
       K1(I) = 0.
       K2(I) = 0.
       DO 10 J=1,4
       C(J,I) = 0.
   10  CONTINUE

       DO 11 I=1,N
   11  C(1,I) = Y(I)

      DX1 = X(2)-X(1)
      DY1 = Y(2)-Y(1)
      DX(1) = DX1
      DY(1) = DY1
      K2(2) = 1.
      K1(1) = 1.5*DY1/DX1
      K2(1) = -0.5

!  Since we have not enough conditions on one border
!  to determine K(1) and K(2), we have to go through all intervls
!  to get the value of K(2)
!  at the other border. Therefore we need two values (first) for each
!  coefficient, one the 'constant' K1, the other the factor in front
!  of K(2) (K2(i,j)).

      N1 = N-1

      DO 12 I=3,N
      DX2 = DX1
      DY2 = DY1
      DX1 = X(I) - X(I-1)
      DY1 = Y(I) - Y(I-1)
      DX(I-1) = DX1
      DY(I-1) = DY1
      DXX  = DX1/DX2
      K1(I) = 3.*(DY1/DX1+DY2/DX2*DXX) -
     +           2.*(1.+DXX)*K1(I-1) - DXX*K1(I-2)
      K2(I) = - 2.*(1.+DXX)*K2(I-1) - DXX*K2(I-2)
   12 CONTINUE



! Now we can determine K(2) = D0
      D0 = (3.*DY1/DX1-K1(N1)-2.*K1(N)) / (K2(N1)+2.*K2(N))
      DO 13 I=1,N
   13 K(I) = K1(I) + D0*K2(I)

!  Determination of the coefficients themselves
      DO 14 I=1,N1
      C(2,I) = K(I)
      C(3,I) = (3.*DY(I)/DX(I) - 2.*K(I) - K(I+1)) / DX(I)
      C(4,I) = (K(I)+K(I+1)-2.*DY(I)/DX(I)) / DX(I) / DX(I)
   14 CONTINUE


!     write (6,3000) d0
!3000 format (' D0 ',5(1x,d12.4)/5(1x,d12.4))
!     write (6,3001) (x(i),y(i),i=1,n)
!3001 format (' X und Y'/2(2x,e12.4))
!     do 99 j=1,n1
!  99 write (6,3002) (c(i,j),i=1,4)
!3002 format (4(1x,d12.4))

      RETURN

      END



!selint                              ***
      subroutine selint(x,y,iy,n1,n2,erg,p,iopt)

! ----------------------------------------------------------------------
!
!__ purpose:
!
!            select the next kp7 grid points, contained in x,
!            which are closest to node  p;
!            save the selected points in workx (x values) and  worky
!            (y-values) which are then used as grid points to
!            interpolate the functions  y(.,x).
!
! ----------------------------------------------------------------------

      parameter (kp3=2,kp4=kp3+1,kp7=kp3+kp4,kpiy=1)
      real x(n2),y(iy,n2),erg(n1),workx(kp7),worky(kpiy,kp7)



      k3 = kp3
      k4 = kp4
      k7 = kp7

      if (iy > kpiy) then
         ikpiy = kpiy
         write (6,32)
   32    format (//,' Warning:  _selint_ --> fatal error detection ',//)
         write (6,42) iy,ikpiy
   42    format (//,' iy=',i3,' ikpiy=',i4,//)
      stop 'error stop'
      end if


      if (n2 <= 2) then

         write (6,32)
         write (6,22) n1,n2,iy,k3,k4
   22    format(//,' n1=',i3,' n2=',i3,' iy=',i3,' k3=',i3,' k4=',i3,//)
         stop 'error stop'


      else if ( (n2 > 2).and.(n2 <= (k3+k4)) ) then


!__ interpolation, restriction to the next kp7 closest grid points is
!   not necessary

      do 222 liy=1,iy
      do 222 k=1,n2
      workx(k) = x(k)
      worky(liy,k) = y(liy,k)
  222 continue
!                                          ***
      call inter(iopt,workx,worky,iy,n1,n2,erg,p)


      else


!__ search for the grid point which is closest to  p  (determine
!   integer number  nmin1)
      abw = abs(x(1)-p)
      do 350 k=1,n2
      diff = abs(x(k)-p)
      if (diff <= abw) then
         nmin1 = k
         abw = diff
      end if
  350 continue

      if (nmin1 <= k3) then
         nmin1 = k4
      else if (nmin1 >= (n2-k3)) then
         nmin1 = n2 - k3
      end if

      ir1 = nmin1 - k3
      ir2 = nmin1 + k3


!__ interpolation

      do 555 liy=1,iy
      do 555 k=ir1,ir2
      workx(k-(nmin1-k4)) = x(k)
      worky(liy,k-(nmin1-k4)) = y(liy,k)
  555 continue




!__ interpolate function  y(.,x) at node  p
!                                          ***
      call inter(iopt,workx,worky,iy,n1,k7,erg,p)

!      do 4 i=1,k3+k4
!      print*,' i=',i,' wkx=',workx(i),' wky=',worky(i),' p=',p,
!     +' erg=',erg(1)
!    4 continue


      end if

      return

      end


!selspl                              **
      subroutine selspl(x,y,nxy,c,ic,yp,p,iref,eass)

! ----------------------------------------------------------------------
!
! purpose:
!
!          compute an interpolatory value to the function  y(x) at the
!          point  x=p by means of a cubic spline interpolation.
!          the definition of  iref  is given in subroutine spline.
!
! return:  yp=y(x=p)
!
! ----------------------------------------------------------------------


      real x(nxy),y(nxy),c(4,ic)

      parameter (kp3=2,kp4=kp3+1,kp7=kp3+kp4,kpxxyy=1)
      dimension workx(kp7),worky(kp7),xxx(kpxxyy),yyy(kpxxyy)


      character * 30 eass


      ixxyy= kpxxyy
      k3   = kp3
      k4   = kp4
      k7   = kp7


      if (ic > 10) then
         write (6,32) ic
   32    format (///,' Warning: _selspl_  ic > 10, ic=', i4, //)
      end if

      xxx(ixxyy)=p
      yyy(ixxyy)=0.



      if (nxy <= 2) then

         write (6,22) nxy,ic,k3,k4
   22    format(//,' nxy=',i3,' ic=',i3,' k3=',i3,' k4=',i3,//)
         stop 'error stop in  _selspl_'

      else if ( (nxy > 2).and.(nxy <= (k3+k4)) ) then


!__ perform interpolation; restriction to the next kp7 closest
!   grid points is not necessary

      do 222 k=1,nxy
      workx(k) = x(k)
      worky(k) = y(k)
  222 continue

!                                           ***
      call spline(workx,worky,nxy,c,nxy,xxx,yyy,ixxyy,iref,eass)


      else


!__ search for the grid point which is closest to  p  (determine
!   integer number  nmin1)
      abw = abs(x(1)-p)
      do 350 k=1,nxy
      diff = abs(x(k)-p)
      if (diff <= abw) then
         nmin1 = k
         abw = diff
      end if
  350 continue


      if (nmin1 <= k3) then
         nmin1 = k4
      else if (nmin1 >= (nxy-k3)) then
         nmin1 = nxy - k3
      end if

      ir1 = nmin1 - k3
      ir2 = nmin1 + k3


!>>>> start  interpolation

      do 555 k=ir1,ir2
      workx(k-(nmin1-k4)) = x(k)
      worky(k-(nmin1-k4)) = y(k)
  555 continue


!__ interpolate function  y(x)  at node  p
!                                         ***
      call spline(workx,worky,k7,c,k7,xxx,yyy,ixxyy,iref,eass)

!      do 4 i=1,k3+k4
!      print*,' i=',i,' wkx=',workx(i),' wky=',worky(i),' p=',p,
!     +' y(p)=',yyy(ixxyy)
!    4 continue

      end if

      yp=yyy(ixxyy)


      return

      end


!selspl8                              **
      subroutine selspl8(x,y,nxy,c,ic,yp,p,iref,eass)

! ----------------------------------------------------------------------
!
! purpose: double precision version of selspl
!
!          compute an interpolatory value to the function  y(x) at the
!          point  x=p by means of a cubic spline interpolation
!          the definition of  iref  is given in subroutine spline.
!
! return:  yp=y(x=p)
!
! ----------------------------------------------------------------------


      implicit real * 8 (a-h,o-z)

      real * 8 c(4,ic)
      real * 4  x(nxy),y(nxy),yp,p

      parameter (kp3=2,kp4=kp3+1,kp7=kp3+kp4,kpxxyy=1)
      real * 8 workx(kp7),worky(kp7),xxx(kpxxyy),yyy(kpxxyy)

      character * 30 eass


      ixxyy= kpxxyy
      k3   = kp3
      k4   = kp4
      k7   = kp7


      if (ic > 10) then
         write (6,32) ic
   32    format (///,' Warning:  _selspl8_  ic > 10, ic=', i4, //)
      end if

      xxx(ixxyy)=p
      yyy(ixxyy)=0.



      if (nxy <= 2) then

         write (6,22) nxy,ic,k3,k4
   22    format(//,' nxy=',i3,' ic=',i3,' k3=',i3,' k4=',i3,//)
         stop 'error stop in  _selspl8_'

      else if ( (nxy > 2).and.(nxy <= (k3+k4)) ) then


!__ perform interpolation; restriction to the next kp7 closest
!   grid points is not necessary

      do 222 k=1,nxy
          dpx = x(k)
          dpy = y(k)
      workx(k) = dpx
      worky(k) = dpy
  222 continue
!                                            ***
      call spline8(workx,worky,nxy,c,nxy,xxx,yyy,ixxyy,iref,eass)


      else


!__ search for the grid point which is closest to  p  (determine
!   integer number  nmin1)
      abw = abs(x(1)-p)
      do 350 k=1,nxy
      diff = abs(x(k)-p)
      if (diff <= abw) then
         nmin1 = k
         abw   = diff
      end if
  350 continue


      if (nmin1 <= k3) then
         nmin1 = k4
      else if (nmin1 >= (nxy-k3)) then
         nmin1 = nxy - k3
      end if

      ir1 = nmin1 - k3
      ir2 = nmin1 + k3


!>>>> start  interpolation

      do 555 k=ir1,ir2
      dpx     = x(k)
      dpy     = y(k)
      workx(k-(nmin1-k4)) = dpx
      worky(k-(nmin1-k4)) = dpy
  555 continue


!__ interpolate function  y(x) at node  p
!                                          ***
      call spline8(workx,worky,k7,c,k7,xxx,yyy,ixxyy,iref,eass)

!      do 4 i=1,k3+k4
!      print*,' i=',i,' wkx=',workx(i),' wky=',worky(i),' p=',p,
!     +' y(p)=',yyy(ixxyy)
!    4 continue

      end if


      yp=yyy(ixxyy)


      return

      end

!zero22                                           *****
      subroutine zero22(x,y,iy,n,xwk,ywk,zwkl,nwk,xzero)


! ----------------------------------------------------------------------
!
! purpose:
!          zero22 calculates the zeroes of iy functions, contained in  y
!
! return: xzero(iy) containes the zeroes of y(.,x)
!
! ----------------------------------------------------------------------


      real x(n),y(iy,n),xwk(nwk),ywk(iy,nwk),zwkl(nwk),xzero(iy),erg(1)

      common/MiscParameters/lim007,icro22,iine9,dzer9l,dzer9p


      mm = nwk

      a   = x(1)
      b   = x(n)
      del = (b-a)/float(mm-1)


!__ compute x-grid
      do 10 i=1,mm
      xwk(i) = a+del*float(i-1)
   10 continue


!__ compute functions  y(.,x) at grid points  x=xwk
      do 20 i=1,mm
      r      = xwk(i)
      iopt   = 1
!                           ***
      call selint(x,y,1,1,n,erg,r,iopt)
      do 25   j=1,iy
   25 ywk(j,i) = erg(j)
   20 continue


!__ test of monotoneous behavior of the interpolated functions y(.,x)
      istop = 1

      do 40  j=1,iy
      do 45  i=1,mm
   45 zwkl(i) = ywk(j,i)
!                     ****
!>>>>>    call sortor(zwkl,mm,-10,ierso)
      ierso=0
!>>>>>
      if (ierso /= 0) then
      write (6,'(//'' _zero22_ elements of ywk not ordered''/)')
      do 41 k=1,mm
   41 write (6,411) j,k,xwk(k),ywk(j,k)
  411 format (' j=',i3,' k=',i3,'  xwk(k)=',e13.6,'  ywk(j,k)=',e13.6)
      istop = -1
      end if
   40 continue
      if (istop < 0) stop 'zero22'


!__   search for zero-points of functions y(.,x)

      do 50 j=1,iy

      do 55 i=1,mm
      iz = i
      if (ywk(j,i) < 0.) go to 9999
      if (i == mm) then
      write (6,'(//'' _zero22_ search for zeroes failed''//)')
      do 51 k=1,mm
   51 write (6,411) j,k,xwk(k),ywk(j,k)
      write (6,'(//'' write input data''/)')
      do 57 k=1,n
   57 write (6,58) k,x(k),y(j,k)
   58 format (2x,'k=',i5,2x,'x=',e12.6,2x,'y=',e12.6)
      istop = -1
      end if
   55 continue
 9999 continue
      if (istop < 0) stop 'zero22'

      xzero(j) = (xwk(iz)+xwk(iz-1))/2.

   50 continue



! **********************************************************************
!__ test result
      do 60 j=1,iy
      if ( (xzero(j) >= x(n-1)) .and. (xzero(j) <= x(n)) ) go to 7777
      if (icro22 <= lim007) then
      write (6,'(//'' _zero22_ use arithmethic average!!'',/,
     +       '' j='',i3,''  iy='',i3,//)') j,iy
      end if
      xzero(j) = (x(n)+x(n-1))/2.
   60 continue
      icro22 = icro22 + 1
! **********************************************************************


!>>>>>>
      print*,' xzero=',xzero(1)
!>>>>>>

 7777 return

      end


!SHOW
      SUBROUTINE SHOW(IPREST,ILL,IUL,ICAB,ICO21A,ICO21B)


! ----------------------------------------------------------------------
!     PURPOSE:
!              CREATE OUTPUT FILE 'prns_showStatus.dat' WHICH CONTAINS
!              INFORMATION ABOUT THE PRESENT STATUS OF THE CALCULATION
! ----------------------------------------------------------------------

      COMMON/CTYPE/TYPEOS,TEOSCB

      CHARACTER * 50 TYPEOS
      CHARACTER * 45 TEOSCB



!__ write information about the calculation's present status to output 
!   file  prns_showStatus.dat
      OPEN (UNIT=17, FILE='prns_showStatus.dat', STATUS='unknown')      

      WRITE (17,'(//)')
      WRITE (17,'(1X,A50)')  TYPEOS
      WRITE (17,'(1X,A45)')  TEOSCB

      WRITE  (17,1) IPREST,ILL,IUL
    1 FORMAT (//,1X,'I_present:',1X,I3,/,1X,'I_lower_limit:',1X,I3,/,
     +        1X,'I_upper_limit:',1X,I3)

      WRITE  (17,2) ICAB,ICO21A,ICO21B
    2 FORMAT (//,1X,'I_present:',1X,I3,/,1X,'I_co21a:',1X,I3,/,
     +        1X,'I_co21b:',1X,I3)



      RETURN

              END


!copy
      subroutine copy(x_in,y_in,n,x_out,y_out)

! ----------------------------------------------------------------------
!
! purpose:
! copy elements contained in x_in, y_in to x_out,y_out  in reverse 
! order
!
! ----------------------------------------------------------------------

      dimension x_in(n),y_in(n),x_out(n),y_out(n)

      do 100 i=1,n
      x_out(n+1-i) = x_in(i)
      y_out(n+1-i) = y_in(i)
 100  continue

!      do 200 i=1,ng
! 200  print*,' x_out=',x_out(i),' y_out=',y_out(i)

      return
             end


!profiles
      subroutine profiles(rrns,rmqoz,irem,ec,igravBary)

! Purpose: write  | r | m(r) | epsilon(r) |  n(r) |  p(r)  to output file

      real * 8 ee14,ee34,ee44,ee18,eem18
      real * 8 ee03,ee13,ee15
      
      parameter (m651=160501)
      common/cnsta1/f251(m651),fout(m651),fst(m651),fst251(m651),feng1
     +              (m651),x251(m651),y251(m651),f251pm(m651),fdrip(m6
     +              51)
      common/bn01/axbn1(m651),axbn2(m651),pdenex(m651),einthx(m651),
     +            deindp(m651),dehdrh(m651),pdent(m651)
      common/RadialScalingFactor/r0r0
      common/PowersTen/ee03,ee15,ee18,eem18,ee13,ee14,ee34,ee44
      common/ctype/typeos,teoscb
      
      real engden(m651),baryden(m651)

      character * 50 typeos
      character * 45 teoscb
      open (unit=17, file='mass_energy_density_profiles.dat', 
!     +                                  access='append', status='old') 
     +                                                status='unknown') 


      do i=1,irem
         engden(i)  = feng1(i)*ec
         baryden(i) = pdenex(i)
      end do

      write (17,'(/,'' jobname = '',a50)') typeos
!      write (17,'( '' igravBary='', I2)') igravBary
      
      if (igravBary == 1) then
         write (17,'( '' igravBary='', I2, '' ( >< stellar model '',
     + ''with a given gravitational mass)'')') igravBary
      else if (igravBary == 2) then
         write (17,'( '' igravBary='', I2, '' ( >< stellar model '',
     + ''with a given baryon mass)'')') igravBary
      else if (igravBary == 3) then
         write (17,'( '' igravBary='', I2, '' ( >< stellar model '',
     + ''with a given baryon number)'')') igravBary
      else
         write (17, '(/, '' Value of igravBary not properly specified'',
     +   /)')
      end if
      
     
      write (17, 199) f251(irem), fst251(irem)*r0r0/ee18, ec, fout(1)*ec
 199  format(1x,'M=', F7.4, ' M_sol,', 2x, 'R=', F8.4, ' km,', 2x,'e_c='
     +         , F8.2, ' MeV/fm^3,', 2x, 'P_c=', F8.2, ' MeV/fm^3') 
      write (17,'(/, ''  radial distance  |     m(r)       |  energy'',
     +               '' density  | number density  |    pressure'')')
      write (17,'(6x, ''(meter)'', 10x, ''(M_sol)'', 12x,''(MeV/fm^3)'',
     +            7x,''(1/fm^3)'', 10x, ''(Mev/fm^3)'')')
      write (17,'('' ----------------------------------------------'',
     +            ''-----------------------------------------'')')
      
      do i=1,irem, 10
      write (17, 10) fst251(i)*r0r0/ee15,f251(i),engden(i),baryden(i),
     +               fout(i)*ec               
      end do
      
 10   format (2x,e14.7,4x,e14.7,4x,e14.7,4x,e14.7,4x,e14.7)
!      write (17, '(''  R='',e14.7,'' m'',4x,''m/m_sun='',e14.7)') 
!     +       rrns,rmqoz
      write (17, '(''   --------- at surface of the star -----------'',
     +             ''-----------------------'')') 

      close (17)

      return 

             end

! last line of file
