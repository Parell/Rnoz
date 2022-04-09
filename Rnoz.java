/*
                Interactive Rocket Thrust Program

     Program to perform one dimensional design and analysis of
                        rocket nozzle
                Derived from turbine nozzle program

                     Version 1.5b   - 13 Dec 05

                         Written by Tom Benson
                       NASA Glenn Research Center

>                              NOTICE
>This software is in the Public Domain.  It may be freely copied and used in
>non-commercial products, assuming proper credit to the author is given.  IT
>MAY NOT BE RESOLD.  If you want to use the software for commercial
>products, contact the author.
>No copyright is claimed in the United States under Title 17, U. S. Code.
>This software is provided "as is" without any warranty of any kind, either
>express, implied, or statutory, including, but not limited to, any warranty
>that the software will conform to specifications, any implied warranties of
>merchantability, fitness for a particular purpose, and freedom from
>infringement, and any warranty that the documentation will conform to the
>program, or any warranty that the software will be error free.
>In no event shall NASA be liable for any damages, including, but not
>limited to direct, indirect, special or consequential damages, arising out
>of, resulting from, or in any way connected with this software, whether or
>not based on warranty, contract, tort or otherwise, whether or not injury
>was sustained by persons or property or otherwise, and whether or not loss
>was sustained from, or arose out of the results of, or use of, the software
>or services provided hereunder.

      New Test : 
                 * change layout slightly
                   add altitude input to flow


                                                     TJB 13 Dec 05
*/

import java.awt.*;
import java.lang.Math ;

public class Rnoz extends java.applet.Applet {

   final double convdr = 3.14515926/180.;
   final double runiv = 1545. ;

   static double rthrt, rexit, rzero ;
   static double press, temp ;
   static double ath, pt, pamb, pexit, tt, machth, mflow ;
   static double athroat, athmx, athmn ;
   static double arat, aexmx, aexmn ;
   static double azrat, azmx, azmn, azero ;
   static double ptin, ptmx, ptmn ;
   static double psin, pemx, pemn ;
   static double alt, altin, altmx, altmn ;
   static double ttin, ttmx, ttmn ;
   static double lngth, lngmx, lngmn ;
   static double fact ;
   static double npr,uex,fgros ;
   static double mweight,rgas,gam0,gamma,tcomb,ofrat,mwtab,ttab ;
   static double mexit,aexit,ans,rns,lns,mnsup ;
   static int counter,flomode ;
   static int fuelopt, fuelold, gamopt, tcomopt, molopt, oxopt, ofopt ;

   static double xg[][]  = new double[11][45] ;
   static double yg[][]  = new double[11][45] ;

   static int antim,ancol,lunits ; 
   static double aconv,lconv,fconv,pconv,tref,tconv,mconv1 ;

   int xt,yt,sldloc,mode ;


   public void setDefaults() {

       lunits = 0 ;
       aconv = 1.;                         /*  area sq inches    */
       lconv = 1.;                         /*  length feet    */
       fconv = 1.0;                        /* pounds   */
       pconv = 1.0  ;                   /* lb/sq in */
       tref =  459.7 ;                 /* zero rankine */
       tconv = 1.0 ;                    /* degrees F */
       mconv1 = 1.0 ;                   /* airflow rate lbs/sec*/

       athmx = 300. ;
       athmn = .1 ;
       aexmx = 100. ;
       aexmn = 1. ;
       azmx = 10. ;
       azmn = 1.1 ;
       ptmx = 3000. ;
       ptmn = 1. ;
       pemx = 15. ;
       pemn = 0.0 ;
       ttmx = 6500. ;
       ttmn = 500. ;
       altmn = 0.0 ;
       altmx = 100000. ;
       altin = 0.0 ;
       alt = 0.0 ;
       flomode = 0 ;

       mwtab = mweight = 16.0 ;
       rgas = 1716. ;                /* air - ft2/sec2 R */
       molopt = 1 ;
       fuelold = fuelopt = 4 ;
       oxopt = 1 ;
       ttab = tcomb = ttin = 5870. ;
       gamopt = 1 ;
       gam0 = 1.32 ;
       gamma = 1.22 ;
       tcomopt = 1 ;
       ofopt = 0;
       ofrat = 8.0 ;

       athroat  = 150.0 ;
       arat = 50.0 ;
       azero  = 600.0 ;
       azrat = 4.0 ;
       ptin = 2000. ;
       psin = 14.7 ;

       mode = 0 ;
       ans = 0.0 ;
       fact = .2 ;
       xt = 40;  yt = 0; 
       sldloc = 40 ;
       lngth = .5 ;
       lngmx = (athroat/144.) * 20.0;
       lngmn = (athroat/144.) * .01;

   }

   public void computeAtmos() {

       if (alt <= 36152.) {           // Troposphere
          temp = 518.6 - 3.56 * alt/1000. ;
          press = 2116. * Math.pow(temp/518.6,5.256) ;
       }
       if (alt >= 36152. && alt <= 82345.) {   // Stratosphere
          temp = 389.98 ;
          press = 2116. * .2236 *
             Math.exp((36000.-alt)/(53.35*389.98)) ;
       }
       if (alt >= 82345.) {
          temp = 389.98 + 1.645 * (alt-82345)/1000. ;
          press = 2116. *.02456 * Math.pow(temp/389.98,-11.388) ;
       }

       psin = pconv * press / 144.  ;
   }


   public void comPute() {
       double gm1,fac1,trat,aircor,g0 ;
       double psub,psup ;
       double mnsub,ptrat,anso,ansn,pso,psn,deriv ;
       double gamtfac ;

       g0    = 32.2 ;
       rgas = runiv * g0 / mweight ;
       alt = altin / lconv ;
       ath = athroat / aconv ;
       rthrt = Math.sqrt(ath / 3.1415926) ;
       aexit = arat * ath ;
       rexit = Math.sqrt(aexit / 3.1415926) ;
       azero = azrat * ath ;
       rzero = Math.sqrt(azero / 3.1415926) ;
       pt    = ptin / pconv ;
       
       if (flomode == 1) computeAtmos () ;

       pamb  = psin / pconv ;
       tt    = (ttin + tref) / tconv ;
       if (gamopt == 1) {
           gamtfac = getGama(tt,gamopt) ;
           gamma = gam0 * gamtfac / 1.4 ;
       }
       gm1   = gamma - 1.0 ;
       fac1  = gm1 / gamma ;

       counter = 0 ;
       machth = 1.0 ;         // assume flow is choked
       aircor = getMfr(1.0,rgas,gamma) ;
       mexit = getMach(2,(aircor/arat),rgas,gamma) ;
       psup = getPisen(mexit,gamma) * pt ;
       mflow = aircor*ath*(pt/14.7)/Math.sqrt(tt/518.) ;

  // under expanded nozzle
       if (pamb <= psup) {   
          mode = 0 ;
          trat = getTisen(mexit,gamma) ;
          uex = mexit * Math.sqrt (gamma * rgas * trat * tt) ;
          pexit = psup ;
       }

  // over expanded nozzle
       if (pamb > psup) {  
           // find exit pressure at which normal shock leaves the nozzle
          mnsup = mexit ;
          psub = psup * getNsps(mnsup,gamma) ;

     // slightly overexpanded .. no shock in nozzle
          if (pamb <= psub) {   
             mode = 1 ;
             pexit = psup ;
             trat = getTisen(mexit,gamma) ;
             uex = mexit * Math.sqrt (gamma * rgas * trat * tt) ;
          }

     // highly overexpanded .. normal shock in nozzle
          if (pamb > psub) {   
             mode = 2 ;
             pexit = pamb ;
             anso  = aexit ;
             mnsup = mexit ;
             mnsub = getNsmac(mnsup,gamma) ;
             ptrat = getNspt(mnsup,gamma) ;
             pso = getPisen(mnsub,gamma) * ptrat * pt ;
             ansn = anso - 1. ;
             while ((Math.abs(pexit - pso) > .001) && (counter < 20)) {
                ++ counter ;
                mnsup = getMach(2,(aircor*ath/ansn),rgas,gamma) ;
                mnsub = getNsmac(mnsup,gamma) ;
                ptrat = getNspt(mnsup,gamma) ;
                mexit = getMach(0,(aircor/arat/ptrat),rgas,gamma) ;
                psn = getPisen(mexit,gamma) * ptrat * pt ;
                deriv = (psn-pso)/(ansn-anso) ;
                pso = psn ;
                anso = ansn ;
                ansn = anso + (pexit -pso)/deriv ;
             }
             ans = anso ;
             rns = Math.sqrt(ans / 3.1415926) ;
             trat = getTisen(mexit,gamma) ;
             uex = mexit * Math.sqrt (gamma * rgas * trat * tt) ;
          }
       }

       if (pamb > .0001) npr = pt / pamb ;
       else npr = 1000.;

       fgros = mflow * uex / g0 + (pexit - pamb) * aexit ;

       loadOut() ;
 
       loadGeom() ;

       view.repaint();
    }
 
    public double getGama(double temp, int opt) {
             // Utility to get gamma as a function of temp
       double number,a,b,c,d ;
       a =  -7.6942651e-13;
       b =  1.3764661e-08;
       c =  -7.8185709e-05;
       d =  1.436914;
       if(opt == 0) {
          number = 1.4 ;
       }
       else {
          number = a*temp*temp*temp + b*temp*temp + c*temp +d ;
       }
       return(number) ;
    } 

    public double getMfr(double mach, double gascon, double gam) {
/* Utility to get the corrected weightflow per area given the Mach number */
       double number,fac1,fac2;

       fac2 = (gam+1.0)/(2.0*(gam-1.0)) ;
       fac1 = Math.pow((1.0+.5*(gam-1.0)*mach*mach),fac2);
       number =  20.78 * Math.sqrt(gam/gascon) * mach / fac1 ;

       return(number) ;
    }

     public double getMach (int sub, double corair, double gascon, double gam) {
/* Utility to get the Mach number given the corrected airflow per area */
         double number,chokair,a,b,k;     /* iterate for mach number */
         double deriv,machn,macho,airo,airn;
         int iter ;

         a = (gam -1)/2.0 ;
         b = -(gam+1.0)/(2.0*(gam-1.0)) ;
         k = 20.78 * Math.sqrt(gam/gascon) ;

         chokair = getMfr(1.0,gascon,gam) ;
         if (corair > chokair) {
           number = 1.0 ;
           return (number) ;
         }
         else {
           if (sub == 1) macho = 1.0 ;   /* sonic */
           else {
              if (sub == 2) macho = 1.703 ; /* supersonic */
              else macho = .25;                /* subsonic */
              airo = getMfr(macho,gascon,gam) ;    /* initial guess */
              iter = 1 ;
              machn = macho - .2  ;
              while (Math.abs(corair - airo) > .0001 && iter < 50) {
                 airn =  getMfr(machn,gascon,gam) ;
                 deriv = (airn-airo)/(machn-macho) ;
                 airo = airn ;
                 macho = machn ;
                 machn = macho + (corair - airo)/deriv ;
                 ++ iter ;
              }
           }
           number = macho ;
         }
         return(number) ;
     }

     public double getPisen(double machin, double gam)  {
       /* Utility to get the isentropic pressure ratio given the mach number */
       double number,fac1,gm1,mach1s;
       mach1s = machin*machin ;
       gm1 = gam - 1.0 ;
       fac1 = 1.0 + .5*gm1*mach1s;
       number = Math.pow(1.0/fac1,gam/gm1) ;

       return(number) ;
    }

    public double getTisen(double machin, double gam)  {
      /* Utility to get the isentropic temperature ratio given the mach number*/       double number,gm1,mach1s;
       mach1s = machin*machin ;
       gm1 = gam - 1.0 ;
       number = 1.0 /(1.0 + .5*gm1*mach1s) ;

       return(number) ;
    }

    public double getNsps (double machin, double gam) {
          // NACA 1135 - normal shock relation ps ratio
       double number, gm1, gp1, msq, fac1, fac2 ;

       msq = machin * machin ;
       gm1 = gam - 1.0 ;
       gp1 = gam + 1.0 ;

       fac2 = (2.0*gam*msq - gm1)/gp1 ;
       number = fac2 ;

       return(number) ;
    }

    public double getNspt (double machin, double gam) {
          // NACA 1135 - normal shock relation pt ratio
       double number, gm1, gp1, msq, fac1, fac2 ;

       msq = machin * machin ;
       gm1 = gam - 1.0 ;
       gp1 = gam + 1.0 ;

       fac2 = (2.0*gam*msq - gm1)/gp1 ;
       fac1 = (gp1*msq)/(gm1*msq + 2.0) ;
       number = (Math.pow(fac1,(gam/gm1)))
               * (Math.pow((1.0/fac2),(1.0/gm1))) ;

       return(number) ;
    }

    public double getNsmac (double machin, double gam) {
          // NACA 1135 - normal shock relation  mach
       double number, gm1, gp1, msq, fac1, fac2 ;

       msq = machin * machin ;
       gm1 = gam - 1.0 ;
       gp1 = gam + 1.0 ;

       fac2 = (2.0*gam*msq - gm1)/gp1 ;
       fac1 = (gp1*msq)/(gm1*msq + 2.0) ;
       number = Math.sqrt(msq / (fac2 * fac1)) ;

       return(number) ;
   }
 
   public void setUnits() {   // Switching Units
       double aths,athms,pts,tts,ttmxs,pss,ptmxs,ttcos,ttabs ;
       double altns, altmxs ;
       int i1,i2,i3,i4 ;
 
       aths = athroat / aconv ;
       athms = athmx / aconv ;
       pts  = ptin / pconv ;
       tts  = ttin / tconv ;
       ttmxs  = ttmx / tconv ;
       ttcos  = tcomb / tconv ;
       ttabs  = ttab / tconv ;
       pss  = psin / pconv ;
       ptmxs  = ptmx / pconv ;
       altns = altin / lconv ;
       altmxs = altmx / lconv ; 

       switch (lunits) {
          case 0: {                             /* English */
            aconv = 1.;                         /*  sq in   */
            lconv = 1.;                         /*  feet    */
            fconv = 1.0;                        /* pounds   */
            pconv = 1.0  ;                   /* lb/sq in */
            tref =  459.7 ;                 /* zero rankine */
            tconv = 1.0 ;                    /* degrees F */
            mconv1 = 1.0 ;                   /* airflow rate lbs/sec*/
            pn.in.geo.l.li1.setText("Ath sq in") ;
            pn.in.flo.l.li1.setText("Pto  psi") ;
            pn.in.flo.l.li2.setText("Tto  F") ;
            if (flomode == 0) pn.in.flo.l.li3.setText("Pfs  psi") ;
            if (flomode == 1) pn.in.flo.l.li3.setText("Alt  ft") ;
            pn.in.chem.r.l6.setText("F") ;
            pn.out.l2.setText("Uex fps") ;
            pn.out.l12.setText("Ath sq in") ;
            pn.out.l14.setText("Pt0 psi") ;
            break;
          }
          case 1: {                             /* Metric */
            aconv = 6.4516;                      /* sq cm */
            lconv = .3048;                      /* meters */
            fconv = 4.448 ;                     /* newtons */
            pconv = 6.891 ;               /* kilo-pascals */
            tref = 273.1 ;                /* zero kelvin */
            tconv = 0.555555 ;             /* degrees C */
            mconv1 = .4536 ;                /* kg/sec */
            pn.in.geo.l.li1.setText("Ath sq cm") ;
            pn.in.flo.l.li1.setText("Pto  kPa") ;
            pn.in.flo.l.li2.setText("Tto  C") ;
            if (flomode == 0) pn.in.flo.l.li3.setText("Pfs  kPa") ;
            if (flomode == 1) pn.in.flo.l.li3.setText("Alt  m") ;
            pn.in.chem.r.l6.setText("C") ;
            pn.out.l2.setText("Uex mps") ;
            pn.out.l12.setText("Ath sq cm") ;
            pn.out.l14.setText("Pt0 kPa") ;
            break ;
          }
       }
 
       athroat = aths * aconv ;
       athmx = athms * aconv ;
       ptin = pts * pconv ;
       ttin = tts * tconv ;
       ttmx = ttmxs * tconv ;
       tcomb = ttcos * tconv ;
       ttab = ttabs * tconv ;
       psin = pss * pconv ;
       ptmx = ptmxs * pconv ;
       altin = altns * lconv ;
       altmx = altmxs * lconv ;
       athmn = .01 * aconv ;
       ptmn = 1. * pconv ;
       pemx = 15. * pconv ;
       pemn = 0.0 * pconv ;
       ttmn = 500. * tconv ;

       pn.in.geo.l.f1.setText(String.valueOf((float) athroat)) ;
       pn.in.flo.l.f1.setText(String.valueOf((float) ptin)) ;
       pn.in.flo.l.f2.setText(String.valueOf((float) ttin)) ;
       if (flomode ==  0) pn.in.flo.l.f3.setText(String.valueOf(filter3(psin))) ;
       if (flomode ==  1) pn.in.flo.l.f3.setText(String.valueOf(filter0(altin))) ;
       pn.in.chem.r.f6.setText(String.valueOf(filter0(tcomb))) ;

       i1 = (int) (((athroat - athmn)/(athmx-athmn))*1000.) ;
       i2 = (int) (((ptin - ptmn)/(ptmx-ptmn))*1000.) ;
       i3 = (int) (((ttin - ttmn)/(ttmx-ttmn))*1000.) ;
       i4 = (int) (((psin - pemn)/(pemx-pemn))*1000.) ;
       if (flomode == 1) i4 = (int) (((altin - altmn)/(altmx-altmn))*1000.) ;

       pn.in.geo.r.s1.setValue(i1) ;
       pn.in.flo.r.s1.setValue(i2) ;
       pn.in.flo.r.s2.setValue(i3) ;
       pn.in.flo.r.s3.setValue(i4) ;

       return ;
    }

    public void loadGeom() {
       double x0,y0,x01,y01,xth,yth,xth1,yth1,xex,yex,xex1,yex1 ;
       double yns,xns ;
       double factor ;
       int i,j ;

       factor = athroat/aconv ;

   // rocket nozzle geometry
       y0   = factor * 1.0 ;                  
       x0   = factor * rzero ;  
       y01  = factor * 20.0 * .5 ;   
       x01  = factor * rzero ;
       yth  = factor * 60.0 * .5 ;   
       xth  = factor * rthrt ;
       yth1 = factor * 80.0 * .5 ;  
       xth1 = factor * rthrt ;
       yex  = factor * 180.0 * (.5 + lngth) ;  
       xex  = factor * rexit ;
       yex1 = factor * 200.0 * (.5 + lngth) ; 
       xex1 = factor * rexit ;

   // shock in nozzle location
       if (mode == 2) {
          xns = factor * rns ;
          lns = lngth * rns / rexit ;
          yns = factor * 180.0 * (.5 + lns) ;
          xg[0][44] = xns ;
          yg[0][44] = yns ;
       }

       for(i=0; i<=5; ++ i) {   //  basic geometry
           xg[0][i] = i * (x01 - x0)/5.0 + x0 ;
           yg[0][i] = i * (y01 - y0)/5.0 + y0 ;
       } 
       for(i=6; i<=14; ++ i) {
           xg[0][i] = (i-5) * (xth - x01)/9.0 + x01 ;
           yg[0][i] = (i-5) * (yth - y01)/9.0 + y01 ;
       } 
       for(i=15; i<=18; ++ i) {
           xg[0][i] = (i-14) * (xth1 - xth)/4.0 + xth ;
           yg[0][i] = (i-14) * (yth1 - yth)/4.0 + yth ;
       } 
       for(i=19; i<=28; ++ i) {
           xg[0][i] = (i-18) * (xex - xth1)/10.0 + xth1 ;
           yg[0][i] = (i-18) * (yex - yth1)/10.0 + yth1 ;
       } 
       for(i=29; i<=34; ++ i) {
           xg[0][i] = (i-28) * (xex1 - xex)/3.0 + xex ;
           yg[0][i] = (i-28) * (yex1 - yex)/3.0 + yex ;
       } 

       for (j=1; j<=5; ++ j) {  // lower flow
           for(i=0; i<=34; ++ i) {
              xg[j][i] = (1.0 - .2 * j) * xg[0][i] ;
              yg[j][i] = yg[0][i] ;
           }
       }

       for (j=6; j<=10; ++ j) {  // mirror
           for(i=0; i<=34; ++ i) {
              xg[j][i] = -xg[10-j][i] ;
              yg[j][i] = yg[0][i] ;
           }
       }
   }
            // change fuel options
                          else {
                             fuelold = fuelopt ;
                             if (fuelopt == 0) {  // air
                                oxopt = 0 ;
                                mwtab = 29.0 ;
                                if (gamopt == 1) gam0 = 1.4 ;
                                ttab = 4000. * tconv ;
                             }
                             if (fuelopt == 1) {  // solid - aluminum
                                oxopt = 0 ;
                                mwtab = 29.4 ;
                                if (gamopt == 1) gam0 = 1.24 ;
                                ttab = 5846. * tconv ;
                             }
                             if (fuelopt == 2) {  // hydrogen peroxide
                                oxopt = 0 ;
                                mwtab = 22.2 ;
                                if (gamopt == 1) gam0 = 1.32 ;
                                ttab = 1450. * tconv ;
                             }
                             if (fuelopt == 3) {  // hydrazine
                                oxopt = 0 ;
                                mwtab = 32.1 ;
                                if (gamopt == 1) gam0 = 1.42 ;
                                ttab = 1125. * tconv ;
                             }
                             if (fuelopt == 4) {  // H2
                                    // default oxidizer is O2 
                                oxopt = 1 ;
                                mwtab = 16.0 ;
                                if (gamopt == 1) gam0 = 1.32 ;
                                ttab = 5870. * tconv ;
                                ofrat = 8.0 ;
                             }
                             if (fuelopt == 5) {  // JP
                                    // default oxidizer is O2
                                oxopt = 1 ;
                                mwtab = 22.0 ;
                                if (gamopt == 1) gam0 = 1.34 ;
                                ttab = 5770. * tconv ;
                                ofrat = 2.3 ;
                             }
                             if (fuelopt == 6) {  // Kerosene
                                    // default oxidizer is O2
                                oxopt = 1 ;
                                mwtab = 22.0 ;
                                if (gamopt == 1) gam0 = 1.34 ;
                                ttab = 5825. * tconv ;
                                ofrat = 2.3 ;
                             }
                             if (fuelopt == 7) {  // UDMH
                                    // default oxidizer is RFNA
                                oxopt = 5 ;
                                mwtab = 22.0 ;
                                if (gamopt == 1) gam0 = 1.34 ;
                                ttab = 5200. * tconv ;
                                ofrat = 2.6 ;
                             }
        // set default mweights, gam0 , tcomb
                      if (oxopt == 1) {       // oxygen
                         if (fuelopt == 3) { // hydrazine
                             mwtab = 18.0 ;
                             if (gamopt == 1) gam0 = 1.35 ;
                             ttab = 5370. * tconv ;
                             ofrat = .75 ;
                         }
                         if (fuelopt == 4) { // hydrogen
                              // default ofrat 8.0
                             ofopt = 0 ;
                             o.y.ofch.select(ofopt) ;
                             mwtab = 16.0 ;
                             if (gamopt == 1) gam0 = 1.32 ;
                             ttab = 5870. * tconv ;
                             ofrat = 8.0 ;
                         }
                         if (fuelopt == 5) { // JP-4
                             mwtab = 22.0 ;
                             if (gamopt == 1) gam0 = 1.34 ;
                             ttab = 5770. * tconv ;
                             ofrat = 2.3 ;
                         }
                         if (fuelopt == 6) { // Kerosene
                             mwtab = 22.0 ;
                             if (gamopt == 1) gam0 = 1.34 ;
                             ttab = 5825. * tconv ;
                             ofrat = 2.3 ;
                         }
                      }

                      if (oxopt == 2) {       // flourine
                         if (fuelopt == 4) { // hydrogen
                              // default ofrat 9.42
                             ofopt = 0 ;
                             o.y.ofch.select(ofopt) ;
                             mwtab = 10.0 ;
                             if (gamopt == 1) gam0 = 1.42 ;
                             ttab = 8100. * tconv ;
                             ofrat = 9.42 ;
                         }
                      }

                      if (oxopt == 3) {       // hydrogen peroxide
                         if (fuelopt == 3) { // hydrazine
                             mwtab = 19.0 ;
                             if (gamopt == 1) gam0 = 1.32 ;
                             ttab = 4690. * tconv ;
                             ofrat = 1.7 ;
                         }
                         if (fuelopt == 5) { // JP-4
                             mwtab = 22.0 ;
                             if (gamopt == 1) gam0 = 1.3 ;
                             ttab = 5770. * tconv ;
                             ofrat = 2.3 ;
                         }
                      }

                      if (oxopt == 4) {       // nitrous oxide
                         if (fuelopt == 6) { // Kerosene
                             mwtab = 25.3 ;
                             if (gamopt == 1) gam0 = 1.34 ;
                             ttab = 5461. * tconv ;
                             ofrat = 6.5 ;
                         }
                      }

                      if (oxopt == 5) {       // RFNA
                         if (fuelopt == 7) { // UDMH
                             mwtab = 22.0 ;
                             if (gamopt == 1) gam0 = 1.34 ;
                             ttab = 5200. * tconv ;
                             ofrat = 2.6 ;
                         }
                      }
        
                                  // oxidizer/fuel ratio  options
                                 ofopt  = ofch.getSelectedIndex() ;

                                 if (oxopt == 1) {  // hydrogen - oxygen
                                   if (ofopt == 0) {   // ofrat = 8
                                      mwtab = 16.0 ;
                                      if (gamopt == 1) gam0 = 1.32 ;
                                      ttab = 5870. * tconv ;
                                      ofrat = 8.0 ;
                                   }
                                   if (ofopt == 1) {   // ofrat = 6
                                      mwtab = 13.0 ;
                                      if (gamopt == 1) gam0 = 1.32 ;
                                      ttab = 5500. * tconv ;
                                      ofrat = 6.0 ;
                                   }
                                   if (ofopt == 2) {   // ofrat = 4
                                      mwtab = 10.0 ;
                                      if (gamopt == 1) gam0 = 1.32 ;
                                      ttab = 5000. * tconv ;
                                      ofrat = 4.0 ;
                                   }
                                 }
                                 if (oxopt == 2) {  // hydrogen - flourine
                                   if (ofopt == 0) {   // ofrat = 9.4
                                      mwtab = 10.0 ;
                                      if (gamopt == 1) gam0 = 1.42 ;
                                      ttab = 8100. * tconv ;
                                      ofrat = 9.42 ;
                                   }
                                   if (ofopt == 1) {   // ofrat = 3.8
                                      mwtab = 7.8 ;
                                      if (gamopt == 1) gam0 = 1.42 ;
                                      ttab = 4600. * tconv ;
                                      ofrat = 3.8 ;
                                   }
                                 }
        
    
     public void start() {
        if (runner == null) {
           runner = new Thread(this) ;
           runner.start() ;
        }
        antim = 0 ;                
        ancol = 1 ;                        
     }

     public void run() {
       int timer ;

       timer = 100 ;
       while (true) {
          ++ antim ;
          try { Thread.sleep(timer); }
          catch (InterruptedException e) {}
          view.repaint() ;
          if (antim == 3) {
             antim = 0;
             ancol = - ancol ;       
          }
       }
     }