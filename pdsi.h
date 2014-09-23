/*
This file was used to calculate the PDSI values online in
the Dynamic Global Vegetation Model (DGVM) LPJ-GUESS for
the publication Steinkamp & Hickler (accepted for 2015).
Is drought-induced forest dieback globally increasing?
Journal of Ecology.

Here a short description on the usage. If something is
unclear do not hesitate to contact me (joerg.steinkamp[at]senckenberg.de)

This file must be included in the guess.h file in the framework directory

Additional changes to guess.h:
The class Patch must have a class pdsi of class PDSI

in driver.cpp in the directory modules
in the routine dailyaccounting_patch:

before "if (date.day==0) {"
  PDSI& pdsi=patch.pdsi;
after "if (date.day==0) {"
  if (date.year==0) pdsi.reset(soil.soiltype.awc);

after "if (date.dayofmonth==0) { [...] }"
pdsi.update(date.month, patch.mintercep[date.month], patch.mevap[date.month], patch.maet[date.month], patch.mpet[date.month], patch.mrunoff[date.month], soil.wcont, patch.stand.climate.prec);

additionally in the used output routine the monthly pdsi
value has to be added use the get_X() method of the PDSI class for that.

 */

class PDSI {
// CAFEC: Climatically appropriate for existing conditions
 private:
  /*
    The formulas are based on Palmer (1965) and Alley (1984)

    The Parameter alpha, beta, gamma, delta in Palmer (1965) [here a, b, c, d]
    are calculated as long term mean values starting with 50% water filled pore 
    space in the first month of the first spinup/simulation year.

    Modifications:
    - the unit is mm instead of inches
    - precipitation is added to the potential runoff
    - instead of applying a spatially averaged K value (17.67 in Palmer, 1965)
      in the calculation of the new monthly K values I applied a running long 
      term annual mean variable called "sdk_avg".
  */
  int nyears; // weighting factor for antecedent years
  int nyears_max; // maximum number of antecedent years for weighting factor

  double a[12];
  double b[12];
  double c[12];
  double d[12];
  int na[12];
  int nb[12];
  int nc[12];
  int nd[12];

  // long term monthly mean values
  double D[12];  // precipitation difference
  double PE[12]; // potential evapotranspiration
  double R[12];  // Recharge
  double RO[12]; // Runoff
  double L[12];  // Loss
  double P[12];  // Precipitation

  // actual monthly PDSI values for one year
  double K[12];
  double Kw[12];
  double X[12];
  double Z[12];

  // long term annual running sum of D*K
  double sdk_avg;

  double p;

  double awc[2];
  double wc[2][2];   // first dim upp(0) and low(1); second dim: new(0) and old(1)

  // DEBUGGING
  double SDK[12];

 public:
  // debugging values for offline pdsi calculation
  double wcontupp[12];
  double wcontlow[12];

  /************
   Constructor
  *************/
  PDSI() {awc[0] = -1;}

  PDSI(double *awcin) {
    awc[0]   = awcin[0];
    awc[1]   = awcin[1];
    wc[0][1] = 0.5 * awc[0];
    wc[1][1] = 0.5 * awc[1];
    p        = 0.0;
    sdk_avg  = 17.67;
    
    nyears     = 0;
    nyears_max = 100;

    for(int m=0; m<12; m++) {
      a[m]  = b[m]  = c[m]  = d[m]  = 0.0;
      na[m] = nb[m] = nc[m] = nd[m] = 0;
      Kw[m] = K[m]  = X[m]  = Z[m]  = 0.0;
      D[m]  = PE[m] = R[m]  = RO[m] = L[m] = P[m] = 0.0;
    }
  }
  // Destructor
  //~PDSI() { }
  /********************************************
   Function to accumulate precipitation and
   calculate monthly PDSI values.

   Call this function on each day of the month
  *********************************************/
  void update(int m, double mintercep, double mevap, double maet, double mpet, double mrunoff, double *wcont, double dprec) {
    double et;
    double pe;
    double r;
    double pr;
    double ro;
    double pro;
    double l;
    double pl;
    double pl_upp;
    double pl_low;

    double p_cafec;
    double dp;

    double sdk;
    double t;

    if (date.dayofmonth==0) {
      wcontupp[m] = wcont[0];
      wcontlow[m] = wcont[1];
    }

    // Accumulate daily values and return if it is not the last day of a month
    if (!date.islastday) {
      // runoff is taken from LPJ directly
      // recharge and loss are calculated as the difference between the
      // water content at the start and end of a month

      //precipitation
      p = p + dprec;
      return;
    }
    p  = p + dprec;

    // Take some values from LPJ (monthly PET, AET, runnoff and actual watercontent)
    pe = mpet;
    et = mintercep + mevap + maet;
    ro = mrunoff;
    wc[0][0] = wcont[0] * awc[0];
    wc[1][0] = wcont[1] * awc[1];

    // Calculate all missing values, which are not already in LPJ
    // recharge
    r  = max(0., (wc[0][0] + wc[1][0]) - (wc[0][1] + wc[1][1]));

    // loss
    l = max(0., (wc[0][1] + wc[1][1]) - (wc[0][0] + wc[1][0]));

    // potential recharge
    pr = awc[0] + awc[1] - (wc[0][1] + wc[1][1]);

    // potential runoff;
    /* added also the precipitation to potential runoff. Otherwise
       the relation r/pr can grow to much too high numbers.*/
    pro = wc[0][1] + wc[1][1] + p;
    //pro = wc[0][1] + wc[1][1];
    //pro = max(0., 3*P[m] - pr);

    //potential loss
    //pl_upp = min(pe, wc[0][1]);
    //pl_low = min(wc[1][1], (pe - pl_upp) * wc[1][1] / (awc[0] + awc[1]));
    //pl = pl_upp + pl_low;

    pl = wc[0][1] + wc[1][1];
    // if necessary values can be converted to inches here

    // Precipitation to maintain lonterm monthly mean values 
    // and difference to actual precipitation
    p_cafec = a[m]*pe + b[m]*pr + c[m]*pro - d[m]*pl;
    dp      = p - p_cafec;

    // this should not be done in the first year.
    // !!! DIVISION BY 0 !!!
    if (P[m] + L[m] > 0) {
      t = (PE[m] + R[m] + RO[m]) / (P[m] + L[m]);
      // substracted 2.107252 (log10(25.4)*1.5) due to the conversion from inches to mm (1in = 25.4mm)
      // R script used to calculate is is attached at the end of this file
      K[m] = 1.5 * log10((t + 2.8) / D[m] * 25.4) +0.5; //- 1.607252;

      if (K[m] < 0.0 && date.year>1000) dprintf("JS_DEBUG_K: %02i %04i %f\n",m+1, date.year, K[m]);

      sdk = 0.0;
      for (int i=0; i<12; i++) {
        sdk = sdk + D[i] * K[i];
      }
      if (sdk < 0.0 && date.year>1000) dprintf("JS_DEBUG_SDK: %02i %04i %f\n",m+1, date.year, sdk);

      // DEBUGGING
      SDK[m] = sdk;

      Kw[m] = fabs(sdk_avg) * K[m] / sdk;

      Z[m] = Kw[m] * dp;

      if (fabs(Z[m]) >= 100.0) {
        dprintf("PDSI: Reset Z: %02i.%04i %f\n", m, date.year, Z[m]);
        Z[m] = Z[m] / fabs(Z[m]) * 99.99;
      }

      if (date.year>1000 && fabs(Z[m]) > 25) {
        dprintf("JS_DEBUG start %i\n", m);
        dprintf("p_cafec = %f * %f + %f * %f + %f * %f - %f * %f\n", a[m], pe, b[m], pr, c[m], pro, d[m], pl);
        dprintf("%f = %f - %f\n", dp, p, p_cafec);
        dprintf("%f = %f * %f\n", Z[m], Kw[m], dp);
        dprintf("%f = abs(%f) * %f / %f\n", Kw[m], sdk_avg, K[m], sdk);
        for (int i=0; i<12; i++) {
          dprintf("   %f = %f * %f\n",D[i] * K[i], D[i], K[i]);
        }

	dprintf("JS_DEBUG end\n");
        //fail("JS_DEBUG end\n");
      }


      if (m==0) {
        X[m] = Z[m]/3.0 + 0.897*X[11];
      } else {
        X[m] = Z[m]/3.0 + 0.897*X[m-1];
      }
    } // END: if (P[m] + L[m]>0)

    // Update the long term mean values !!! check for division by 0 !!!
    if (isnan(et/pe) || isinf(et/pe) || et/pe>=1 ) {
      a[m] = (a[m]*na[m] + 1.0) / (na[m] + 1.0);
      if (na[m] < nyears_max) na[m]++;
    } else {
      a[m] = (a[m]*na[m] + et/pe) / (na[m] + 1.0);
      if (na[m] < nyears_max) na[m]++;
    }

    if (r/pr<=1.0) {
      b[m] = (b[m]*nb[m] + r/pr) / (nb[m] + 1.0);
      if (nb[m] < nyears_max) nb[m]++;
    }
    if (ro/pro<=1.0) {
      c[m] = (c[m]*nc[m] + ro/pro) / (nc[m] + 1.0);
      if (nc[m] < nyears_max) nc[m]++;
    }
    if (l/pl<=1.0) {
      d[m] = (d[m]*nd[m] + l/pl) / (nd[m] + 1.0);
      if (nd[m] < nyears_max) nd[m]++;
    }

    // above if - else can be replaced by min function

    D[m]  = (D[m]  * nyears + fabs(dp)) / (nyears + 1.0);
    PE[m] = (PE[m] * nyears + pe)       / (nyears + 1.0);
    R[m]  = (R[m]  * nyears + r)        / (nyears + 1.0);
    RO[m] = (RO[m] * nyears + ro)       / (nyears + 1.0);
    L[m]  = (L[m]  * nyears + l)        / (nyears + 1.0);
    P[m]  = (P[m]  * nyears + p)        / (nyears + 1.0);

    // !!! this is a crude workaround to avoid division by 0 !!!
    if (D[m] == 0.0) D[m]=0.1;

    // set the water content to the new value and reset the precipitation.
    // This happens on each last day of the month (no check needed, see return above)
    wc[0][1] = wc[0][0];
    wc[1][1] = wc[1][0];
    p = 0.0;
    // increase the number of years "to keep in mind"
    if (m==11 && nyears<nyears_max) nyears++;
  } // END: void update (...)

  /************************************
    Functions to return monthly values
  *************************************/
  double get_X(int m) {
    return(X[m]);
  }

  double get_Z(int m) {
    return(Z[m]);
  }

  double get_K(int m) {
    return(K[m]);
  }
  // DEBUGGING
  double get_Kw(int m) {
    return(Kw[m]);
  }
  double get_SDK(int m) {
    return(SDK[m]);
  }

  /*******************************************
    Reset global values for each new location
  ********************************************/
  void reset(double *awcin) {
    awc[0]   = awcin[0];
    awc[1]   = awcin[1];
    wc[0][1] = 0.5 * awc[0];
    wc[1][1] = 0.5 * awc[0];
    p        = 0.0;
    sdk_avg  = 17.67;

    nyears = 0;
    nyears_max = 100;

    for(int m=0; m<12; m++) {
      a[m]  = b[m]  = c[m]  = d[m]  = 0.0;
      na[m] = nb[m] = nc[m] = nd[m] = 0;
      Kw[m] = K[m]  = X[m]  = Z[m]  = 0.0;
      D[m]  = PE[m] = R[m]  = RO[m] = L[m] = P[m] = 0.0;
    }
  }
};
/* LITERATURE:

   Alley, W. M., 1984: The Palmer drought severity index: Limitations and assumptions. 
     J. Climate Appl. Meteor., 23, 1100–1109.

   Palmer, W. C., 1965: Meteorological drought. 
     U.S. Weather Bureau Research Paper 45, Washington, DC, 58 pp. 
*/

