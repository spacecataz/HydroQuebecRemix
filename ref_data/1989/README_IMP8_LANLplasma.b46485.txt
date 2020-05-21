Subject: Data Format Description for the LANL 2-min Plasma Data on IMPs 6,7,8
From: J. H. King
Organization: National Space Science Data Center

Entered as B46485 in NSSDC's Technical Reference File, January 12, 1998.
Modified April 7, 2003, HKH.

=============================================================================
IMP 6,7,8 LANL19
2-MINUTE PLASMA DATA
Data items are represented in ASCII.
=============================================================================
   ITEM              FORMAT     DESCRIPTION            UNITS/NOTES
   ====              ======     ===========            =====
1  Space Craft ID    A4         = "IMP8"
   space             X1
2  Date              I6         YYMMDD
   space             X1
3  Seconds of day    F8.1
   space             X1
4  Utime             F8.1       HHMMSS.0               See Footnote below
   space             X1
5  Radius            E11.4                             Re
   space             X1
6  SELAT             E11.4                             deg
   space             X1
7  SELONG            E11.4                             deg
   space             X1
8  Density           E11.4      Proton Numerical       cm-3
   space             X1
9  Flow Speed        E11.4      Proton Numerical       km/sec
   space             X1
10 Flow Azimuth      E11.4      Proton Numerical       deg
   space             X1
11 TPAR (MAX)        E11.4      Proton Numerical       K
   space             X1
12 TPER (MIN)        E11.4      Proton Numerical       K
   space             X1
13 PSIT              E11.4      Pressure Tensor        deg
                                Symmetry Axis, AZ
   space             X1
14 ALPHA/PROTON      E11.4      Flux Ratio
   space             X1
15 Flow Speed AP     E11.4      Relative Flow Speed    km/sec
   space             X1
16 Flow Azimuth AP   E11.4      Relative Flow Angle    deg
   space             X1
17 TRatio            E11.4      Alpha Temp/Proton Temp
   space             X1
18 Temp              E11.4      Alpha Temp Anisotropy
   space             X1
19 PSIT              E11.4      Alpha Pressure         deg
                                Asymmetry Axis

Footnote: the MM values of this word are incorrect for IMP 8 for 1973-1985 
and will be corrected shortly (April 7, 2003). In the meantime, note that 
word 3 gives correct time tag information.
=============================================================================
