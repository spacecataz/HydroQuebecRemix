     ISEE-3 ONE-MINUTE AVERAGED MAGNETIC FIELD DATA (1984-1990)
       COVERING THE HELIOSPHERIC PHASE OF THE ISEE 3 MISSION

Title: 1min_hgi_1984_1990
NSSDC ID: SPHE-00844
Old ID: 78-079A-02U

This data set contains 1-minute magnetic field data in simple ASCII records. 
It was created at NSSDC from a more complex, multi-resolution data set
(NSSDC ID = SPHE-00673; Old ID = 78-079A-02D) provided by the Principal 
Investigator team and now available from
ftp://nssdcftp.gsfc.nasa.gov/spacecraft_data/isee/isee3/magnetic_fields/1min_ascii/ 

For years 1984-1990 we added spacecraft position in HGI coordinate (see note2).

Note that a similar data set, covering the Magnetospheric phase of the ISEE 3 mission,
is also available: 1min_1978_1983, NSSDC ID: SPHE-00843, Old ID: 78-0709A-02T.

              Data format for years 1984-1990


 1   IYR                        I4            Last 2 digits of year
 2   IDAY                       I4            Day of year (Jan. 1 = Day 1)
 3   HOUR                       I3            Hour of day (0-23) 
 4   MIN                        I3            Min of hour  (0-59)
 5  <Bx>                       F6.1           SE,nT 
 6  <By>                       F6.1           SE,nT
 7  <Bz>                       F6.1           SE,nT
 8  <Bx^2>                     F8.1 
 9  <BxBy>                     F8.1 
10  <BxBz>                     F8.1
11  <By^2>                     F8.1
12  <ByBz>                     F8.1
13  <Bz^2>                     F8.1
14  <cos alpha> = <Bx/|B|>     F7.3
15  <cos beta>  = <By/|B|>     F7.3
16  <cos gamma> = <Bz/|B|>     F7.3
17  <|B|>                      F6.1            nT
18  <|B|^2>                    F8.1 
19   X                        F10.2          GSE, Re
20   Y                        F10.2          GSE, Re
21   Z                        F10.2          GSE, Re
22   R                         F7.2          HGI, Au
23   latitude                  F6.1          HGI, deg.
24   longitude                 F7.1          HGI, deg.

----------------------------------------------------

Note1:
Values of <Bx>, <By>, <Bz>, and <|B|> are given in nanoteslas. The
coordinate system for the B-field components is the JPL-defined I,S
coordinate system (origin at the spacecraft): I is the unit vector in
the direction of the ISEE-3 spin axis (positive in the northward
direction), and S is the unit vector from the spacecraft to the sun.
The z-axis is parallel to to I, the y-axis to the cross-product I x S,
and the x-axis to Y x Z. The I,S coordinate system is approximately
the same as the Solar Ecliptic (SE) system since the spacecraft z-axis
(the spin axis) is maintained within 0.5 degree of perpendicular to
the ecliptic plane. (SE is defined the same way as GSE, but with the
spacecraft [point of observation] substituted for Earth). 

*Note2:
2 Heliographic Inertial Coordinate System (HGI)
The HGI coordinates are Sun-centered and inertially fixed with respect to an 
X-axis directed along the intersection line of theecliptic and solar equatorial
 planes, and defines zero of the longitude, HGI_LONG. The solar equator plane 
is inclined at 7.25degrees from the ecliptic. This direction was towards ecliptic
 longitude of 74.367 deg on 1 January 1900 at 12:00 UT; because of
the precession of the Earth's equator, this longitude increases by 1.4 deg/century.
 The Z-axis is directed perpendicular to and
northward of the solar equator, and the Y-axis completes the right-handed set. 
The longitude, HGI_LONG increase from zero in
the X-direction towards Y-direction.The latitude HG_LAT increases to +90 deg towards
 the north pole, and to -90 deg towardsm south pole.

Updated:" 8/14/02 HKH
