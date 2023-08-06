# HITRAN Application Programming Interface (HAPI)
===============================================

Current version: 1.2.2.1

## Version history

  1) FIXED GRID BUG (ver. 1.1.0.1)
  2) FIXED OUTPUT FORMAT FOR CROSS-SECTIONS (ver. 1.1.0.1)
  3) ADDED CPF BY SCHREIER (JQSRT_112_2011) (ver. 1.1.0.2)
  4) OPTIMIZED EXPRESSION EVALUATIONS FOR SELECT (ver. 1.1.0.3)
  5) ADDED SUPPORT FOR MIXTURES (ver. 1.1.0.4)
  6) ADDED SUPPORT FOR USER-DEFINED ENV DEPENDENCES (ver. 1.1.0.5)
  7) ADDED PROFILE SELECTION (ALPHA) (ver. 1.1.0.6)
  8) ADDED METADATA FOR HTP, FIXED NORMALIZATION IN CONVOLVESPECTRUMSAME (ver. 1.1.0.7)
  9) FIXED A "LONELY HEADER" BUG IN CACHE2STORAGE (ver. 1.1.0.7.1)
  10) ADDED SUPPORT FOR PHOSGENE AND CYANOGEN (ver. 1.1.0.7.2)
  11) OPTIMIZED STORAGE2CACHE (by Nils-Holger Loeber) (ver. 1.1.0.7.3)
  12) ADDED SKIPPABLE PARAMETERS IN HEADERS (ver. 1.1.0.7.4)
  13) ADDED SUPPORT FOR FORTRAN D-NOTATION (ver. 1.1.0.7.5)
  14) ADDED SUPPORT FOR WEIRD-FORMATTED INTENSITY VALUES E.G. "2.700-164" (ver. 1.1.0.7.6)
  15) ADDED TIPS-2017 (ver. 1.1.0.8)
  16) ADDED SUPPORT FOR CUSTOM EXTENSIONS OF THE DATA FILES (ver. 1.1.0.8.1)
  17) FIXED LINK TO (2,0) ISOTOPOLOGUE IN TIPS-2017 (ver. 1.1.0.8.2)
  18) ADDED SAVEHEADER FUNCTION (ver. 1.1.0.8.3)
  19) ADDED METADATA FOR SF6 (ver. 1.1.0.8.4)
  20) ADDED D2O ISOTOPOLOGUE OF WATER TO DESCRIPTION (ver. 1.1.0.8.5)
  21) FIXED LINE ENDINGS IN STORAGE2CACHE AND QUERYHITRAN (ver. 1.1.0.8.6)
  22) ADDED SUPPORT FOR NON-INTEGER LOCAL ISO IDS (ver. 1.1.0.8.7)
  23) FIXED PARAMETER NAME CASE BUG (by Robert J. Hargreaves) (ver. 1.1.0.8.8)
  24) CAST LOCAL_ISO_ID=0 TO 10 FOR CARBON DIOXIDE (ver. 1.1.0.8.9)
  25) USING NUMPY.ARRAYS FOR NUMERIC COLUMNS OF LOCAL_TABLE_CACHE (ver. 1.1.0.9.0)
  26) ADDED DESCRIPTIONS FOR BROADENING BY H2O (ver. 1.1.0.9.1)
  27) ADDED PROXY SUPPORT IN FETCH AND FETCH_BY_IDS (ver. 1.1.0.9.2)
  28) ADDED LIMIT FOR NUMBER OF LINES DURING TABLE READ (ver. 1.1.0.9.3)
  29) FIXED ABSOLUTE PATH BUG IN TABLE NAMES (ver. 1.1.0.9.4)
  30) CORRECTED ABUNDANCE OF THE HD ISOTOPOLOGUE (ver. 1.1.0.9.5)
  31) ADDED UNIFIED INTERFACES FOR ABSCOEF AND XSC CALCULATIONS (ver. 1.1.0.9.6)
  32) ADDED PARLISTS FOR LINE MIXING (VOIGT AND SDVOIGT) (ver. 1.1.0.9.7)
  33) ADDED SUPPORT FOR ROSENKRANZ LM PARAMETERS TO PCQSDHC AND LORENTZ (ver. 1.1.1.0)
  34) FIXED THE TYPEERROR IN ARANGE (ver. 1.1.2.0)
  35) ADDED NEW FUNCTIONAL INTERFACES FOR ALL CROSS-SECTION CALCULATING ROUTINES (ver. 1.2.0.0)
  36) ADDED CALCULATION OF THE ISO_ID TABLE ON STARTUP (ver. 1.2.1.0)
  37) ADDED SUPPORT FOR TIPS-2021 (ver. 1.2.2.0)
  38) FIXED BUG WITH WAVENUMBERGRID (ver. 1.2.2.1)

## Introduction

The HITRAN Application Programming Interface (HAPI) [1] is a set of routines in Python which aims to provide remote access to functionality and data given by the HITRANonline. At the present time, the API can download, filter and process line-by-line transition data.

The main purpose of this API is to extend the functionality of the main site, in particular, in the calculation of spectra using several types of line shape, including the flexible HT (Hartmann-Tran) profile [2] and optionally accounting for instrumental functions. Each feature of the API is represented by a Python function taking a set of arguments which describe the parameters defining the task.

The current version is in the beta stage. All comments and suggestions are welcome: please email [rkochanov@cfa.harvard.edu](mailto:rkochanov@cfa.harvard.edu).

## Features

Features

Some of the prominent current features of HAPI are:

    1) Downloading line-by-line data from the HITRANonline site to a local machine;
    2) Filtering and processing the data in SQL-like fashion;
    3) Conventional Python structures (lists, tuples, and dictionaries) for representing spectroscopic data;
    4) Compatibility with a large set of third-party Python libraries to work with the data;
    5) Python implementation of the HT profile [2,3,4,5] which can be used in spectral simulations. This line shape can also be reduced to a number of conventional line profiles such as Gaussian (Doppler), Lorentzian, Voigt, Rautian, Speed-dependent Voigt and speed-dependent Rautian;
    6) Python implementation of the total internal partition sums algorithm, TIPS-2017[6] which is used in the calculation of the temperature dependence of HITRAN[7] transition intensities. The older software TIPS-2011[8] is also available;
    7) High-resolution spectral simulation accounting for pressure, temperature and optical path length. The following spectral functions can be calculated:
        a) absorption coefficient
        b) absorption spectrum
        c) transmittance spectrum
        d) radiance spectrum
    8) Spectral calculation using a number of instrumental functions to simulate experimental spectra;
    9) Possibility to extend the user's functionality by adding custom line shapes, partition sums and apparatus functions.

## Citation

It is free to use HAPI. If you use HAPI in your research or software development, please cite it using the following reference:

R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016) [Link to article](http://dx.doi.org/10.1016/j.jqsrt.2016.03.005).

To make a reference to particular version of HAPI, use corresponding DOI from the [Zenodo](https://zenodo.org/collection/user-hapi) community in addition to the reference given above. 

## References

[1] R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016) [Link to article](http://dx.doi.org/10.1016/j.jqsrt.2016.03.005).

[2] N. H. Ngo, D. Lisak, H. Tran, J.-M. Hartmann, An isolated line-shape model to go beyond the Voigt profile in spectroscopic databases and radiative transfer codes, J. Quant. Spectrosc. Radiat. Transfer 129, 89-100 (2013). [Link to article](http://www.sciencedirect.com/science/article/pii/S0022407313002598)

[3] H. Tran, N. H. Ngo, J.-M. Hartmann, Efficient computation of some speed-dependent isolated line profiles, J. Quant. Spectrosc. Radiat. Transfer 129, 199-203 (2013) [Link to article](http://www.sciencedirect.com/science/article/pii/S0022407313002598)

[4] H. Tran, N. H. Ngo, J.-M. Hartmann, Erratum to "Efficient computation of some speed-dependent isolated line profiles", J. Quant. Spectrosc. Radiat. Transfer 134, 104 (2014) [Link to article](http://www.sciencedirect.com/science/article/pii/S0022407313004445)

[5] J. Tennyson, P. F. Bernath, A. Campargue et al., Recommended isolated-line profile for representing high-resolution spectroscopic transitions (IUPAC Technical Report), Pure Appl. Chem. 86, 1931-1943 (2014) [Link to article](http://www.degruyter.com/view/j/pac.2014.86.issue-12/pac-2014-0208/pac-2014-0208.xml)

[6] R. R. Gamache, C. Roller, E. Lopes, I. E. Gordon, L. S. Rothman, et al., Total internal partition sums for 166 isotopologues of 51 molecules important in planetary atmospheres: Application to HITRAN2016 and beyond, J. Quant. Spectrosc. Radiat. Transfer 203, 70-87 (2017). [Link to article](https://www.sciencedirect.com/science/article/pii/S0022407317301516)

[7] I. E. Gordon, L. S. Rothman, C. Hill, R. V. Kochanov, Y. Tan, et al., The HITRAN2016 molecular spectroscopic database, J. Quant. Spectrosc. Radiat. Transfer 203, 3-69 (2017). [link to article](https://www.sciencedirect.com/science/article/pii/S0022407317301073)

[8] A. L. Laraia, R. R. Gamache, J. Lamouroux, I. E. Gordon, L. S. Rothman, Total internal partition sums to support planetary remote sensing, Icarus 215, 391-400 (2011). [Link to article](http://www.sciencedirect.com/science/article/pii/S0019103511002132)

