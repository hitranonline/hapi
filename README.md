# HITRAN Application Programming Interface (HAPI)
===============================================

Current version: 1.1.0.8.7

## Version history

FIXED GRID BUG (ver. 1.1.0.1)
FIXED OUTPUT FORMAT FOR CROSS-SECTIONS (ver. 1.1.0.1)
ADDED CPF BY SCHREIER (JQSRT_112_2011) (ver. 1.1.0.2)
OPTIMIZED EXPRESSION EVALUATIONS FOR SELECT (ver. 1.1.0.3)
ADDED SUPPORT FOR MIXTURES (ver. 1.1.0.4)
ADDED SUPPORT FOR USER-DEFINED ENV DEPENDENCES (ver. 1.1.0.5)
ADDED PROFILE SELECTION (ALPHA) (ver. 1.1.0.6)
ADDED METADATA FOR HTP, FIXED NORMALIZATION IN CONVOLVESPECTRUMSAME (ver. 1.1.0.7)
FIXED A "LONELY HEADER" BUG IN CACHE2STORAGE (ver. 1.1.0.7.1)
ADDED SUPPORT FOR PHOSGENE AND CYANOGEN (ver. 1.1.0.7.2)
OPTIMIZED STORAGE2CACHE (by Nils-Holger Löber) (ver. 1.1.0.7.3)
ADDED SKIPABLE PARAMETERS IN HEADERS (ver. 1.1.0.7.4)
ADDED SUPPORT FOR FORTRAN D-NOTATION (ver. 1.1.0.7.5)
ADDED SUPPORT FOR WEIRD-FORMATTED INTENSITY VALUES E.G. "2.700-164" (ver. 1.1.0.7.6)
ADDED TIPS-2017 (ver. 1.1.0.8)
ADDED SUPPORT FOR CUSTOM EXTENSIONS OF THE DATA FILES (ver. 1.1.0.8.1)
FIXED LINK TO (2,0) ISOTOPOLOGUE IN TIPS-2017 (ver. 1.1.0.8.2)
ADDED SAVEHEADER FUNCTION (ver. 1.1.0.8.3)
ADDED METADATA FOR SF6 (ver. 1.1.0.8.4)
ADDED D2O ISOTOPOLOGUE OF WATER TO DESCRIPTION (ver. 1.1.0.8.5)
FIXED LINE ENDINGS IN STORAGE2CACHE AND QUERYHITRAN (ver. 1.1.0.8.6)
ADDED SUPPORT FOR NON-INTEGER LOCAL ISO IDS (ver. 1.1.0.8.7)

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

