# Weighted-Delaunay-Constrictions
Calculates constriction size distribution for DEM simulations using the weighted Delaunay method

Implementation of Reboul et al 2008 (https://doi.org/10.1007/s10035-008-0111-5)  and Reboul et al., 2010 (https://doi.org/10.1016/j.compgeo.2009.09.002)
Please cite this Implementation as Shire et al 2013 (https://doi.org/10.1680/geot.11.T.020).
 
The implementation has been developed for a linux OS and uses the Matlab command ‘unix’ . To run in windows this stage will need to be adapted.  
For details on how to compile the CGAL code see readme in Weighted_Delaunay_CGAL folder.

There are four folders folders, three containing the same matlab code:

- Code_only (self explanatory)
- Validation_simple (validation using a simple cubic packing, for which analytical solutions are known - see Shire PhD thesis 2014 page 139 for more details)
- example_input (contains input file for the Cu = 1.5 sample from Shire & O'Sullivan 2016)

There is also a reference papers folder containng Shire et al 2013 and Shire & O'Sullivan 2016. 
