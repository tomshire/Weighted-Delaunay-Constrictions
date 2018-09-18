Constriction size distribution calculation.  
Implementation of Reboul et al 2008 (https://doi.org/10.1007/s10035-008-0111-5)  and Reboul et al., 2010 (https://doi.org/10.1016/j.compgeo.2009.09.002)
Please cite this Implementation as Shire et al 2013 (https://doi.org/10.1680/geot.11.T.020). 
The implementation has been developed for a linux OS and uses the Matlab command ‘unix’ . To run in windows this stage will need to be adapted.  
For details on how to compile the CGAL code see readme in Weighted_Delaunay_CGAL folder.


#######  Input  ########
The input to the code is a text file named ‘Input_Data.txt’  (case sensitive). 
The file is a space delimited text file with no headers columns:
x y z radius

########  Running the code  ########
The code is run from Runner.m
The overlap parameter (see Stage 2) and a filename identifier must be specified. 
There are three stages to the implementation:

1.	Weighted Delaunay triangulation of particle centres using CGAL “Regular Triangulation”  https://doc.cgal.org/latest/Triangulation_3/index.html#Triangulation3secclassRegulartriangulation 
As noted above, this stage has been developed for a Debian-based linux OS.  For details on how to compile the CGAL code see readme in Weighted_Delaunay_CGAL folder.
TRI02_Read_CGAL reads the output from this file. 

2.	Merge Delaunay tetrahedra to create larger voids
TRI03_cellmerge.m
Here the overlap parameter must be input be specified as a fraction of the smaller particle (see Shire et al 2013 for details). 
Note that a proportion of Delaunay cells are deleted to avoid edge effects (currently any within 5% of the edge of the sample – this can be altered as necessary..

3.	Constriction calculation
TRI04_constrictionsize_weighted_delaunay.m
This takes the output from the void merging calculation and calculates constrictions on the face of each void, with appropriate checks as set out by Reboul et al 2010. 

###### Output ########
All output uses the matlab .mat file type

From Stage 2 (void merging):

voids_{filename}.mat:
‘dt’: contains the Delaunay tetrahedra verticies and vertex coordinates. Can be used with functions related to the matlab triangulation class
‘merge’: Each row corresponds to a tetrahedron whose Vertices are listed in ‘dt’.  Note a zero value here means that the tetrahedron was deleted to avoid edge effects
‘posrad’: The position and radius of each DEM particle

From Stage 3 (constriction calculation):

constrictions_{filename}.mat:
Rows 1:3 – coordinates of circle centre
Row 4 – constriction diameter
Rows 5:7 – normal vector to constriction face

particles_{filename}.mat:
Row 1 - Particle ID
Rows 2:4 – Particle centroid
Row 5 – Particle radius
