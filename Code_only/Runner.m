clear;

% ===== CSD Calculation Software ========
% Implementation based on the work of Reboul et al 2008 (https://doi.org/10.1007/s10035-008-0111-5) ...
% and Reboul et al., 2010 (https://doi.org/10.1016/j.compgeo.2009.09.002)

% Please cite this implementation as Shire et al 2013
% (https://doi.org/10.1680/geot.11.T.020) 

% Developed for Linux OS (Ubunutu 16.04). CGAL calculation will need
% adapting for windows OS. 

%% Input parameters

overlap = 0.5;  %% Overlap parameter as a fraction (Shire and O'Sullivan 2016  (doi:10.1680/jgeot.15.P.215))

filename = 'example_code';

%%

if exist('weighted_delaunay', 'file') == 0
    fprintf('\n Error: Please compile the CGAL executable "weighted_delaunay" and place in this folder" (see README) \n')
    return;
    
elseif exist('Input_Data.txt', 'file') == 0
    fprintf('\n Error: Please create "Input_Data.txt" file input with particle coordinates and radius (see README) \n')
    return;
    
else

    if isunix == 0;
       fprintf('\n Warning: CGAL code developed for linux systems! \n');
    end
    
unix ./weighted_delaunay; % Carry out CGAL regular triangulation weighted by particle radius 

TRI02_Read_CGAL % Read CGAL output

[dt, merge] = TRI03_cellmerge(overlap, filename); 

TRI04_constrictionsize_weighted_delaunay(filename);

    fprintf('\n Complete. Please see README for details of output \n')

end

