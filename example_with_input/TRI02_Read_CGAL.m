function TRI02_Read_CGAL

%% Read input from CGAL regular delaunay and save .mat file
%containing points and vertices

% Read input (created by CGAL)
A = dlmread('Output_Data.txt');

% Split off points
num_points = A(2,1); % row in text file which gives number of points

Points = A(3:num_points+2, :); % POints follow this line

% CHANGE WEIGHT FROM radius^2 to radius!!!

Points(:,4) = sqrt(Points(:,4));

num_vert = A(num_points + 3, 1); % row in text file which gives number of cells


Vertices = A((num_points + 4):(num_points + num_vert + 3),:); % IDs of vertex points

clear A;

save('Del_points_and_cells.mat', 'Points', 'Vertices');

end