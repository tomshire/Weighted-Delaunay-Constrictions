% constrictsizefinalversion.m
% 21/04/2011
% gordon gaudray

% this program find the distribution of the contrictions in a packing
%It loads a file which contains the coordinates and the radius of the
%particles of the studied sample.Then it applies a Delaunay triangulation
%on teh sample and find the faces of the tetrahedra which are created.For
%each face, the algorithm find the largest disc which can fit between the three particles of the face thanks to a
%optimization method. It checks then if some particles are overlapped by it
%and in this case it finds the constriction size for all the combinations
%of the three particles  and all the overlapped ones and checks again if it
%meets others particles.
%It keeps at the end the largest constriction which belongs to the void
%space.
% Finally the results are stored in an array and the programm plot the
% distribution size of the constrictions and of the particles.

% fname is a string of the .mat file output from TRI03

function TRI04_constrictionsize_weighted_delaunay(fname)



%========================= PARAMETERS ==============================
% choice of the sample 

filename = ['voids_' fname];

info = load(filename);  % The loaded text file is stored in Q . Then the first column is the ID of the particle,columns 2 to 4 are the coordinates and the fifth column is the radius of the particle.
dt = info.dt;
merge = info.merge;
% face = info.listLVL2;
[out_faces face]=csd_facelist_creation (dt,merge);

% clear dt merge

Q=info.posrad;  % %ID, coords, radius
Q = [(1:size(Q,1))' Q];

clear info

%%optimization

%%starting point   (from which point the optimization algorithm will start) 

start=1; %isobarycenter


%%function choice

%  fun=1;  %fminunc  (based on gradient method) - seems to be no advantage
%  to using this approach
 fun=2;  %fminsearch (Nelder�Mead method) We have chosen this algorithm for the optimization because it is the fastest.It does not need to compute the gradient of the function to optimizate.
 
% % options (adjusted)

% options = [];
% options = optimset('MaxFunEvals',500);
% options = optimset('MaxIter',);
% options = optimset('TolFun',1e-006);
% options = optimset('TolX',1e-006);
options = optimset('TolX',1e-006,'MaxFunEvals',10000, 'Display','off');

%%tolerance for the test if particles are overlapped y a constriction.It
%%allows to avoid considering tangential particles as overlapped particles
%%because of the bad accuracy of the test '=0'.
tolerance=0.01;



%====================== DEFINE FUNCTIONS FOR LATER USE ==========

function bar=barycenter(p1,m1,p2,m2,p3,m3)
bar=(m2/(m1+m2+m3))*(p2-p1)+(m3/(m1+m2+m3))*(p3-p1)+p1;
end


function [outside_faces list2]=csd_facelist_creation (dt,merge)

    roughlist=zeros(0,3);

for i=1:max(merge)
    rows=find(merge==i);
    dtemp=zeros(0,4);
    for j=1:length(rows)
        dtemp=[dtemp; dt(rows(j),:)];
    end
    
    %particles defining that void
    parts=unique(dtemp);
    
    subpts=zeros(0,3);
    for k=1:length(parts)
        subpts=[subpts; dt.Points(parts(k),:)];
    end
    
    dtemp2=DelaunayTri(subpts);
    FF=freeBoundary(dtemp2);
    for l=1:size(FF,1)
        for m=1:size(FF,2)
            FF(l,m)=parts(FF(l,m));
        end
    end
    
    roughlist=[roughlist; FF];
end

% LIST considering stuff on the edges also


list=unique(roughlist,'rows');
list = sort(list,2,'ascend');

[list2, I,J] = unique(list, 'rows');

nn = hist(J, length(I)); % J gives the IDs of the first occurance of a 
                        % given Delaunay face in the original fullface
                        % array. From histogram of this in
                        % (length(facecheck)) number
                        % bins (i.e. one bin per unique face), can get
                        % whether a face appears once or twice!!!

nn = nn'; % If n == 1, outside face, if n == 2, inside face!!
% Facecheck indicies which correspond to inside faces
nn2 = find(nn == 1);


outside_faces = nn2;

% inside_faces = sort(inside_faces,2,'ascend');
% list=unique(inside_faces, 'rows');


end




%====================== CREATION OF THE SAMPLE =====================================
%this part create an array with the coordinates and the radius of the particles from the
%loaded file and uses it to create an array with the tetrahedra of the
%delaunay triangulation. then, it plots the packing.Note : at the end of
%this part , it is possible to display or not the triangulation and the
%axis




% coordinates of the spheres
x = Q(:,2);
y = Q(:,3);
z = Q(:,4);
r = Q(:,5);




%%
%===================== OPTIMIZATION =====================================



%% initialise sphere grids for later efficiency

% create array of which sphere is in which grid

% number of grids in each dimension
nboxx=5;
nboxy=5;
nboxz=5;    

%Min and max coordinates
xmin=min(Q(:,2)); 
ymin=min(Q(:,3));
zmin=min(Q(:,4));
xmax=max(Q(:,2));
ymax=max(Q(:,3));
zmax=max(Q(:,4));

%bound of the sample (to delete constrictions on the sides)
xborder=xmax; % warning, the box in wich will fit the particles has to be [1,xborder]*[1,yborder]*[1,zborder]
yborder=ymax;
zborder=zmax;

no_boxes = nboxx*nboxy*nboxz;
spherebox=zeros(no_boxes,1);

xindex = zeros(size(Q,1),1);
yindex = zeros(size(Q,1),1);
zindex = zeros(size(Q,1),1);
boxindex = zeros(size(Q,1),1);

    xindex=ceil(((Q(:,2)-xmin)./(xmax-xmin))*nboxx);
    yindex=ceil(((Q(:,3)-ymin)./(ymax-ymin))*nboxy);   
    zindex=ceil(((Q(:,4)-zmin)./(zmax-zmin))*nboxz);     
    %which box a given particle is in
    boxindex(:)=(zindex(:)-1)*nboxx*nboxy+nboxx*(yindex(:)-1)+xindex(:);


    Q_10 = uint16((length(Q))/10); %divide Q by 10 to get max size of array to store IDs
    sphere_ids = zeros(Q_10,no_boxes);
    
for k = 1:no_boxes
    num_terms = 0;
    num_terms = find(boxindex(:)==k);
    
    spherebox(k,1)=length(num_terms); % no of spheres in each box
    sphere_ids(1:length(num_terms),k)=num_terms ;
end

%find maximum number of particles in any box
[row_id, col_id]  = find(sphere_ids>0);  
max_row = max(row_id);
%Delete rows of sphere_ids with only zero terms
sphere_ids((max_row+1):Q_10,:) = [];




%%
%initialisation bis
coords=zeros(length(face),3); %coordinates of the center of the constriction
radius=zeros(length(face),1); % radius of the constriction
orientation=zeros(length(face),3);%orientation of the constriction
eval=zeros(length(face),1);  % value of the function to optimizate


%it  is realized here a loop on the faces of the tetrahedra to
%find the bicounter constriction of this face. 
for i=1:length(face);
    
     %define the three particles of the triangle
    apex1=Q(face(i,1),2:4);  % coordinates of the particle
    apex2=Q(face(i,2),2:4);
    apex3=Q(face(i,3),2:4);
    r1=Q(face(i,1),5); %radius of the particle
    r2=Q(face(i,2),5);
    r3=Q(face(i,3),5);
        
     %create a basis of the triangle to ensure the result would be in the
     %plane
     
    
     
    
     %define the new orthonormal basis (origin is at apex 1)
     u=(apex2-apex1)/norm(apex2-apex1);  % u,v orthogonal unit vectors in plane
     v=(u-(dot(u,u)/dot(apex3-apex1,u))*(apex3-apex1))/norm(u-(dot(u,u)/dot(apex3-apex1,u))*(apex3-apex1));
     w=cross(u,v)/norm(cross(u,v));  % normal to plane
     
    
     PP=[u(1),v(1),w(1);u(2),v(2),w(2);u(3),v(3),w(3)]; % transfer matrix
     Ref=[apex1(1);apex1(2);apex1(3)]; % reference point
    
   
     %starting point
       % barycenter is a matlab function defined in another matlab file.
     if start==1
        x0=barycenter(apex1,1,apex2,1,apex3,1); 
     elseif start==2
        x0=barycenter(apex1,r1,apex2,r2,apex3,r3); 
     elseif start==3
        x0=barycenter(apex1,1/r1,apex2,1/r2,apex3,1/r3);
     elseif start==4
        x0=barycenter(apex1,r1^3,apex2,r2^3,apex3,r3^3); 
     elseif start==5
        x0=barycenter(apex1,1/r1^3,apex2,1/r2^3,apex3,1/r3^3); 
     end
     
     % starting points in the new basis
     lambda01=dot(x0-apex1,u)/dot(u,u);
     lambda02=dot(x0-apex1,v)/dot(v,v);
     lambda0=[lambda01,lambda02];
     
     
     % create the function to optimizate : this function is equal to zero
     % when the constriction is touching the three particles which are
     % defining the triangle.
     %this is the function p.198 (3) in the Reboul paper
     
     % variable change xs=apex1+lamda1*u+lambda2*v
        constrcenter2 = @(lambda) (        sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex1(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex1(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex1(3))^2)-r1 ...
                                         - sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex3(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex3(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex3(3))^2)+r3)^2 ...
                                   +      (sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex1(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex1(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex1(3))^2)-r1 ...
                                         - sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex2(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex2(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex2(3))^2)+r2)^2 ...
                                   +      (sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex2(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex2(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex2(3))^2)-r2 ...
                                         - sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex3(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex3(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex3(3))^2)+r3)^2;
    
     if fun==1
         [mymin,fval] = fminunc(constrcenter2,lambda0,options);
     elseif fun==2
         [mymin,fval] = fminsearch(constrcenter2,lambda0,options);
     end
     
     % mymin(1)=x-value in local coordinate system
     % mymin(2)=y-value in local coordinate system     
     %%%%%%%%%%%%%%%%
%              %coordinates and radius of the solution in the original basis   
           

         coords(i,:)=(Ref+PP*[mymin(1);mymin(2);0])';  % Ref is the ref piont in the new system ie apex1
         %the radius is equal to the minimal distance from the center of
         %the constriction to one of the three others particle.
        % the radius is defined as the minimal distance from the center of
        % the constriction to the center of one of the apex minus the
        % radius of this apex
         radius(i)=min([sqrt((apex1(1)+mymin(1)*u(1)+mymin(2)*v(1)-apex1(1))^2+(apex1(2)+mymin(1)*u(2)+mymin(2)*v(2)-apex1(2))^2+(apex1(3)+mymin(1)*u(3)+mymin(2)*v(3)-apex1(3))^2)-r1, ...
                        sqrt((apex1(1)+mymin(1)*u(1)+mymin(2)*v(1)-apex3(1))^2+(apex1(2)+mymin(1)*u(2)+mymin(2)*v(2)-apex3(2))^2+(apex1(3)+mymin(1)*u(3)+mymin(2)*v(3)-apex3(3))^2)-r3, ...
                        sqrt((apex1(1)+mymin(1)*u(1)+mymin(2)*v(1)-apex2(1))^2+(apex1(2)+mymin(1)*u(2)+mymin(2)*v(2)-apex2(2))^2+(apex1(3)+mymin(1)*u(3)+mymin(2)*v(3)-apex2(3))^2)-r2]);    
         orientation(i,:)=w;
         %%%%%%%%%%%%%%%%%%%%%%
         
   %==============FIND THE OVERLAPPED PARTICLES==========
  %this part is creating two arrays, one with all the particles of the
  %sample which meet the plane defined by the triangle called
  %possiblepoints, and an other one called overlappedparticle which
  %is an array included in possiblepoints but with only the particles which
  %are overlapped by the initial constriction
   
   
   
     % points which could encounter the constriction
     count=1; % count the number of time that the piece of code between line 309 and 311 is executed
     overlappedparticle=zeros(5,5);   %array of th particle which are overlapped by the initial constriction
  
     possiblepoints=zeros(length(face),4);  % represents the circles defined by the intersection between the particles which encounter the face of the tetrahedron and this face
    % columns 1 to 3 are the coordinates of the center of this circle and
    % the fourth column is its radius
    t=0;
    tic;
    




%% Just loop through those spheres in grid boxes adjacent to the current grid box

% find box location of centre of delaunay triangle, and look in this and all
% adjacent boxes


% Find centre of triangle
tri_coord = zeros(3,1);

tri_coord(1) = (apex1(1)+apex2(1)+apex3(1))/3;
tri_coord(2) = (apex1(2)+apex2(2)+apex3(2))/3;
tri_coord(3) = (apex1(3)+apex2(3)+apex3(3))/3;

% Find box indexes
    box_xindex=ceil(((tri_coord(1)-xmin)/(xmax-xmin))*nboxx);
    box_yindex=ceil(((tri_coord(2)-ymin)/(ymax-ymin))*nboxy);   
    box_zindex=ceil(((tri_coord(3)-zmin)/(zmax-zmin))*nboxz);     
    %which box a given particle is in
    box=0;
    box=(box_zindex-1)*nboxx*nboxy+nboxx*(box_yindex-1)+box_xindex;
  
    % loop through this and adjacent boxes
    for h=-1:1
        for m=-1:1
            for n=-1:1
                
                searchbox=(box_zindex+n-1)*nboxx*nboxy+nboxx*(box_yindex+m-1)+box_xindex+h;
                % test against spheres stored in searchbox
                
                if searchbox>0 && searchbox <= no_boxes % the boxes that exist
                    %search all particles in box
                    lastpoint = 0;
                    lastpoint = find(sphere_ids(:,searchbox)>0,1,'last');
                    %loop over all the particle IDs stored in the box
                    for p = 1:lastpoint
                        % j is the particle ID in question
                        j = sphere_ids(p,searchbox);
                        %
                        
                        
                        
                        
                        
                        % distance from the particles to the plane defined by the triangle
                        testedpoint=Q(j,2:4);   % tested particles
                        rr=Q(j,5);
                        
                        distance= abs(dot(w,testedpoint)-dot(apex1,w))/norm(w);
                        % belongs to the void space?
                        if distance<rr
                            % this is a matrix computation which allows to determine the
                            % center of the projected circle (cf Reboul thesis p183 (D.12) )
                            M1=[w(1),w(2),w(3);u(1),u(2),u(3);v(1),v(2),v(3)];
                            M2=[dot(apex1,w);dot(testedpoint,u);dot(testedpoint,v)];
                            invertM1=[cross(M1(:,2),M1(:,3))';cross(M1(:,3),M1(:,1))';cross(M1(:,1),M1(:,2))']/det(M1);
                            % coordinates of the center of the projected circle
                            possiblecenter=(M1\M2)';
                            possiblepoints(j,:)=[possiblecenter,sqrt(rr^2-distance^2)];
                            
                            %test if the particle is overlapped by the constriction
                            if norm(possiblecenter-coords(i,:))+tolerance<radius(i)+sqrt(rr^2-distance^2)
                                overlappedparticle(count,:)=[3+count,possiblecenter,sqrt(rr^2-distance^2)]; % add the particle in the array
                                count=count+1;
                            end
                            
                            
                        end
                        
                        
                    end
                end
            end
        end
        
        
        
        
    end
    t=toc;
     
%%

     %===========================IMPROVED REBOUL ALGORITHM=========================
     %In this part a combination of three elements among the three apex of
     %the triangle and all the overlapped particle is done to create new triangles and the code
     %computes the new constrictions defined by these triangles and then it
     %compares them and keeps the bigger one as the real constriction if it
     %exists.
     
     

  
   if overlappedparticle(1,5)~= 0  % test if it exists overlappedparticles
         overlappedparticle(count:5,:)=[]; %delete the unused rows of overlapped particle
        index=zeros(1,count-1); % gives an index to the overlapped particle in order to do the combination
	newpoints=zeros(2+count,5);  % create an array with the three initilas particles of the triangle and the overlapped particles
        possiblepoints2=unique(possiblepoints,'rows');  % sort the array and delete all the rows which are equal to zero
        possiblepoints2(1,:)=[];
         for mm=1:(count-1);
             newpoints(3+mm,:)= overlappedparticle(mm,:);  %it fills the array newpoints with the overlapped particles
             index(1,mm) = overlappedparticle(mm,1);  %link the particle to a number
            
        end
            
         
  % the three initial points
         
            
	    newpoints(1,:)=Q(face(i,1),:); % fill the array newpoints with the initial points
	    newpoints(2,:)=Q(face(i,2),:);
        newpoints(3,:)=Q(face(i,3),:);

            Combinations=nchoosek([1,2,3,index],3); % realize the combinations of the particles in newpoints
            
          
        
% recompute the constriction for the new triangle (same thing that before)
    
            for ff=1:size(Combinations,1);
                apex1=newpoints(Combinations(ff,1),2:4);
                apex2=newpoints(Combinations(ff,2),2:4);
                apex3=newpoints(Combinations(ff,3),2:4);
         
               
            r1=newpoints(Combinations(ff,1),5);
            r2=newpoints(Combinations(ff,2),5);
            r3=newpoints(Combinations(ff,3),5);
            
           
           
            %new orthonormal basis 
     u=(apex2-apex1)/norm(apex2-apex1);
     v=(u-(dot(u,u)/dot(apex3-apex1,u))*(apex3-apex1))/norm(u-(dot(u,u)/dot(apex3-apex1,u))*(apex3-apex1));
    
         
     PP=[u(1),v(1),w(1);u(2),v(2),w(2);u(3),v(3),w(3)]; % transfer matrix
     Ref=[apex1(1);apex1(2);apex1(3)]; 
            
            newcoords=zeros(size(Combinations,1),3);
            newradius=zeros(size(Combinations,1),1);
            newmymin=zeros(size(Combinations,1),2);
            
            
                   
             if start==1
                x0=barycenter(apex1,1,apex2,1,apex3,1); 
            elseif start==2
                x0=barycenter(apex1,r1,apex2,r2,apex3,r3); 
            elseif start==3
                x0=barycenter(apex1,1/r1,apex2,1/r2,apex3,1/r3);
            elseif start==4
                x0=barycenter(apex1,r1^3,apex2,r2^3,apex3,r3^3); 
            elseif start==5
                x0=barycenter(apex1,1/r1^3,apex2,1/r2^3,apex3,1/r3^3); 
             end
     
            % starting points in the new basis
            lambda01=dot(x0-apex1,u)/dot(u,u);
            lambda02=dot(x0-apex1,v)/dot(v,v);
            lambda0=[lambda01,lambda02];
     
     
            % create the function to optimizate
            % variable change xs=apex1+lamda1*u+lambda2*v
             constrcenter2 = @(lambda) ( sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex1(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex1(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex1(3))^2)-r1 ...
                                         - sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex3(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex3(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex3(3))^2)+r3)^2 ...
                                   +      (sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex1(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex1(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex1(3))^2)-r1 ...
                                         - sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex2(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex2(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex2(3))^2)+r2)^2 ...
                                   +      (sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex2(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex2(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex2(3))^2)-r2 ...
                                         - sqrt((apex1(1)+lambda(1)*u(1)+lambda(2)*v(1)-apex3(1))^2+(apex1(2)+lambda(1)*u(2)+lambda(2)*v(2)-apex3(2))^2+(apex1(3)+lambda(1)*u(3)+lambda(2)*v(3)-apex3(3))^2)+r3)^2;
    
            if fun==1
                [dmymin,ffval] = fminunc(constrcenter2,lambda0,options);
            elseif fun==2
                [dmymin,ffval] = fminsearch(constrcenter2,lambda0,options);
            end
            
            %coordinates and radius of the solution in the original basis
         newcoords(ff,:)=(Ref+PP*[dmymin(1);dmymin(2);0])';
         newmymin(ff,:)=dmymin; % result in the new basis
         newradius(ff)=min([sqrt((apex1(1)+mymin(1)*u(1)+mymin(2)*v(1)-apex1(1))^2+(apex1(2)+mymin(1)*u(2)+mymin(2)*v(2)-apex1(2))^2+(apex1(3)+mymin(1)*u(3)+mymin(2)*v(3)-apex1(3))^2)-r1, ...
                        sqrt((apex1(1)+mymin(1)*u(1)+mymin(2)*v(1)-apex3(1))^2+(apex1(2)+mymin(1)*u(2)+mymin(2)*v(2)-apex3(2))^2+(apex1(3)+mymin(1)*u(3)+mymin(2)*v(3)-apex3(3))^2)-r3, ...
                        sqrt((apex1(1)+mymin(1)*u(1)+mymin(2)*v(1)-apex2(1))^2+(apex1(2)+mymin(1)*u(2)+mymin(2)*v(2)-apex2(2))^2+(apex1(3)+mymin(1)*u(3)+mymin(2)*v(3)-apex2(3))^2)-r2]);  
           
        
         
         
         % test here if the new constriction is overlapping one of the
         % initial point or one of the points which are in possiblepoints,
         % it means the particles which met the plane defined by the three
         % initial particles.
         for hh=1:size(possiblepoints2,1);
         
            %test if the new constriction is overlapping one of the particles which encounter the face of the tetrahedron defined by the three initial particle
             if  norm(possiblepoints2(hh,1:3)-newcoords(ff,:))+tolerance<newradius(ff)+possiblepoints2(hh,4)
                 newradius(ff)=0;
             end
         
         end
        
         
             end
             
          [emptydisc,I]=max(newradius);  % keep the bigger constriction among all the constrictions found by the combinaison
          
          coords(i,:)=newcoords(I,:);
          mymin=newmymin(I,:);
          radius(i)=newradius(I);
    end
         
        
         
     % delete solutions which are in the edges
     if coords(i,1)<xborder && coords(i,2)<yborder && coords(i,3)<zborder && xmin<coords(i,1) && ymin<coords(i,2) && zmin<coords(i,3) %&& radius(i)<tolerance2 
         
     
     else
    
         radius(i)=0;
         coords(i,:)=[0,0,0];
         orientation(i,:)=[0,0,0];
         
     
     end
   
   
     
end

    
CM_CONSTRIC=[coords,radius,orientation];
CM_CONSTRIC([out_faces],:)=[];

constrictions=unique(CM_CONSTRIC,'rows'); 
constrictions(1,:)=[];

particledata = Q;
%==============================DISTRIBUTION OF THE CONSTRICTIONS==========

% plot the distribution of the constrictions' radius and the
% particles'radius
% 
% figure
% hold on
% [alln,xout] = hist(SOL2(:,4),size(SOL2,1));
% 
% alldata=cumsum(alln);
% 
% plot(xout,100*(alldata/alldata(length(alldata))),'k');
% 
% clear alln xout;
% 
% [alln,xout] = hist(Q(:,5),size(Q,1));
% 
% alldata=cumsum(alln);
% 
% plot(xout,100*(alldata/alldata(length(alldata))),'g');
% 
% 
% xlabel('Constriction Size; Particle Radius');
% 
% xlim([0 max(r)]);
% 
% ylabel('% Smaller');
% 
% ylim([0 100]);
% 
% legend('Constriction size','Particle radius');


%================SAVE THE RESULTS==========



savefile1 = ['constrictions_' fname '.mat'];
savefile2 = ['particles_' fname '.mat'];
% savefile3 = ['histogram-' fname 'new.mat'];

save(savefile1,'constrictions');
save(savefile2,'particledata');
% save(savefile3,'oppositeconstr');
void = 1;

end
  
