function [dt, merge] = TRI03_cellmerge(overlaptol, filename)
% <overlap tolerance>, <filename prefix
% for output>

%% Efficient Reboul merging process
% note that overlaptol is NOT % but in per one (50% -> 0.5)

tic
%% STEP 1: Delaunay tessellation

%Import triangulation from CGAL

tri_data = load('Del_points_and_cells');
% overlaptol = 0.5;
% filename = 'CU6_6000_PART';

posrad = tri_data.Points;
posrad = posrad;
Vertices = tri_data.Vertices;

clear tri_data

dt= triangulation(Vertices, posrad(:,1), posrad(:,2), posrad(:,3));
dtemp=dt.ConnectivityList;


%% STEP 2: Delete tetrahedra on the edges

% I'm going to delete the tetrahedra whose circumcentre is not in the 90%
% central box.

incentres=zeros(size(dt,1),3);
for a=1:size(incentres,1)
    partcoordaux=[posrad(dt(a,1),1:3); posrad(dt(a,2),1:3); posrad(dt(a,3),1:3); posrad(dt(a,4),1:3)];
    incentres(a,1)=mean(partcoordaux(:,1));
    incentres(a,2)=mean(partcoordaux(:,2));
    incentres(a,3)=mean(partcoordaux(:,3));
end


% max_part_diam = min(posrad(:,4))*2; %% Use minimum radius, as sample is so small! (normally use maximum!)

cell_length_x = max(posrad(:,1)) - min(posrad(:,1));
cell_length_y = max(posrad(:,2)) - min(posrad(:,2));
cell_length_z = max(posrad(:,3)) - min(posrad(:,3));



minacceptx = min(posrad(:,1))+  0.05*cell_length_x;
minaccepty = min(posrad(:,2))+  0.05*cell_length_y;
minacceptz = min(posrad(:,3))+  0.05*cell_length_z;

maxacceptx= max(posrad(:,1))-  0.05*cell_length_x;
maxaccepty= max(posrad(:,2))-  0.05*cell_length_y;
maxacceptz= max(posrad(:,3))-  0.05*cell_length_z;

dtdel=zeros(0,1);
for b=1:size(incentres,1)
    if incentres(b,1)<minacceptx || incentres(b,1)>maxacceptx || incentres(b,2)<minaccepty || incentres(b,2)>maxaccepty || incentres(b,3)<minacceptz || incentres(b,3)>maxacceptz
        % then it's out of edges, delete
        dtdel=[dtdel;b];
    end
end


%% STEP 3: For each tetrahedra calculate the inscribed sphere

disp('LVL2_Calculate inscribed spheres')
% have to produce the centre and radius
inscribedsph=zeros(size(dtemp,1),4);

for c=1:size(inscribedsph,1)
    if mod(c,2500)==0
        c;
        fprintf('\n %i of %i inscribed spheres calculated\n', c, size(inscribedsph,1));
    end
    if size(find(dtdel==c),1)==0
        % this thing here should solve the inscribed sphere to the tetrahedron
        % 4 unknowns 4 equations (one for each sphere)
        
        % NOTE: if you write %f instead of %d the solver inputs a float
        % number. USE THAT for non-vox data.
        
         syms r x y z;       
        
        var1 = posrad(dtemp(c,1),1:4);        
        [strg1]=((var1(1)-x)^2+(var1(2)-y)^2+(var1(3)-z)^2 == (r+var1(4))^2);
        
        var2 = posrad(dtemp(c,2),1:4);        
        [strg2]=((var2(1)-x)^2+(var2(2)-y)^2+(var2(3)-z)^2 == (r+var2(4))^2);        

        var3 = posrad(dtemp(c,3),1:4);        
        [strg3]=((var3(1)-x)^2+(var3(2)-y)^2+(var3(3)-z)^2 == (r+var3(4))^2); 
        
        
        var4 = posrad(dtemp(c,4),1:4);        
        [strg4]=((var4(1)-x)^2+(var4(2)-y)^2+(var4(3)-z)^2 == (r+var4(4))^2);
        

        
        [R,Cx,Cy,Cz]= solve (strg1,strg2,strg3,strg4);
        
        resultsofsph=[eval(Cx) eval(Cy) eval(Cz) eval(R)];
        
        %sort in descending radii
        resultsofsph=sortrows(resultsofsph,-4);
        inscribedsph(c,1:4)=resultsofsph(1,:);
       
    end
end

disp ('LVL2_Inscribed spheres done')
%% STEP 4: merge the tetrahedra seing the overlap

merge=zeros(size(dtemp,1),1);
overlaps_save=zeros(0,11);

disp('Merge process begun')

count=1;
for d=1:size(merge,1)
    if mod(d,2500)==0
        d;
        fprintf('\n %i of %i merge checks carried out\n', d, size(merge,1));
    end
    if size(find(dtdel==d),1)==0
        neigths1=neighbors(dt,d);
        
        neigths1(find(isnan(neigths1)==1)) = [];
        
        % level 2: nieghbours of the neighbours
        neights2=neigths1;
        for z=1:length(neights2)
            neights2=[neights2, neighbors(dt,neigths1(z))];
        end
        
        neights2(find(isnan(neights2)==1)) = [];
        
        neigths=unique(neights2);
        %delete d
        position=find(neigths==d);
        if position==1
            neigths=neigths(2:length(neigths));
        elseif position==length(neigths)
            neigths=neigths(1:length(neigths)-1);
        else
            neigths=[neigths(1:position-1), neigths(position+1:length(neigths))];
        end
        %done!
        for  e=1:length(neigths)
            if size(find(dtdel==neigths(e)),1)==0
                % if the neighbour should be considered
                overlap=(-1)*((sqrt((inscribedsph(d,1)-inscribedsph(neigths(e),1))^2+(inscribedsph(d,2)-inscribedsph(neigths(e),2))^2+(inscribedsph(d,3)-inscribedsph(neigths(e),3))^2)-inscribedsph(d,4)-inscribedsph(neigths(e),4))/(min(inscribedsph(d,4),inscribedsph(neigths(e),4))));
                overlaps_save=[overlaps_save; d, inscribedsph(d,:), neigths(e), inscribedsph(neigths(e),:), overlap];
                if overlap>=overlaptol
                    % then they belong to the same void
                    %check if any of them already labeled: can be none, 1 or
                    %both labelled
                    if merge(d)==0 && merge(neigths(e))==0
                        merge(d)=count;
                        merge(neigths(e))=count;
                        count=count+1;
                    elseif merge(d)==0
                        merge(d)=merge(neigths(e));
                    elseif merge(neigths(e))==0
                        merge(neigths(e))=merge(d);
                    else
                        % replace all the voids labelled with the 'e' label by
                        % the 'd' one
                        poss=find(merge==merge(neigths(e)));
                        for f=1:length(poss)
                            merge(poss(f))=merge(d);
                        end
                    end
                else
                    if merge(d)==0
                        merge(d)=count;
                        count=count+1;
                    end
                end
            end
        end
    end
end


% last: re-label the voids from 1 to n

hss=histc(merge,(1:max(merge)));
relbl=zeros(max(merge)-sum(~histc(merge,(1:max(merge)))));
count2=1;
for g=1:length(hss)
    if hss(g)~=0
        relbl(count2)=g;
        count2=count2+1;
    end
end


%first I create a vector with the positions
rlbl_oposite=zeros(max(max(relbl)));

%then I fill that vector as follows:
count3=1;
for h=1:length(rlbl_oposite)
    if h==relbl(count3)
        rlbl_oposite(h)=count3;
        count3=count3+1;
    end
end

for h=1:length(merge)
    if merge(h)~=0
        merge(h)=rlbl_oposite(merge(h));
    end
end

% end

%% STEP 5: save data

disp('Saving data')

save (['voids_' filename],'posrad','dt','merge');

toc

end