clear;

%% Read .dat file
data = read_mumat('Solid sphere\Spheres\sphere_mu.dat');
ntet = data.ntet;
maxNB = 10;

%% Calculate tetrahedron centers
tet_cen = zeros(3,ntet);
for i = 1:ntet
    tet_cen(:,i) = sum(data.coords(:,data.tet(:,i)),2)/4;
end

%% Find nearest neighbors
neighbours = zeros(maxNB,ntet);
for i = 1:ntet
    dx = tet_cen(:,:)-tet_cen(:,i);
    dist = vecnorm(dx);
    dist(i) = nan;
    for j = 1:maxNB

        [M, I] = min(dist,[],'omitnan');
        neighbours(j,i) = I;
        dist(I) = nan;
    end
end

%% Read
fileID = fopen('test.txt', 'r');
formatSpec = '%i';
sizeA = [11 Inf];
simdata = fscanf(fileID,formatSpec,sizeA)';

lerr = false;
for i = 1:length(simdata)
    ind = simdata(i,1);
    dif = sum(simdata(i,2:end)-neighbours(:,ind)');
    if dif ~= 0
        disp(ind)
        lerr = true;
    end
end
if lerr
    disp('Difference detected')
end
%B = sortrows(simdata,1,'ascend');
%B = B(:,2:end);

return
%% Visualize

tet_cen_x = tet_cen(1,:);
tet_cen_y = tet_cen(2,:);
tet_cen_z = tet_cen(3,:);

tet_test = 1;
test_cen_x = tet_cen_x(tet_test); 
test_cen_y = tet_cen_y(tet_test);
test_cen_z = tet_cen_z(tet_test);

tet_verts = data.coords(:,data.tet(:,tet_test))';
tet_verts_x = tet_verts(:,1);
tet_verts_y = tet_verts(:,2);
tet_verts_z = tet_verts(:,3);

test_neighbours = neighbours(:,tet_test);
test_n_x = tet_cen_x(test_neighbours);
test_n_y = tet_cen_y(test_neighbours);
test_n_z = tet_cen_z(test_neighbours);



%%
figure;
scatter3(tet_cen_x, tet_cen_y, tet_cen_z,'b')
hold on

scatter3(test_cen_x, test_cen_y, test_cen_z, 'r','MarkerFaceColor','r')
scatter3(tet_verts_x, tet_verts_y, tet_verts_z, 'r.')
plot3(tet_verts_x, tet_verts_y, tet_verts_z,'r-')
plot3([tet_verts_x(1) tet_verts_x(3)], ...
    [tet_verts_y(1) tet_verts_y(3)], ...
    [tet_verts_z(1) tet_verts_z(3)],'r-')
plot3([tet_verts_x(1) tet_verts_x(4)], ...
    [tet_verts_y(1) tet_verts_y(4)], ...
    [tet_verts_z(1) tet_verts_z(4)],'r-')
plot3([tet_verts_x(2) tet_verts_x(4)], ...
    [tet_verts_y(2) tet_verts_y(4)], ...
    [tet_verts_z(2) tet_verts_z(4)],'r-')

scatter3(test_n_x, test_n_y, test_n_z, 'g','MarkerFaceColor','g')

d = 100;
xlim([test_cen_x-d test_cen_x+d])
ylim([test_cen_y-d test_cen_y+d])
zlim([test_cen_z-d test_cen_z+d])

% r = norm(tet_cen(:,tet_test)-tet_cen(:,test_neighbours(end-1)));
% [x,y,z] = sphere;
% x=x*r; y = y*r; z = z*r;
% surf(x+test_x,y+test_y, z+test_z)