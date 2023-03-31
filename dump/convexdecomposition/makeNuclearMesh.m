clear all;

% A triangulation of m faces and n vertex is represented through:
% * a set of faces which is a (3,m) matrix where face(:,i) are the vertices indexes of the ith face.
% * a set of vertex which is a (d,n) matrix. 

addpath(genpath('toolbox_graph'));
setenv('PATH', ['triangle' ':' getenv('PATH') ]);

I = imread('rice.png');
mask = bwareaopen(imbinarize(I),10);
mask = imclearborder(mask);

[allBoundaryVertices, L] = bwboundaries(mask,'noholes');
allVertices = regionprops(L, 'PixelList');

i = 1;

% boundary vertices
Bverts = allBoundaryVertices{i};
bvidx = 1:size(Bverts,1);
Bedges = [bvidx', circshift(bvidx',1)];

% add inside vertices
% first removing with setdiff and then adding with cat is to ensure
% order, so Bedges indices are correct
% it also allows subsampling the inside without removing the outside
verts = fliplr(allVertices(i).PixelList);
verts = verts(1:8:end,:); % this may not be the smartest way to subsample
Iverts = setdiff(verts, Bverts, 'rows');
verts = cat(1, Bverts, Iverts);

holes = [];

% visualize:
imshow(L==i,[0 2])
hold on
scatter(Iverts(:,2),Iverts(:,1));
line(Bverts(:,2),Bverts(:,1),'Color','g')
hold off

%allvertices = cat(1, allBoundaryVertices{i, insideVertices);

% for this situation we may just be able to construct a triangulation
% without this costly general method that triangulates any point cloud
% which would be faster, but first we need proof of principle and worry
% about investing in speed later
tic
[V, F] = triangulate(verts, Bedges, holes);
toc

plot_mesh(V, F);

V3 = [V' 0*V(1,:)'];
write_ply(V3, F, 'nucleus.ply');