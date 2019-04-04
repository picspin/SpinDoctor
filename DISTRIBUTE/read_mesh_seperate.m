function [mymesh] = read_mesh_seperate(fname, sigma0)

% create FE mesh on canonical configuration; bend and twist the FE mesh nodes by analytical transformation
% 
% Input:
%     1. fname
%     
% Output:
%     1. mymesh is a structure with 10 elements:
%         Nnode
%         Nele
%         Nface
%         Pts_cmpt_reorder
%         Ele_cmpt_reorder
%         Pts_ind
%         Pts_boundary_reorder
%         Fac_boundary_reorder
%         Nboundary
%         Ncmpt   


% [Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder, Nboundary,Ncmpt]

disp(['Reading from neuron FE mesh from ', fname]);

mymodel1 = createpde();
mymodel2 = createpde();
elements1 = load([fname, '_dendrites_elements.txt']);
nodes1 = load([fname, '_dendrites_nodes.txt']);
elements2 = load([fname, '_soma_elements.txt']);
nodes2 = load([fname, '_soma_nodes.txt']);

geometryFromMesh(mymodel1,nodes1',elements1');
geometryFromMesh(mymodel2,nodes2',elements2');
%pdeplot3D(mymodel);

%%%%%%%%%%%%%
specifyCoefficients(mymodel1,'m',0,'d',1,'c',sigma0,'a', 0,'f',0);
nface = mymodel1.Geometry.NumFaces; %number of face
%     normal dot (c grad u) + q*u = g
for iface = 1:nface
    applyBoundaryCondition(mymodel1,'neumann','face',iface,'g',1,'q',1,'Vectorized','on'); % 3-D geometry
end
%%%%%%%%%%%%%%
%%%%%%%%%%%%%
specifyCoefficients(mymodel2,'m',0,'d',1,'c',sigma0,'a', 0,'f',0);
nface = mymodel2.Geometry.NumFaces; %number of face
%     normal dot (c grad u) + q*u = g
for iface = 1:nface
    applyBoundaryCondition(mymodel2,'neumann','face',iface,'g',1,'q',1,'Vectorized','on'); % 3-D geometry
end
%%%%%%%%%%%%%%


[elems2faces1, faces2nodes1] = get_faces(elements1);
faces2elems1=entryInWhichRows(elems2faces1); 
bindex1=find(faces2elems1(:, 2)==0);           %index of boundary facets - belong to one element only
boundary1=faces2nodes1(bindex1, :);             %nodes of boundary facets 
facets1=sort(boundary1, 2);

[elems2faces2, faces2nodes2] = get_faces(elements2);
faces2elems2=entryInWhichRows(elems2faces2); 
bindex2=find(faces2elems2(:, 2)==0);           %index of boundary facets - belong to one element only
boundary2=faces2nodes2(bindex2, :);             %nodes of boundary facets 
facets2=sort(boundary2, 2);


mymesh.Ncmpt = 2;
mymesh.Nboundary = 2;

mymesh.Nnode = [size(nodes1,1); size(nodes2,1)];
mymesh.Nele = [size(elements1,1); size(elements2,1)];
mymesh.Nface = [mymodel1.Geometry.NumFaces; mymodel2.Geometry.NumFaces];

mymesh.Pts_cmpt_reorder{1} = mymodel1.Mesh.Nodes;
mymesh.Ele_cmpt_reorder{1} = mymodel1.Mesh.Elements;
mymesh.Pts_ind{1} = [1:size(mymodel1.Mesh.Nodes,2)]';

mymesh.Pts_cmpt_reorder{2} = mymodel2.Mesh.Nodes;
mymesh.Ele_cmpt_reorder{2} = mymodel2.Mesh.Elements;
mymesh.Pts_ind{2} = [1:size(mymodel2.Mesh.Nodes, 2)]';

mymesh.Pts_boundary_reorder{1}{1} = [1:size(mymodel1.Mesh.Nodes,2)]';
mymesh.Fac_boundary_reorder{1}{1} = facets1';
mymesh.Pts_boundary_reorder{1}{2} = [];
mymesh.Fac_boundary_reorder{1}{2} = [];

mymesh.Pts_boundary_reorder{2}{2} = [1:size(mymodel2.Mesh.Nodes,2)]';
mymesh.Fac_boundary_reorder{2}{2} = facets2';

end
