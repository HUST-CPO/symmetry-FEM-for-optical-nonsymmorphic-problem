function mesh=GetEdge(mesh)

%边界索引
el2no=mesh.Elements';
n1=el2no([1 1 2],:);
n2=el2no([2 3 3],:);
el_ed2no_array=[n1(:) n2(:)];
[mesh.Edge,~,mesh.EdgeOfElements]=unique(el_ed2no_array,'rows');%EdgesOfElements网格边索引
mesh.NbrEdge=length(mesh.Edge);

end