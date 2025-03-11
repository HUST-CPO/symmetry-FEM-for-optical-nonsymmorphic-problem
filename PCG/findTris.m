function EdgeIndex=findTris(targetflag1,Mesh)
index=find(Mesh.Edgesflag==targetflag1(1,1));
for i=2:size(targetflag1,2)
    index=[index;find(Mesh.Edgesflag==targetflag1(1,i))];
end
index=sort(index);
Boundary1=Mesh.Edges(index,:);
[~ ,PBCEdge] = ismember(Boundary1,Mesh.Edges,'rows');
EdgeIndex=sort(PBCEdge);