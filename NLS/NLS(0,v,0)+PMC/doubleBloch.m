function[SrcNodeIndex,DstNodeIndex,edgephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,edgephi,Mesh)

index=[];
for i = 1:size(DstNodeIndex,1)
    rows = [1,2];
    dif=Mesh.Nodes(rows,Mesh.Edge(SrcNodeIndex(i,1),:))-[2.25e-3,2.25e-3;0,0];
    if(all(all(abs(dif)<0.000000001)))
        index=[index,i];
    end
    rows = [1,2];
    dif=Mesh.Nodes(rows,Mesh.Edge(SrcNodeIndex(i,1),:))-[-2.25e-3,-2.25e-3;2.25e-3,2.25e-3];
    if(all(all(abs(dif)<0.000000001)))
        index=[index,i];
    end
end
SrcNodeIndex(index)=[];
DstNodeIndex(index)=[];
edgephi(index,:)=[];
