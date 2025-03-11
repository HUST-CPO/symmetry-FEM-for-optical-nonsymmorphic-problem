function[SrcNodeIndex,DstNodeIndex,edgephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,edgephi,Mesh)

index1=[];
index2=[];
for i = 1:size(DstNodeIndex,1)
    rows = [1,2];
    dif=Mesh.Nodes(rows,Mesh.Edge(SrcNodeIndex(i,1),:))-[0,0;0,0];
    for j = 1:size(DstNodeIndex,1)
        if(DstNodeIndex(j,1)==SrcNodeIndex(i,1))
            index1=[index1,i];
            index2=[index2,j];
        end
    end
end
for i = 1:size(index1,1)
    SrcNodeIndex(index1(i),1)=SrcNodeIndex(index2(i),1);
    edgephi(index1(i),1)=edgephi(index1(i),1)*edgephi(index2(i),1);
end

