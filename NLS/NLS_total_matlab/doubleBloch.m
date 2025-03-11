function[SrcNodeIndex,DstNodeIndex,edgephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,edgephi,Mesh)
num = 0;
B = [];
C=[];
D=[];
for i = 1:size(DstNodeIndex,1)
    rows = [1,2];
    dif=Mesh.Nodes(rows,Mesh.Edge(SrcNodeIndex(i,1),:))-[-2.25e-3,-2.25e-3;-2.25e-3,-2.25e-3];
    if(all(all(abs(dif)<0.000000001)))
        for j = 1:size(DstNodeIndex,1)
            if(DstNodeIndex(i,1)==SrcNodeIndex(j,1))
                SrcNodeIndex(j,1)=SrcNodeIndex(i,1);
                edgephi(j,:)=edgephi(j,:).*edgephi(i,:);
                B=[B;SrcNodeIndex(i,1),DstNodeIndex(i,1);SrcNodeIndex(j,1),DstNodeIndex(j,1)];
            end
        end
    end
    rows = [1,3];
    dif=Mesh.Nodes(rows,Mesh.Edge(SrcNodeIndex(i,1),:))-[-2.25e-3,-2.25e-3;-1e-3,-1e-3];
    if(all(all(abs(dif)<0.000000001)))
        for j = 1:size(DstNodeIndex,1)
            if(DstNodeIndex(i,1)==SrcNodeIndex(j,1))
                SrcNodeIndex(j,1)=SrcNodeIndex(i,1);
                edgephi(j,:)=edgephi(j,:).*edgephi(i,:);
                C=[C;SrcNodeIndex(i,1),DstNodeIndex(i,1);SrcNodeIndex(j,1),DstNodeIndex(j,1)];
            end
        end
    end
    rows = [2,3];
    dif=Mesh.Nodes(rows,Mesh.Edge(SrcNodeIndex(i,1),:))-[-2.25e-3,-2.25e-3;-1e-3,-1e-3];
    if(all(all(abs(dif)<0.000000001)))
        for j = 1:size(DstNodeIndex,1)
            if(DstNodeIndex(i,1)==SrcNodeIndex(j,1))
                SrcNodeIndex(j,1)=SrcNodeIndex(i,1);
                edgephi(j,:)=edgephi(j,:).*edgephi(i,:);
                num = num+1;
                D=[D;SrcNodeIndex(i,1),DstNodeIndex(i,1);SrcNodeIndex(j,1),DstNodeIndex(j,1)];
            end
        end
    end
end
A=[SrcNodeIndex,DstNodeIndex,edgephi];
A=unique(A,"rows");
SrcNodeIndex=A(:,1);
DstNodeIndex=A(:,2);
edgephi=A(:,3);
