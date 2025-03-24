function[SrcNodeIndex,DstNodeIndex,edgephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,edgephi,Mesh)
num = 0;
C=[];
for i = 1:size(DstNodeIndex,1)
    dif=Mesh.Nodes(Mesh.Edge(SrcNodeIndex(i,1),:),:);
    for j = 1:size(DstNodeIndex,1)
        if(all(all(abs(dif(:,1))<0.000000001)))
            if(DstNodeIndex(i,1)==SrcNodeIndex(j,1))
                SrcNodeIndex(j,1)=SrcNodeIndex(i,1);
                edgephi(j,:)=edgephi(j,:).*edgephi(i,:);
                num = num+1;
                C=[C;SrcNodeIndex(i,1),DstNodeIndex(i,1);SrcNodeIndex(j,1),DstNodeIndex(j,1)];
            end
        end
        if(all(all(abs(dif(:,2)+[0.008;0.008])<0.000000001)))
            if(DstNodeIndex(i,1)==DstNodeIndex(j,1)&&SrcNodeIndex(i,1)~=SrcNodeIndex(j,1))
                DstNodeIndex(j,1)=SrcNodeIndex(j,1);
                SrcNodeIndex(j,1)=SrcNodeIndex(i,1);
                edgephi(j,:)=(1./edgephi(j,:)).*edgephi(i,:);
                num = num+1;
                C=[C;SrcNodeIndex(i,1),DstNodeIndex(i,1);SrcNodeIndex(j,1),DstNodeIndex(j,1)];
            end
        end
    end
end
A=[SrcNodeIndex,DstNodeIndex];
[A,index,~]=unique(A,"rows");
SrcNodeIndex=A(:,1);
DstNodeIndex=A(:,2);
edgephi=edgephi(index,:);
