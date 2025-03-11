function[SrcNodeIndex,DstNodeIndex,edgephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,edgephi,Mesh)
for i = 1:size(DstNodeIndex,1)
    for j = 1:size(DstNodeIndex,1)
        if(SrcNodeIndex(i,1)==SrcNodeIndex(j,1))
            tempSrcNodeIndex=SrcNodeIndex(i,1);
            tempedgephi=edgephi(i,:);
        end
    end
end
for i = 1:size(DstNodeIndex,1)
    for j = 1:size(DstNodeIndex,1)
        if(DstNodeIndex(i,1)==DstNodeIndex(j,1)&&i~=j)
            ii=i;
            jj=j;
        end
    end
end
SrcNodeIndex(jj,1)=tempSrcNodeIndex;
edgephi(jj,:)=edgephi(jj,:).*tempedgephi;
SrcNodeIndex(ii,1)=tempSrcNodeIndex;
edgephi(ii,:)=edgephi(ii,:).*tempedgephi;
A=[SrcNodeIndex,DstNodeIndex];
[A,index,~]=unique(A,"rows");
SrcNodeIndex=A(:,1);
DstNodeIndex=A(:,2);
edgephi=edgephi(index,:);
