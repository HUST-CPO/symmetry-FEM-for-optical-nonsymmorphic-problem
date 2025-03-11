function[SrcNodeIndex,DstNodeIndex,edgephi]=glidechange(SrcNodeIndex,DstNodeIndex,NbrEdges,glideindex,edgephi)
tempDstNodeIndex1=DstNodeIndex;
tempDstNodeIndex2=DstNodeIndex+NbrEdges;

tempSrcNodeIndex1=SrcNodeIndex;
tempSrcNodeIndex2=SrcNodeIndex+NbrEdges;

tempedgephi1=edgephi(:,1);
tempedgephi2=edgephi(:,2);
for i = 1:size(DstNodeIndex,1)
    for j = 1:size(glideindex,1)
        if tempDstNodeIndex1(i,1)==glideindex(j,1)
            tempDstNodeIndex1(i,1)=glideindex(j,2);
        end
    end
end
for i = 1:size(DstNodeIndex,1)
    for j = 1:size(glideindex,1)
        if tempDstNodeIndex2(i,1)==glideindex(j,2)
            tempDstNodeIndex2(i,1)=glideindex(j,1);
        end
    end
end
DstNodeIndex=[tempDstNodeIndex1;tempDstNodeIndex2];
SrcNodeIndex=[tempSrcNodeIndex1;tempSrcNodeIndex2];
edgephi=[tempedgephi1;tempedgephi2];