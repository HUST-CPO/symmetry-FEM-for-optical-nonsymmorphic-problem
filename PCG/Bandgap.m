clc
clear all

load MeshData.mat
eigen=[];
Time1=[];
Time2=[];
ni = 0:0.1:3;

%数据
nbrVertex=Mesh.NbrNodes;


targetflag1=[7,5,3,1];
targetflag2=[19,20,21,22];
find1=findTris(targetflag1,Mesh);
find2=findTris(targetflag2,Mesh);
targetflag1=[2];
targetflag2=[9];
find3=findTris(targetflag1,Mesh);
find4=findTris(targetflag2,Mesh);

for ii = 1:length(ni)
    tic
    if (0<=ni(ii))&&(ni(ii)<=1)
        ns.u=ni(ii)/2;ns.v=0;ns.w=0;
        [~,copyOfBvertex1,~,nodephi1]=GetcopyOfBedge(find1,find2,Mesh,ns,1);
        [~,copyOfBvertex2,~,nodephi2]=GetcopyOfBedge2(find3,find4,Mesh,ns);
        DstNodeIndex = [copyOfBvertex1(:,2);copyOfBvertex2(:,2)];
        SrcNodeIndex = [copyOfBvertex1(:,1);copyOfBvertex2(:,1)];
        nodephi = [nodephi1;nodephi2];
        [SrcNodeIndex,DstNodeIndex,nodephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,nodephi,Mesh);
        martix_tgae=1;
    elseif (1<ni(ii))&&(ni(ii)<=2)
        ns.u=1/2;ns.v=(ni(ii)-1)/2;ns.w=0;
        [~,copyOfBvertex1,~,nodephi1]=GetcopyOfBedge(find1,find2,Mesh,ns,2);
        glideindex(:,1)=copyOfBvertex1(:,2);
        glideindex(:,2)=copyOfBvertex1(:,2)+Mesh.NbrNodes;
        [~,copyOfBvertex2,~,nodephi2]=GetcopyOfBedge2(find3,find4,Mesh,ns);
        DstNodeIndex = [copyOfBvertex1(:,2);copyOfBvertex2(:,2)];
        SrcNodeIndex = [copyOfBvertex1(:,1);copyOfBvertex2(:,1)];
        nodephi = [nodephi1;nodephi2];
        [SrcNodeIndex,DstNodeIndex,nodephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,nodephi,Mesh);
        [SrcNodeIndex,DstNodeIndex,nodephi]=glidechange(SrcNodeIndex,DstNodeIndex,Mesh.NbrNodes,glideindex,nodephi);
        martix_tgae=2;
    elseif (2<ni(ii))&&(ni(ii)<=3)
        ns.u=(3-ni(ii))/2;ns.v=1/2;ns.w=0;
        [~,copyOfBvertex1,~,nodephi1]=GetcopyOfBedge(find1,find2,Mesh,ns,3);
        [~,copyOfBvertex2,~,nodephi2]=GetcopyOfBedge2(find3,find4,Mesh,ns);
        DstNodeIndex = [copyOfBvertex1(:,2);copyOfBvertex2(:,2)];
        SrcNodeIndex = [copyOfBvertex1(:,1);copyOfBvertex2(:,1)];
        nodephi = [nodephi1;nodephi2];
        [SrcNodeIndex,DstNodeIndex,nodephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,nodephi,Mesh);
        martix_tgae=1;
    end
    tic;
     if martix_tgae==1
        DeleteIndex=unique(sort(DstNodeIndex));
        Index=1:(Mesh.NbrNodes);
        Index(DeleteIndex)=[];
        for j=1:2
            P=speye(Mesh.NbrNodes);
            for i = 1:size(DstNodeIndex,1)
                P(DstNodeIndex(i,1),SrcNodeIndex(i,1))=nodephi(i,j);
            end
            P=P(:,Index);
            eigen1((j-1)+1:j,1)=Assembel(Mesh,P);
        end
        eigen=[eigen,eigen1];
    elseif martix_tgae==2
        DeleteIndex=unique(sort(DstNodeIndex));
        Index=1:(Mesh.NbrNodes*2);
        Index(DeleteIndex)=[];
        P=speye(Mesh.NbrNodes*2);
        for i = 1:size(DstNodeIndex,1)
            P(DstNodeIndex(i,1),SrcNodeIndex(i,1))=nodephi(i);
        end
        P=P(:,Index);
        ff=Assembel2(Mesh,P);
        ff=[ff;ff];
        eigen=[eigen,ff];
    end
    Time2=[Time2,toc];
    nn=ni(ii);
    save('Data_2D_1300.mat','eigen','Time1','Time2','nn');
end