clc
clear all

load MeshData.mat
Ai=load('Ai.mat');
Aj=load('Aj.mat');
Av=load('Av.mat');
Bi=load('Bi.mat');
Bj=load('Bj.mat');
Bv=load('Bv.mat');
Ai=Ai.Ai;Bi=Bi.Bi;Aj=Aj.Aj;Bj=Bj.Bj;Av=Av.Av;Bv=Bv.Bv;
A = sparse(Ai,Aj,Av);
B = sparse(Bi,Bj,Bv);

clear Ai Aj Av Bi Bj Bv 
c_const=299792458;%真空光速
lam0=c_const/1e8;%波长
ff0=c_const/lam0;
k0=2*pi/lam0;

eigen=[];
Time1=[];
Time2=[];
ni = 0:0.1:1;

%数据
nbrVertex=Mesh.NbrNodes;


targetflag1=1;
targetflag2=6;
find1=findTris(targetflag1,Mesh);
find2=findTris(targetflag2,Mesh);
targetflag1=2;
targetflag2=5;
find3=findTris(targetflag1,Mesh);
find4=findTris(targetflag2,Mesh);
targetflag1=3;
targetflag2=7;
find5=findTris(targetflag1,Mesh);
find6=findTris(targetflag2,Mesh);

for ii = 1:length(ni)
    if (0<=ni(ii))&&(ni(ii)<=1)
        ns.u=ni(ii)/2;ns.v=0;ns.w=0;
        [~,copyOfBvertex1,~,nodephi1]=GetcopyOfBedge(find1,find2,Mesh,ns,1);
        [~,copyOfBvertex2,~,nodephi2]=GetcopyOfBedge2(find3,find4,Mesh,ns);
        [~,copyOfBvertex3,~,nodephi3]=GetcopyOfBedge3(find5,find6,Mesh,ns,1);
        copyOfBvertex2(end-1,:)=[];
        nodephi2(end-1,:)=[];
        DstNodeIndex = [copyOfBvertex1(:,2);copyOfBvertex2(:,2);copyOfBvertex3(1:end-1,2)];
        SrcNodeIndex = [copyOfBvertex1(:,1);copyOfBvertex2(:,1);copyOfBvertex3(1:end-1,1)];
        nodephi = [nodephi1;nodephi2(:,:);nodephi3(1:end-1,:)];
        % [SrcNodeIndex,DstNodeIndex,nodephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,nodephi,Mesh);
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
    elseif (2<ni(ii))&&(ni(ii)<=3)
        ns.u=(3-ni(ii))/2;ns.v=(3-ni(ii))/2;ns.w=0;
        [~,copyOfBvertex1,~,nodephi1]=GetcopyOfBedge(find1,find2,Mesh,ns,3);
        [~,copyOfBvertex2,~,nodephi2]=GetcopyOfBedge2(find3,find4,Mesh,ns);
        [~,copyOfBvertex3,~,nodephi3]=GetcopyOfBedge3(find5,find6,Mesh,ns,3);
        copyOfBvertex2(end-1,:)=[];
        nodephi2(end-1,:)=[];
        DstNodeIndex = [copyOfBvertex1(:,2);copyOfBvertex2(:,2);copyOfBvertex3(1:end-1,2)];
        SrcNodeIndex = [copyOfBvertex1(:,1);copyOfBvertex2(:,1);copyOfBvertex3(1:end-1,1)];
        nodephi = [nodephi1;nodephi2(:,:);nodephi3(1:end-1,:)];
        % [SrcNodeIndex,DstNodeIndex,nodephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,nodephi,Mesh);
    end
    DeleteIndex=unique(sort(DstNodeIndex));
    Index=1:(Mesh.NbrNodes);
    Index(DeleteIndex)=[];
    tic;
    for j=1:2
        P=speye(Mesh.NbrNodes);
        for i = 1:size(DstNodeIndex,1)
            P(DstNodeIndex(i,1),SrcNodeIndex(i,1))=nodephi(i,j);
        end
        P=P(:,Index);
        temp_A=P'*A*P;
        temp_B=P'*B*P;

        [~,D] = eigs(temp_A,temp_B,4,k0*k0);

        k00=diag(sqrt(D));
        lam_solve=2*pi./k00;
        ff=sort(abs(c_const./lam_solve));

        eigen1(4*(j-1)+1:4*j,1)=ff;
    end
    eigen=[eigen,eigen1];
    Time2=[Time2,toc];
    nn=ni(ii);
    % save('Data_2D_1300.mat','eigen','Time1','Time2','nn');
end
save('Data1.mat','eigen','Time1','Time2','nn');