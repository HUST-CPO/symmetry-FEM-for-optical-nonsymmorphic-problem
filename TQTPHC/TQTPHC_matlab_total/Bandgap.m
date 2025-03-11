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
ni = 0:0.1:3;

%数据
nbrVertex=Mesh.NbrNodes;


targetflag1=2;
targetflag2=3;
find1=findTris(targetflag1,Mesh);
find2=findTris(targetflag2,Mesh);
targetflag1=1;
targetflag2=4;
find3=findTris(targetflag1,Mesh);
find4=findTris(targetflag2,Mesh);


for ii = 1:length(ni)
    if (0<=ni(ii))&&(ni(ii)<=1)
        ns.u=ni(ii)/2;ns.v=0;ns.w=0;
    elseif (1<ni(ii))&&(ni(ii)<=2)
        ns.u=1/2;ns.v=(ni(ii)-1)/2;ns.w=0;
    elseif (2<ni(ii))&&(ni(ii)<=3)
        ns.u=(3-ni(ii))/2;ns.v=(3-ni(ii))/2;ns.w=0;
    end
    [copyOfBvertex1,nodephi1]=GetcopyOfBedge(find1,find2,Mesh,ns);
    [copyOfBvertex2,nodephi2]=GetcopyOfBedge2(find3,find4,Mesh,ns);
    DstNodeIndex = [copyOfBvertex1(:,2);copyOfBvertex2(:,2)];
    SrcNodeIndex = [copyOfBvertex1(:,1);copyOfBvertex2(:,1)];
    nodephi = [nodephi1;nodephi2(:,:)];
    [SrcNodeIndex,DstNodeIndex,nodephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,nodephi,Mesh);
    DeleteIndex=unique(sort(DstNodeIndex));
    Index=1:(Mesh.NbrNodes);
    Index(DeleteIndex)=[];
    tic
    P=speye(Mesh.NbrNodes);
    for i = 1:size(DstNodeIndex,1)
        P(DstNodeIndex(i,1),SrcNodeIndex(i,1))=nodephi(i);
    end
    P=P(:,Index);

    temp_A=P'*A*P;
    temp_B=P'*B*P;

    [~,D] = eigs(temp_A,temp_B,8,k0*k0);

    k00=diag(sqrt(D));
    lam_solve=2*pi./k00;
    ff=sort(abs(c_const./lam_solve));
    eigen=[eigen,ff];
    Time2=[Time2,toc];
    nn=ni(ii);
end
save('Data_total.mat','eigen','Time1','Time2','nn');