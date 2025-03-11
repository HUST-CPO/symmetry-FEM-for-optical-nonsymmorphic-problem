clc
clear all

%初始化
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
c_const=299792458;%真空光速
lam0=c_const/3e10;%波长
% ff0=c_const/lam0;
% lam0=c_const/ff0;
k0=2*pi/lam0;

% load('Data_slurm_total_test.mat')
Time1=[];
Time2=[];
nn=[];
eigen=[];
%% 
% ni =0:0.1:3;
ni =2.05;
%边排序
el2no=Mesh.Elements;
n1=el2no([1 1 1 2 2 3],:);
n2=el2no([2 3 4 3 4 4],:);
el_ed2no_array=[n1(:) n2(:)];
%按1-2 1-3 1-4 2-3 2-4 3-4排序
%Edges 全局边n*2  EdgesOfElements四面体边索引 6N*1 sign四面体边方向 6N*1
[Mesh.Edge,~,Mesh.EdgesOfElements]=unique(el_ed2no_array,'rows');
NbrEdges=max(size(Mesh.Edge));%边数
targetflag1=[1,5,9,16,21];
targetflag2=[88,89,90,91,92];
find1=findTris(targetflag1,Mesh);
find2=findTris(targetflag2,Mesh);
targetflag1=[2,6,31,48,72];
targetflag2=[28,29,33,52,75];
find3=findTris(targetflag1,Mesh);
find4=findTris(targetflag2,Mesh);
targetflag1=[3];
targetflag2=[8];
find5=findTris(targetflag1,Mesh);
find6=findTris(targetflag2,Mesh);
targetflag1=1:1:92;
index1=[11,12,23,24,42,43,73,74,1,5,9,16,21,88,89,90,91,92,2,6,31,48,72,28,29,33,52,75,3,8];
targetflag1(index1)=[];
find7=findTris(targetflag1,Mesh);


for ii = 1:length(ni)
    tic
    if (0<=ni(ii))&&(ni(ii)<=1)
        ns.u=(1-ni(ii))/2;ns.v=(1-ni(ii))/2;ns.w=0;
    elseif (1<ni(ii))&&(ni(ii)<=2)
        ns.u=(ni(ii)-1)/2;ns.v=0;ns.w=0;
    elseif (2<ni(ii))&&(ni(ii)<=3)
        ns.u=1/2;ns.v=(ni(ii)-2)/2;ns.w=0;
    end
    [copyOfBedge1,edgephi1]=GetcopyOfBedge1(find1,find2,Mesh,ns);
    [copyOfBedge2,edgephi2]=GetcopyOfBedge2(find3,find4,Mesh,ns);
    [copyOfBedge3,edgephi3]=GetcopyOfBedge3(find5,find6,Mesh,ns);
    DstNodeIndex = [copyOfBedge1(:,2);copyOfBedge2(:,2);copyOfBedge3(:,2)];
    SrcNodeIndex = [copyOfBedge1(:,1);copyOfBedge2(:,1);copyOfBedge3(:,1)];
    edgephi = [edgephi1;edgephi2;edgephi3];
    [SrcNodeIndex,DstNodeIndex,edgephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,edgephi,Mesh);
    DeleteIndex=unique(sort([DstNodeIndex;find7]));
    Index=1:(NbrEdges);
    Index(DeleteIndex)=[];
    P=speye(NbrEdges);
    for i = 1:size(DstNodeIndex,1)
        P(DstNodeIndex(i,1),SrcNodeIndex(i,1))=edgephi(i);
    end
    P=P(:,Index);
    %% 求解
    Time1=[Time1,toc];
    tic
    A_temp=P'*A*P;
    B_temp=P'*B*P;

    r=symrcm(B_temp);
    B_temp=B_temp(r,r);
    A_temp=A_temp(r,r);

    [~,D] = eigs(A_temp,B_temp,10,k0*k0);
    % [V,D] = eig(A,B);
    k_solve=diag(sqrt(D));

    lam_solve=2*pi./k_solve;
    ff=sort(abs(c_const./lam_solve));
    eigen=[eigen,ff];
    Time2=[Time2,toc];
    nn=[nn,ni(ii)];
    save('Data_total.mat','eigen','Time1','Time2','nn');
end
