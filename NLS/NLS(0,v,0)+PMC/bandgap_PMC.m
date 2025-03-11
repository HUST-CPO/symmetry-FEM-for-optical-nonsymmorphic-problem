clc
clear all

%初始化
load MeshData.mat
A1i=load('A1i.mat');
A1j=load('A1j.mat');
A1v=load('A1v.mat');
B1i=load('B1i.mat');
B1j=load('B1j.mat');
B1v=load('B1v.mat');
A1i=A1i.Ai;B1i=B1i.Bi;A1j=A1j.Aj;B1j=B1j.Bj;A1v=A1v.Av;B1v=B1v.Bv;
A = sparse(A1i,A1j,A1v);
B = sparse(B1i,B1j,B1v);

clear A1i A1j A1v B1i B1j B1v 

c_const=299792458;%真空光速
lam0=c_const/2.5e10;%波长
k0=2*pi/lam0;

% load('Data_matlab_test.mat')
Time1=[];
Time2=[];
nn=[];
eigen=[];

ni =1.1:0.1:2;
% ni =1.1:0.1:2;

%边排序
el2no=Mesh.Elements;
n1=el2no([1 1 1 2 2 3],:);
n2=el2no([2 3 4 3 4 4],:);
el_ed2no_array=[n1(:) n2(:)];
%按1-2 1-3 1-4 2-3 2-4 3-4排序
%Edges 全局边n*2  EdgesOfElements四面体边索引 6N*1 sign四面体边方向 6N*1
[Mesh.Edge,~,Mesh.EdgesOfElements]=unique(el_ed2no_array,'rows');
NbrEdges=max(size(Mesh.Edge));%边数
targetflag1=[36,52,65];
targetflag2=[49,50,63];
find1=findTris(targetflag1,Mesh);
find2=findTris(targetflag2,Mesh);
targetflag1=[1,4,10];
targetflag2=[69,70,71];
find3=findTris(targetflag1,Mesh);
find4=findTris(targetflag2,Mesh);
targetflag1=[2,5,29];
targetflag2=[16,18,31];
find5=findTris(targetflag1,Mesh);
find6=findTris(targetflag2,Mesh);
% targetflag1=[3,7,8,9,11,12,13,16,17,19,20,22,23,25,26,30,32,33,35,36,37,38,39,40,41,42,43,44,45,47,48,49,51,54,55,56,57,59,60,61,62,63,64,65];
targetflag1=1:1:71;
% index1=[6,13,28,35,39,42,45,48,62,36,52,65,49,50,63,1,4,10,69,70,71,2,5,29,16,18,31];
index1=[6,13,28,35,39,42,45,48,62,36,52,65,49,50,63,1,4,10,69,70,71,2,5,29,16,18,31,3,7,12,27,38,41,47,61];
targetflag1(index1)=[];
find7=findTris(targetflag1,Mesh);

for ii = 1:length(ni)
    tic
    if (0<=ni(ii))&&(ni(ii)<=1)
        ns.u=ni(ii)/2;ns.v=0;ns.w=0;
        [copyOfBedge1,edgephi1]=GetcopyOfBedge1(find1,find2,Mesh,ns,1);
        [copyOfBedge2,edgephi2]=GetcopyOfBedge2(find3,find4,Mesh,ns);
        [copyOfBedge3,edgephi3]=GetcopyOfBedge3(find5,find6,Mesh,ns,1);
        DstNodeIndex = [copyOfBedge1(:,2);copyOfBedge2(:,2);copyOfBedge3(:,2)];
        SrcNodeIndex = [copyOfBedge1(:,1);copyOfBedge2(:,1);copyOfBedge3(:,1)];
        edgephi = [edgephi1;edgephi2;edgephi3];
        [SrcNodeIndex,DstNodeIndex,edgephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,edgephi,Mesh);
        martix_tgae=1;
    elseif (1<ni(ii))&&(ni(ii)<=2)
        ns.u=1/2;ns.v=(ni(ii)-1)/2;ns.w=0;
        [copyOfBedge1,edgephi1]=GetcopyOfBedge1(find1,find2,Mesh,ns,2);
        [copyOfBedge2,edgephi2]=GetcopyOfBedge2(find3,find4,Mesh,ns);
        [copyOfBedge3,edgephi3]=GetcopyOfBedge3(find5,find6,Mesh,ns,2);
        DstNodeIndex = [copyOfBedge1(:,2);copyOfBedge2(:,2);copyOfBedge3(:,2)];
        SrcNodeIndex = [copyOfBedge1(:,1);copyOfBedge2(:,1);copyOfBedge3(:,1)];
        edgephi = [edgephi1;edgephi2;-edgephi3];
        [SrcNodeIndex,DstNodeIndex,edgephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,edgephi,Mesh);
        martix_tgae=1;
    elseif (2<ni(ii))&&(ni(ii)<=3)
        ns.u=0;ns.v=0;ns.w=(ni(ii)-2)/2;
        [copyOfBedge1,edgephi1]=GetcopyOfBedge1(find1,find2,Mesh,ns,3);
        glideindex(:,1)=copyOfBedge1(:,2);
        glideindex(:,2)=copyOfBedge1(:,2)+NbrEdges;
        [copyOfBedge2,edgephi2]=GetcopyOfBedge2(find3,find4,Mesh,ns);
        [copyOfBedge3,edgephi3]=GetcopyOfBedge3(find5,find6,Mesh,ns);
        DstNodeIndex = [copyOfBedge1(:,2);copyOfBedge2(:,2);copyOfBedge3(:,2)];
        SrcNodeIndex = [copyOfBedge1(:,1);copyOfBedge2(:,1);copyOfBedge3(:,1)];
        edgephi = [-edgephi1;edgephi2;edgephi3];
        [SrcNodeIndex,DstNodeIndex,edgephi]=doubleBloch(SrcNodeIndex,DstNodeIndex,edgephi,Mesh);
        [SrcNodeIndex,DstNodeIndex,edgephi]=glidechange(SrcNodeIndex,DstNodeIndex,NbrEdges,glideindex,edgephi);
        martix_tgae=2;
    end
    DeleteIndex=unique(sort([DstNodeIndex;find7]));
    Index=1:(NbrEdges);
    Index(DeleteIndex)=[];
    P1=speye(NbrEdges);
    P2=speye(NbrEdges);
    for i = 1:size(DstNodeIndex,1)
        P1(DstNodeIndex(i,1),SrcNodeIndex(i,1))=edgephi(i,1);
        P2(DstNodeIndex(i,1),SrcNodeIndex(i,1))=edgephi(i,2);
    end
    P1=P1(:,Index);
    P2=P2(:,Index);
    Time1=[Time1;toc];
    tic;
    A1=P1'*A*P1;
    B1=P1'*B*P1;

    r=symrcm(B1);
    B1=B1(r,r);
    A1=A1(r,r);

    [~,D] = eigs(A1,B1,3,k0*k0);
    k1=diag(sqrt(D));

    lam_solve=2*pi./k1;
    ff1=sort(abs(c_const./lam_solve));
    ff3=sort((c_const./lam_solve));

    A2=P2'*A*P2;
    B2=P2'*B*P2;

    r=symrcm(B2);
    B2=B2(r,r);
    A2=A2(r,r);

    [~,D] = eigs(A2,B2,3,k0*k0);
    k2=diag(sqrt(D));

    lam_solve=2*pi./k2;
    ff2=sort(abs(c_const./lam_solve));
    ff4=sort((c_const./lam_solve));
    ff=([ff1;ff2]);

    eigen=[eigen,ff];

    Time2=[Time2,toc];
    nn=[nn,ni(ii)];
    save('T1&T2','eigen','Time1','Time2',"nn");
end
