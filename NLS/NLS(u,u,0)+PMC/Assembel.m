clc
clear all

load MeshData.mat
epsilonr=[2.1*(1-1i*tan(0.00028)),1,1,1,2.1*(1-1i*tan(0.00028)),1,1,2.1*(1-1i*tan(0.00028)),1,1];%相对介电常数
% epsilonr=[1,2.1*(1-1i*tan(0.00028))];
el2no=Mesh.Elements;
n1=el2no([1 1 1 2 2 3],:);
n2=el2no([2 3 4 3 4 4],:);
el_ed2no_array=[n1(:) n2(:)];
%按1-2 1-3 1-4 2-3 2-4 3-4排序
%Edges 全局边n*2  EdgesOfElements四面体边索引 6N*1 sign四面体边方向 6N*1
[Mesh.Edge,~,Mesh.EdgesOfElements]=unique(el_ed2no_array,'rows');

%四面体体积分 高斯积分节点
IntegralOrder=5;
xt = [0.25,0.166666666667,0.166666666667,0.166666666667,0.5];
yt = [0.25,0.166666666667,0.166666666667,0.5,0.166666666667];
zt = [0.25,0.166666666667,0.5,0.166666666667,0.166666666667];
pt = [-0.133333333333,0.075,0.075,0.075,0.075];
% IntegralOrder=15;
% xt = [0.25,...
%     0.0919710780526,0.0919710780526,0.0919710780526,0.724086765841,...
%     0.319793627829,0.319793627829,0.319793627829,0.0406191165118,...
%     0.0563508326895,0.0563508326895,0.44364916731,0.0563508326895,...
%     0.44364916731,0.44364916731];
% yt = [0.25,...
%     0.0919710780526,0.0919710780526,0.724086765841,0.0919710780526,...
%     0.319793627829,0.319793627829,0.0406191165118,0.319793627829,...
%     0.0563508326895,0.44364916731,0.0563508326895,0.44364916731,...
%     0.0563508326895,0.44364916731];
% zt = [0.25,...
%     0.0919710780526,0.724086765841,0.0919710780526,0.0919710780526,...
%     0.319793627829,0.0406191165118,0.319793627829,0.319793627829,...
%     0.44364916731,0.0563508326895,0.0563508326895,0.44364916731,...
%     0.44364916731,0.0563508326895];
% pt = [0.0197530864198,...
%     0.0119895139632,0.0119895139632,0.0119895139632,0.0119895139632,...
%     0.011511367871,0.011511367871,0.011511367871,0.011511367871,...
%     0.00881834215168,0.00881834215168,0.00881834215168,0.00881834215168,...
%     0.00881834215168,0.00881834215168];

%矩阵初始化
NbrNon=Mesh.NbrElements*6*6;
NumA=0;
Ai=zeros(NbrNon,1);
Aj=zeros(NbrNon,1);
Av=complex(zeros(NbrNon,1));
NumB=0;
Bi=zeros(NbrNon,1);
Bj=zeros(NbrNon,1);
Bv=complex(zeros(NbrNon,1));

for n=1:Mesh.NbrElements
    %网格坐标
    x=zeros(4,1);y=zeros(4,1);z=zeros(4,1);
    x(1)=Mesh.Nodes(1,Mesh.Elements(1,n));y(1)=Mesh.Nodes(2,Mesh.Elements(1,n));z(1)=Mesh.Nodes(3,Mesh.Elements(1,n));
    x(2)=Mesh.Nodes(1,Mesh.Elements(2,n));y(2)=Mesh.Nodes(2,Mesh.Elements(2,n));z(2)=Mesh.Nodes(3,Mesh.Elements(2,n));
    x(3)=Mesh.Nodes(1,Mesh.Elements(3,n));y(3)=Mesh.Nodes(2,Mesh.Elements(3,n));z(3)=Mesh.Nodes(3,Mesh.Elements(3,n));
    x(4)=Mesh.Nodes(1,Mesh.Elements(4,n));y(4)=Mesh.Nodes(2,Mesh.Elements(4,n));z(4)=Mesh.Nodes(3,Mesh.Elements(4,n));
    %边长度
    l=zeros(6,1);
    l(1)=sqrt((x(1)-x(2))*(x(1)-x(2))+(y(1)-y(2))*(y(1)-y(2))+(z(1)-z(2))*(z(1)-z(2)));
    l(2)=sqrt((x(1)-x(3))*(x(1)-x(3))+(y(1)-y(3))*(y(1)-y(3))+(z(1)-z(3))*(z(1)-z(3)));
    l(3)=sqrt((x(1)-x(4))*(x(1)-x(4))+(y(1)-y(4))*(y(1)-y(4))+(z(1)-z(4))*(z(1)-z(4)));
    l(4)=sqrt((x(2)-x(3))*(x(2)-x(3))+(y(2)-y(3))*(y(2)-y(3))+(z(2)-z(3))*(z(2)-z(3)));
    l(5)=sqrt((x(2)-x(4))*(x(2)-x(4))+(y(2)-y(4))*(y(2)-y(4))+(z(2)-z(4))*(z(2)-z(4)));
    l(6)=sqrt((x(3)-x(4))*(x(3)-x(4))+(y(3)-y(4))*(y(3)-y(4))+(z(3)-z(4))*(z(3)-z(4)));
    %雅可比矩阵
    Jac=zeros(3,3);
    Jac(1,1)=x(1)-x(4);Jac(1,2)=y(1)-y(4);Jac(1,3)=z(1)-z(4);
    Jac(2,1)=x(2)-x(4);Jac(2,2)=y(2)-y(4);Jac(2,3)=z(2)-z(4);
    Jac(3,1)=x(3)-x(4);Jac(3,2)=y(3)-y(4);Jac(3,3)=z(3)-z(4);
    %     InvJac=inv(Jac);
    DetJac=abs(det(Jac));
    TJac=Jac'/DetJac;

    tempE=zeros(6,3,IntegralOrder);E=zeros(6,3,IntegralOrder);curlE=zeros(6,3,IntegralOrder);
    %计算高斯积分点对应的基函数值
    for i=1:IntegralOrder
        tempE(1,:,i)=BF_E4(xt(i),yt(i),zt(i))*l(1);tempE(4,:,i)=BF_E6(xt(i),yt(i),zt(i))*l(4);
        tempE(2,:,i)=BF_E5(xt(i),yt(i),zt(i))*l(2);tempE(5,:,i)=-BF_E2(xt(i),yt(i),zt(i))*l(5);
        tempE(3,:,i)=-BF_E1(xt(i),yt(i),zt(i))*l(3);tempE(6,:,i)=-BF_E3(xt(i),yt(i),zt(i))*l(6);
        E(1,:,i)=(Jac\tempE(1,:,i)')';E(4,:,i)=(Jac\tempE(4,:,i)')';
        E(2,:,i)=(Jac\tempE(2,:,i)')';E(5,:,i)=(Jac\tempE(5,:,i)')';
        E(3,:,i)=(Jac\tempE(3,:,i)')';E(6,:,i)=(Jac\tempE(6,:,i)')';

        tempE(1,:,i)=BF_CurlE4(xt(i),yt(i),zt(i))*l(1);tempE(4,:,i)=BF_CurlE6(xt(i),yt(i),zt(i))*l(4);
        tempE(2,:,i)=BF_CurlE5(xt(i),yt(i),zt(i))*l(2);tempE(5,:,i)=-BF_CurlE2(xt(i),yt(i),zt(i))*l(5);
        tempE(3,:,i)=-BF_CurlE1(xt(i),yt(i),zt(i))*l(3);tempE(6,:,i)=-BF_CurlE3(xt(i),yt(i),zt(i))*l(6);
        curlE(1,:,i)=(TJac*tempE(1,:,i)')';curlE(4,:,i)=(TJac*tempE(4,:,i)')';
        curlE(2,:,i)=(TJac*tempE(2,:,i)')';curlE(5,:,i)=(TJac*tempE(5,:,i)')';
        curlE(3,:,i)=(TJac*tempE(3,:,i)')';curlE(6,:,i)=(TJac*tempE(6,:,i)')';
    end
    %双循环计算子矩阵
    domain=Mesh.Domains(n);
    epsilon=epsilonr(domain);


    Ae=zeros(6,6);Be=zeros(6,6);
    for i=1:6
        for j=1:6
            for k=1: IntegralOrder
                Ae(i,j)=Ae(i,j)+pt(k)*DetJac*(curlE(i,1,k)*curlE(j,1,k)+curlE(i,2,k)*curlE(j,2,k)+curlE(i,3,k)*curlE(j,3,k));
                Be(i,j)=Be(i,j)+pt(k)*DetJac*epsilon*(E(i,1,k)*E(j,1,k)+E(i,2,k)*E(j,2,k)+E(i,3,k)*E(j,3,k));
            end
        end
    end
    %矩阵组装
    for i=1:6
        for j=1:6
            MappingIndexSi=Mesh.EdgesOfElements((n-1)*6+i);
            MappingIndexSj=Mesh.EdgesOfElements((n-1)*6+j);

            NumA=NumA+1;
            Ai(NumA)=MappingIndexSi;
            Aj(NumA)=MappingIndexSj;
            Av(NumA)=Ae(i,j);

            NumB=NumB+1;
            Bi(NumB)=MappingIndexSi;
            Bj(NumB)=MappingIndexSj;
            Bv(NumB)=Be(i,j);
        end
    end
end

save('A1i.mat','Ai');
save('A1j.mat','Aj');
save('A1v.mat','Av');
save('B1i.mat','Bi');
save('B1j.mat','Bj');
save('B1v.mat','Bv');

A = sparse(Ai,Aj,Av);
B = sparse(Bi,Bj,Bv);

NbrEdges=max(size(Mesh.Edge));%边数
% targetflag1=[1,2,3,4,5,7,8,9,10,11,12,13,14,16,17,18,19,20,22,23,24,25,26,28,29,30,32,33,35,36,37,38,39,40,41,42,43,44,45,47,48,49,50,51,53,54,55,56,57,59,60,61,62,63,64,65,66,67,68,69,70,71];
% targetflag1=[1,2,3,4,5,12];
targetflag1=1:1:71;
index1=[6,15,21,27,31,34,46,52,58];
targetflag1(index1)=[];
find7=findTris(targetflag1,Mesh);
DeleteIndex=unique(sort(find7));
Index=1:(NbrEdges);
Index(DeleteIndex)=[];
P=speye(NbrEdges);
P=P(:,Index);

A=P'*A*P;
B=P'*B*P;


% r=symrcm(B);
% B=B(r,r);
% A=A(r,r);

c_const=299792458;%真空光速
lam0=c_const/1e11;%波长
k0=2*pi/lam0;


[~,D] = eigs(A,B,6,k0*k0);
k1=diag(sqrt(D));

lam_solve=2*pi./k1;
ff1=sort(abs(c_const./lam_solve));
ff2=sort((c_const./lam_solve));

% A1=P1'*A*P1;
% B1=P1'*B*P1;
% 
% r=symrcm(B1);
% B1=B1(r,r);
% A1=A1(r,r);
% 
% [~,D] = eigs(A1,B1,6,k0*k0);
% k1=diag(sqrt(D));
% 
% lam_solve=2*pi./k1;
% ff1=sort(abs(c_const./lam_solve));
% 
% A2=P2'*A*P2;
% B2=P2'*B*P2;
% 
% r=symrcm(B2);
% B2=B2(r,r);
% A2=A2(r,r);
% 
% [~,D] = eigs(A2,B2,6,k0*k0);
% % [V,D] = eig(A,B);
% k2=diag(sqrt(D));
% 
% lam_solve=2*pi./k2;
% ff2=sort(abs(c_const./lam_solve));
% ff=sort([ff1;ff2]);