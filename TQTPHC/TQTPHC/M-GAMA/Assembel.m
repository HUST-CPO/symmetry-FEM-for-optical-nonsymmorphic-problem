clc
clear all

load MeshData.mat

epsilonr=[1,1,6.2,6.2,6.2];
mur=[1,1,1,1,1];

%高斯积分点
[xt,yt,pt,IntegralOrder]=GetGuassPoints(2);

%矩阵初始化
Dof=Mesh.NbrNodes;
NbrNon=Mesh.NbrElement*3*3;
NumA=0;
Ai=zeros(NbrNon,1);
Aj=zeros(NbrNon,1);
Av=complex(zeros(NbrNon,1));
NumB=0;
Bi=zeros(NbrNon,1);
Bj=zeros(NbrNon,1);
Bv=complex(zeros(NbrNon,1));

%循环网格 矩阵组装
for n=1:Mesh.NbrElement
    %节点坐标
    %节点坐标
    xx=Mesh.Nodes(Mesh.Elements(n,:),1);
    yy=Mesh.Nodes(Mesh.Elements(n,:),2);

    %雅可比矩阵
    Jac=zeros(3,3);
    Jac(1,1)=xx(2)-xx(1);Jac(1,2)=yy(2)-yy(1);
    Jac(2,1)=xx(3)-xx(1);Jac(2,2)=yy(3)-yy(1);Jac(3,3)=1;
    InvJac=inv(Jac);
    JacS(1,1)=InvJac(2,2);JacS(1,2)=-InvJac(2,1);
    JacS(2,1)=-InvJac(1,2);JacS(2,2)=InvJac(1,1);JacS(3,3)=1;
    DetJac=abs(det(Jac));
    DetJacx=det(Jac);
    TJac=Jac'/DetJacx;

    %basis functions
    Ez=zeros(3,3,IntegralOrder);curlEz=zeros(3,3,IntegralOrder);
    temp=zeros(3,1);
    for i=1:IntegralOrder
        for j=1:3
            temp(3)=BF_Ez(j,xt(i),yt(i));temp(1)=0;temp(2)=0;
            Ez(j,:,i)=temp;
            [temp(1),temp(2)]=BF_curlEz(j,xt(i),yt(i));temp(3)=0;
            curlEz(j,:,i)=JacS*temp;
        end
    end

    %材料参数
    domain=Mesh.Domains(n);
    mu=mur(domain);
    epsilon=epsilonr(domain);

    %子矩阵
    Sz=zeros(3,3);Tz=zeros(3,3);
    for i=1:3
        for j=1:3
            for k=1:IntegralOrder
                Sz(i,j)=Sz(i,j)+pt(k)*DetJac/mu*dot(curlEz(i,:,k),curlEz(j,:,k));
                Tz(i,j)=Tz(i,j)+pt(k)*DetJac/mu*epsilon*dot(Ez(i,:,k),Ez(j,:,k));
            end
        end
    end
    %存到三元组
    for i=1:3
        for j=1:3
            MappingIndexSi=Mesh.Elements(n,i);
            MappingIndexSj=Mesh.Elements(n,j);

            NumA=NumA+1;
            Ai(NumA)=MappingIndexSi;
            Aj(NumA)=MappingIndexSj;
            Av(NumA)=Sz(i,j);


            NumB=NumB+1;
            Bi(NumB)=MappingIndexSi;
            Bj(NumB)=MappingIndexSj;
            Bv(NumB)=Tz(i,j);
        end
    end
end

save('Ai.mat','Ai');
save('Aj.mat','Aj');
save('Av.mat','Av');
save('Bi.mat','Bi');
save('Bj.mat','Bj');
save('Bv.mat','Bv');

