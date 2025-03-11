function [Ex,Ey]=PostProce_GetExEy(node,elem,edgeOfTri,Et)

Elements = elem';
NbrElement=max(size(Elements));
%绘制 
u=[0,1,0];v=[0,0,1];
NbrNodes=length(node);
Ex=zeros(NbrNodes,1);
Ey=zeros(NbrNodes,1);
numEt=zeros(NbrNodes,1);
for n=1:NbrElement
    %节点坐标
    xx=node(Elements(:,n),1);
    yy=node(Elements(:,n),2);
    %雅可比矩阵
    Jac=zeros(3,3);
    Jac(1,1)=xx(2)-xx(1);Jac(1,2)=yy(2)-yy(1);
    Jac(2,1)=xx(3)-xx(1);Jac(2,2)=yy(3)-yy(1);Jac(3,3)=1;
    InvJac=inv(Jac);
    %基函数
    et=zeros(3,3,3);
    for i=1:3
       for j=1:3
            et(i,:,j)=InvJac*BF_Triangle(3,1,j,u(i),v(i));
       end
    end
    %数值解
    for i=1:3
        for j=1:3
            Ex(Elements(i,n))=Ex(Elements(i,n))+et(i,1,j)*Et(edgeOfTri(n,j));
            Ey(Elements(i,n))=Ey(Elements(i,n))+et(i,2,j)*Et(edgeOfTri(n,j));
        end
        numEt(Elements(i,n))=numEt(Elements(i,n))+1;
    end
end
Ex=Ex./numEt;
Ey=Ey./numEt;

end