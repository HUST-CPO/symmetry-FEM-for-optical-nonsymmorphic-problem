function [copyOfBvertex,NodePhi]=GetcopyOfBedge(PBCEdge1,PBCEdge2,Mesh,ns)

dif0 = [0,-0.02;0,-0.02];
cols=[1,2];

for i = 1:size(PBCEdge1,1)
    for j = 1:size(PBCEdge2,1)
        dif = Mesh.Nodes(Mesh.Edges(PBCEdge1(i),:),:)-Mesh.Nodes(Mesh.Edges(PBCEdge2(j),:),:);
        if(all(all(abs(dif(:,cols)-dif0)<0.0001)))
            copyOfBvertex(i,1)=Mesh.Edges(PBCEdge1(i),1);
            copyOfBvertex(i,2)=Mesh.Edges(PBCEdge2(j),2);
            copyOfBvertex(size(PBCEdge1,1)+i,1)=Mesh.Edges(PBCEdge1(i),2);
            copyOfBvertex(size(PBCEdge1,1)+i,2)=Mesh.Edges(PBCEdge2(j),1);
            NodePhi(i,:) = Getphi(0,1,0,ns.u,ns.v,ns.w);
            NodePhi(size(PBCEdge1,1)+i,:) = Getphi(0,1,0,ns.u,ns.v,ns.w);
        end
        dif = Mesh.Nodes(Mesh.Edges(PBCEdge1(i),:),:)-flipud(Mesh.Nodes(Mesh.Edges(PBCEdge2(j),:),:));
        if(all(all(abs(dif(:,cols)-dif0)<0.0001)))
            copyOfBvertex(i,1)=Mesh.Edges(PBCEdge1(i),1);
            copyOfBvertex(i,2)=Mesh.Edges(PBCEdge2(j),2);
            copyOfBvertex(size(PBCEdge1,1)+i,1)=Mesh.Edges(PBCEdge1(i),2);
            copyOfBvertex(size(PBCEdge1,1)+i,2)=Mesh.Edges(PBCEdge2(j),1);
            NodePhi(i,:) = Getphi(0,1,0,ns.u,ns.v,ns.w);
            NodePhi(size(PBCEdge1,1)+i,:) = Getphi(0,1,0,ns.u,ns.v,ns.w);
        end
    end
end
[copyOfBvertex,index,~] = unique(copyOfBvertex,"rows","stable");
NodePhi=NodePhi(index,:);