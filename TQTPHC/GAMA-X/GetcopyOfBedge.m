function [copyOfBedge,copyOfBvertex,EdgePhi,NodePhi]=GetcopyOfBedge(PBCEdge1,PBCEdge2,Mesh,ns,tage)

dif0 = [-0.01;-0.01];
cols=2;

for i = 1:size(PBCEdge1,1)
    for j = 1:size(PBCEdge2,1)
        dif = Mesh.Nodes(Mesh.Edges(PBCEdge1(i),:),:)+Mesh.Nodes(Mesh.Edges(PBCEdge2(j),:),:);
       if(all(all(abs(dif(:,cols)-dif0)<0.000000001)))
            copyOfBedge(i,1)=PBCEdge1(i);
            copyOfBedge(i,2)=PBCEdge2(j);
            copyOfBvertex(i,1)=Mesh.Edges(PBCEdge1(i),1);
            copyOfBvertex(i,2)=Mesh.Edges(PBCEdge2(j),1);
            copyOfBvertex(size(PBCEdge1,1)+i,1)=Mesh.Edges(PBCEdge1(i),2);
            copyOfBvertex(size(PBCEdge1,1)+i,2)=Mesh.Edges(PBCEdge2(j),2);
            EdgePhi(i,:) = Getglidephi(ns.u,tage);
            NodePhi(i,:) = Getglidephi(ns.u,tage);
            NodePhi(size(PBCEdge1,1)+i,:) = Getglidephi(ns.u,tage);
        end
        dif = Mesh.Nodes(Mesh.Edges(PBCEdge1(i),:),:)+flipud(Mesh.Nodes(Mesh.Edges(PBCEdge2(j),:),:));
        if(all(all(abs(dif(:,cols)-dif0)<0.00000001)))
            copyOfBedge(i,1)=PBCEdge1(i);
            copyOfBedge(i,2)=PBCEdge2(j);
            copyOfBvertex(i,1)=Mesh.Edges(PBCEdge1(i),1);
            copyOfBvertex(i,2)=Mesh.Edges(PBCEdge2(j),2);
            copyOfBvertex(size(PBCEdge1,1)+i,1)=Mesh.Edges(PBCEdge1(i),2);
            copyOfBvertex(size(PBCEdge1,1)+i,2)=Mesh.Edges(PBCEdge2(j),1);
            EdgePhi(i,:) = -Getglidephi(ns.u,tage);
            NodePhi(i,:) = Getglidephi(ns.u,tage);
            NodePhi(size(PBCEdge1,1)+i,:) = Getglidephi(ns.u,tage);
        end
    end
end
[copyOfBvertex,index,~] = unique(copyOfBvertex,"rows","stable");
NodePhi=NodePhi(index,:);