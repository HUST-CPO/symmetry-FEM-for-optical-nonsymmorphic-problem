function [copyOfBedge,EdgePhi]=GetcopyOfBedge1(PBCEdge1,PBCEdge2,Mesh,ns,tage)

dif0=[2.25e-3,2.25e-3;2.25e-3,2.25e-3;0,0];
rows = [1,2,3];

for i = 1:size(PBCEdge1,1)
    for j = 1:size(PBCEdge2,1)
        dif = [Mesh.Nodes(1,Mesh.Edge(PBCEdge1(i),:))+Mesh.Nodes(1,Mesh.Edge(PBCEdge2(j),:));Mesh.Nodes(2,Mesh.Edge(PBCEdge1(i),:))+Mesh.Nodes(2,Mesh.Edge(PBCEdge2(j),:));Mesh.Nodes(3,Mesh.Edge(PBCEdge1(i),:))-Mesh.Nodes(3,Mesh.Edge(PBCEdge2(j),:))];
        if(all(all(abs(dif(rows,:)-dif0)<0.000000001)))
            copyOfBedge(i,1)=PBCEdge1(i);
            copyOfBedge(i,2)=PBCEdge2(j);
            EdgePhi(i,:) = -Getglidephi(ns.v,tage);
        end
        dif = [Mesh.Nodes(1,Mesh.Edge(PBCEdge1(i),:))+fliplr(Mesh.Nodes(1,Mesh.Edge(PBCEdge2(j),:)));Mesh.Nodes(2,Mesh.Edge(PBCEdge1(i),:))+fliplr(Mesh.Nodes(2,Mesh.Edge(PBCEdge2(j),:)));Mesh.Nodes(3,Mesh.Edge(PBCEdge1(i),:))-fliplr(Mesh.Nodes(3,Mesh.Edge(PBCEdge2(j),:)))];
        if(all(all(abs(dif(rows,:)-dif0)<0.000000001)))
            copyOfBedge(i,1)=PBCEdge1(i);
            copyOfBedge(i,2)=PBCEdge2(j);
            EdgePhi(i,:) = Getglidephi(ns.v,tage);
        end
    end
end