function [copyOfBedge,EdgePhi]=GetcopyOfBedge2(PBCEdge1,PBCEdge2,Mesh,ns)

dif0=[0,0;0,0];
rows = [1,3];

for i = 1:size(PBCEdge1,1)
    for j = 1:size(PBCEdge2,1)
        dif = Mesh.Nodes(:,Mesh.Edge(PBCEdge1(i),:))-Mesh.Nodes(:,Mesh.Edge(PBCEdge2(j),:));
        if(all(all(abs(dif(rows,:)-dif0)<0.000000001)))
            copyOfBedge(i,1)=PBCEdge1(i);
            copyOfBedge(i,2)=PBCEdge2(j);
            EdgePhi(i,:) = Getphi(0,1,0,ns.u,ns.v,ns.w);
        end
        dif = Mesh.Nodes(:,Mesh.Edge(PBCEdge1(i),:))-fliplr(Mesh.Nodes(:,Mesh.Edge(PBCEdge2(j),:)));
        if(all(all(abs(dif(rows,:)-dif0)<0.000000001)))
            copyOfBedge(i,1)=PBCEdge1(i);
            copyOfBedge(i,2)=PBCEdge2(j);
            EdgePhi(i,:) = -Getphi(0,1,0,ns.u,ns.v,ns.w);
        end
    end
end