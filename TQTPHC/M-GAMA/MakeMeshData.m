clc
clear all

model=mphload("TQTPHC_M_GAMA.mph");
[~, m2]=mphmeshstats(model);

%网格数据初步处理
Mesh.Nodes=m2.vertex';%节点坐标
Mesh.NbrNodes=length(Mesh.Nodes);
Elements=m2.elem{2}+1;
Mesh.Elements = sort(Elements)';%网格索引
Mesh.NbrElement=length(Mesh.Elements);
Mesh.Domains=m2.elementity{2};%区域
Mesh=GetEdge(Mesh);%网格边编码
Mesh.Edges=m2.elem{1}'+1;%边界边
Mesh.Edgesflag=m2.elementity{1};%边界标记
Mesh.coonOfEdges=findEdge(1:15,Mesh);
Mesh.domainOfEdges=Mesh.Domains(Mesh.coonOfEdges(:,1));


save('MeshData.mat','Mesh');