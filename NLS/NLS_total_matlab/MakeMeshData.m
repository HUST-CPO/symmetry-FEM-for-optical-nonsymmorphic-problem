clc
clear all

model=mphload('NLS_total.mph');
[m1, m2]=mphmeshstats(model);
Mesh.NbrNodes=max(size(m2.vertex));%顶点数
Mesh.NbrElements=max(size(m2.elem{2}));%四面体数
Mesh.Nodes=m2.vertex; %顶点坐标 3*n
Mesh.Domains=m2.elementity{2};%四面体区域索引 n*1
Mesh.Elements=m2.elem{2}+1;%四面体顶点索引 4*n
Mesh.Tris=m2.elem{3}+1;%边界面
Mesh.Trisflag=m2.elementity{3};%边界面所在区域
Mesh.Elements = sort(Mesh.Elements);

save('MeshData.mat','Mesh');