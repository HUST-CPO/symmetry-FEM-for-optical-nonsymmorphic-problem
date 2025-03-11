function plotTri(node,elem)
%绘制网格
figure(1);
subplot(1,1,1);
h=trisurf(elem,node(:,1),node(:,2),zeros(size(node,1),1));
set(h,'facecolor',[0.5 0.9 0.45],'edgecolor','black');
view(2); axis equal; axis tight; axis off;
end