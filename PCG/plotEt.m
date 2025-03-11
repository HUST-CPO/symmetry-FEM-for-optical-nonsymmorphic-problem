function  plotEt(node,elem,Ex,Ey,Ez,alpha)
%绘制Ex Ey
figure(2);
subplot(1,3,1);
h1=trisurf(elem,node(:,1),node(:,2),real(Ex));
shading interp
set(h1,'edgecolor','black','edgealpha',alpha);
view(2); axis equal; axis tight; axis off;
colorbar;
colormap jet;
caxis([min(real(Ex)) max(real(Ex))]);

subplot(1,3,2);
h2=trisurf(elem,node(:,1),node(:,2),real(Ey));
shading interp
set(h2,'edgecolor','black','edgealpha',alpha);
view(2); axis equal; axis tight; axis off;
colorbar;
colormap jet;
caxis([min(real(Ey)) max(real(Ey))]);

subplot(1,3,3);
h2=trisurf(elem,node(:,1),node(:,2),real(Ez));
shading interp
set(h2,'edgecolor','black','edgealpha',alpha);
view(2); axis equal; axis tight; axis off;
colorbar;
colormap jet;
caxis([min(real(Ez)) max(real(Ez))]);

figure(3);
normE=sqrt(abs(Ex).*abs(Ex)+abs(Ey).*abs(Ey)+abs(Ez).*abs(Ez));
h2=trisurf(elem,node(:,1),node(:,2),normE);
shading interp
set(h2,'edgecolor','black','edgealpha',alpha);
view(2); axis equal; axis tight; axis off;
colorbar;
colormap jet;
caxis([min(normE) max(normE)]);

end