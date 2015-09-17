fig=gcf;
fig.Position(3)=450;
fig.Position(4)=450;

fa=fig.CurrentAxes;

delete(findall(gcf,'Type','light')) % remove existing lights
light; lightangle(260,-45) % add 2 lights
lighting gouraud % preferred lighting for a curved surface
view(-100,-30) % set viewpoint

maxa=1.0; % max of what radius above can be
axis([-maxa maxa -maxa maxa -maxa maxa]);
axis off

set(fa,'CameraViewAngle',9) % zoom into scene
axes(fa);
rotate3d on % setup for for gui interaction

shg;
