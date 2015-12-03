tv = linspace(0.25, 1.5, 3);
pv = linspace(0.3, 1.2, 7);
tic;
Ylm=sphHarmGrid(15,2,tv,pv);
toc;
Ylm'
