[thetaVec,phiVec,maskR,thetaR,phiR]=coastRegion(3977:4210,0.5,0.5);
close all
subplot(1,2,1)
spy(maskR)
axis tight
subplot(1,2,2)
plot(phiR,thetaR)
set(gca,'YDir','reverse');
axis equal
axis tight
shg

% Region 001:	 0001:0612
% Region 013:	 0741:2386
% Region 088:	 3977:4210
% Region 092:	 4281:4794
% Region 093:	 4796:6667
% Region 180:	 8296:8604
