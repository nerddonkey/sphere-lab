[thetaVec,phiVec,maskR,thetaR,phiR]=coastRegion(8296:8604,0.5,0.5);
subplot(1,2,1)
spy(maskR)
axis tight
subplot(1,2,2)
plot(phiR,thetaR)
set(gca,'YDir','reverse');
axis equal
axis tight
shg
