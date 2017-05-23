load coast2
coastline=[90-lat mod(long,360)]*pi/180;
dlmwrite('coastlinerad.txt',coastline,'precision','%20.13f','delimiter',' ')
