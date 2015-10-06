function coastReg(minCoast)
if nargin<1
	minCoast=200;
end
load coast
if ~isnan(lat(end))
	lat(end+1)=nan;
	long(end+1)=nan;
end
idx2=find(isnan(lat))';
idx1=[1 idx2(1:end-1)+1];
n=numel(idx1);
reg_lat=cell(n,1);
reg_long=cell(n,1);
fprintf('\nMajor Regions and Index Range\n\n')
for k=1:n
	reg_lat{k}=lat(idx1(k):idx2(k)-1);
	reg_long{k}=long(idx1(k):idx2(k)-1);
	if (idx2(k)-idx1(k))>=minCoast
		plot(reg_long{k},reg_lat{k}); shg; pause(1)
			fprintf('Region %03d:	 %04d:%04d\n',k,idx1(k),idx2(k)-1)
	end
end
