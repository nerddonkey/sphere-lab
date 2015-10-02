view(80,30)
for i=0:5:358
	%view(80,30)
	delete(findall(gcf,'Type','light')) % remove existing lights
	lightangle(i-45,0)
	%lightangle(i+135,-45)
	lighting flat;%gouraud
	%camorbit(5.0,0.0);
	drawnow;
end