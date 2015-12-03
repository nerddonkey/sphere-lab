tic
sims=2000;
n=80;
x=zeros(n*sims,1);
y=zeros(n*sims,1);
count=0;
scale=1 / sqrt(n);

for i = 1:sims
    [V,D] = eig(scale * randn(n, n));
    for j = 1:n
        count = count+1;
        x(count)=real(D(j,j));
        y(count)=imag(D(j,j));
    end
end

count
toc

plot(x,y,'r.')
axis('equal')
shg