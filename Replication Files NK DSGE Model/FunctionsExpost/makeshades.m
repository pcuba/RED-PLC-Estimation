function hf = makeshades(x,y,z,c1)

x = x(:)';y=y(:)';z=z(:)';
X1 = [x,fliplr(x)]; 
Y1 = [y,fliplr(z)]; 

t = 0.45; % set transparency [0, 1]


hf = fill(X1,Y1,c1);
set(hf,'FaceAlpha',t,'EdgeColor','none')

