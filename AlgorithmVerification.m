c= [1 1; 4 2; -1 2; -4 4; -3 -1; -1 -4; 3 -3; 4 -3.5;];  %centres of rbfs
G2d = @(x1,x2,y1,y2) exp(-((x1-x2).^2+(y1-y2).^2)./0.1); %2d Gaussian
rbf = @(x,y) [G2d(x,c(1,1),y,c(1,2)) G2d(x,c(2,1),y,c(2,2)) G2d(x,c(3,1),y,c(3,2)) G2d(x,c(4,1),y,c(4,2)) G2d(x,c(5,1),y,c(5,2)) G2d(x,c(6,1),y,c(6,2)) G2d(x,c(7,1),y,c(7,2)) G2d(x,c(8,1),y,c(8,2))];
a = [1 2 3 2.4 1.3 1.8 3.1 1.8]; %actual value of parameter
phi = @(x,y) (a*rbf(x,y)');      %Deansity function
k = 5;                           %P- Gain
ux = @(x1,x2) (k*(x2-x1));
uy = @(y1,y2) (k*(y2-y1));

L_dot = @(x,y) (rbf(x,y)'*rbf(x,y));  %Capital lambda
l_dot = @(x,y) (rbf(x,y)'*phi(x,y));  %Small lambda
a_cap = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]; %estimate of parameters
a_cap = a_cap';                            %to make it column vector
Nvert = size(c,1);                         %No. of gaussian centres
Y = eye(Nvert)*20;
disp(Nvert);
tvi = 1;         %target vertex index
x = 1;
y = 0;
delt = 0.01;    %Integration timestep
L = 0;          %Capital lambda  
l = 0;          %Small lambda
while(tvi <= Nvert)
    L = L + delt*L_dot(x,y);  
    l = l + delt*l_dot(x,y);
    b  = delt*Y*(L*a_cap-l);
    a_cap = a_cap - b;
    %disp(b);
    %disp(a_cap1)
    %disp("update");
    xt = c(tvi,1);    %Target co-ordinate
    yt = c(tvi,2);
    x=x+delt*ux(x,c(tvi,1));    %Position update
    y=y+delt*uy(y,c(tvi,2));
    disp((x-xt)^2 + (y-yt)^2);
    if(((x-xt)^2 + (y-yt)^2)<0.00005)
        tvi = tvi+1;
    end
    if(tvi == (Nvert+1))
        disp(a_cap);
    end    
end
