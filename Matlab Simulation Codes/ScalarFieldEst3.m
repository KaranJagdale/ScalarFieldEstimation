
c= [1 1; -4 4; -3 -1; 4 2; -1 2; -1 -4; 3 -3; 4 -3.5;];  %centres of rbfs
bot_locx = [-1 3 1.5 0];
bot_locy = [1.3 2 3 -2.5];
nG = length(c);
n = length(bot_locx);
edgeWeight = [];
for i = 1:n
    for j = 1:n
        if(adjacentvertex(bot_locx,bot_locy,i,j))
            %disp('edgeweight');
            
            %edgeWeight(i,j) = sqrt((bot_locx(i)-bot_locx(j))^2 + (bot_locy(i)-bot_locy(j))^2);
            edgeWeight(i,j) = 0;
            %disp(edgeWeight(i,j));
        else
            edgeWeight(i,j) = 0;
        end
    end
end    
for i = 1:n
    %[xbord1, ybord1]=compute_voronoi(i,[0 3 3 0], [0 0 3 3], bot_locx, bot_locy);
    %disp(i);
    %disp(xbord1);
    %disp(ybord1);
end 
[xbord1, ybord1] = compute_voronoi(1, [0 3 3 0], [0 0 3 3], bot_locx, bot_locy);
[xbord2, ybord2] = compute_voronoi(2, [0 3 3 0], [0 0 3 3], bot_locx, bot_locy);
%disp('Karan');
for i=1:length(xbord1)
   for j = 1:length(xbord2)
       xbord1(i) = round(10.^4.*xbord1(i))./(10.^4);
       xbord2(j) = round(10.^4.*xbord2(j))./(10.^4);
       if (xbord1(i) == xbord2(j))
           %disp('got');
           %disp(xbord1(i));
       end    
   end    
end 

%disp(intersect(xbord1,xbord2));
%disp('check');
adjacentvertex(bot_locx,bot_locy,1,2)
disp(edgeWeight);
a_cap = ones(n,nG)*0.1; 

c1 = GaussianDistribution(c,bot_locx,bot_locy,1);
c2 = GaussianDistribution(c,bot_locx,bot_locy,2);
c3 = GaussianDistribution(c,bot_locx,bot_locy,3);
c4 = GaussianDistribution(c,bot_locx,bot_locy,4);
disp('c1');
disp(c1);
%disp('c2');
disp(c2);
%disp('c3');
disp(c3);
%disp('c4')
disp(c4);
G2d = @(x1,x2,y1,y2) exp(-((x1-x2).^2+(y1-y2).^2)./0.1); %2d Gaussian
rbf = @(x,y) [G2d(x,c(1,1),y,c(1,2)) G2d(x,c(2,1),y,c(2,2)) G2d(x,c(3,1),y,c(3,2)) G2d(x,c(4,1),y,c(4,2)) G2d(x,c(5,1),y,c(5,2)) G2d(x,c(6,1),y,c(6,2)) G2d(x,c(7,1),y,c(7,2)) G2d(x,c(8,1),y,c(8,2))];
a = [1 2 3 2.4 1.3 1.8 3.1 1.8]; %actual value of parameter
phi = @(x,y) (a*rbf(x,y)');      %Deansity function
k = 5;                           %P- Gain
ux = @(x1,x2) (k*(x2-x1));
uy = @(y1,y2) (k*(y2-y1));

L_dot = @(x,y) (rbf(x,y)'*rbf(x,y));  %Capital lambda
l_dot = @(x,y) (rbf(x,y)'*phi(x,y));  %Small lambda
%a_cap = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]; %estimate of parameters
a_cap = a_cap';                            %to make it column vector
                        %No. of gaussian centres
Y = eye(nG)*50;

tvi = ones(n,1);         %target vertex index
x = bot_locx';
y = bot_locy';
delt = 0.001;    %Integration timestep
L = zeros(nG,nG,n);          %Capital lambda  
l = zeros(nG,n);          %Small lambda
X1 = [];
Y1 = [];
T = [];
t = 0;
loop= 1;
update = [];
epsilon = 5;
%perpbisector(1,2,3,4)

while(loop)
    L(:,:,1) = L(:,:,1) + delt*L_dot(x(1),y(1)); 
    L(:,:,2) = L(:,:,2) + delt*L_dot(x(2),y(2));
    L(:,:,3) = L(:,:,3) + delt*L_dot(x(3),y(3));
    L(:,:,4) = L(:,:,4) + delt*L_dot(x(4),y(4));
    l(:,1) = l(:,1) + delt*l_dot(x(1),y(1));
    l(:,1) = l(:,2) + delt*l_dot(x(2),y(2));
    l(:,1) = l(:,3) + delt*l_dot(x(3),y(3));
    l(:,1) = l(:,4) + delt*l_dot(x(4),y(4));
    %disp('1');
    %disp(delt*Y*((L(:,:,1)*a_cap(:,1)-l(:,1))));
    %disp('2');
    %disp(epsilon*consensus(a_cap,edgeWeight,1));
    update1 = delt*Y*((L(:,:,1)*a_cap(:,1)-l(:,1))+epsilon*consensus(a_cap,edgeWeight,1));
    a_cap(:,1) = a_cap(:,1) - update1;
    update2 = delt*Y*((L(:,:,2)*a_cap(:,2)-l(:,2))+epsilon*consensus(a_cap,edgeWeight,2));
    a_cap(:,2) = a_cap(:,2) - update2;
    update3 = delt*Y*((L(:,:,3)*a_cap(:,3)-l(:,3))+epsilon*consensus(a_cap,edgeWeight,3));
    a_cap(:,3) = a_cap(:,3) - update3;
    update4 = delt*Y*((L(:,:,4)*a_cap(:,4)-l(:,4))+epsilon*consensus(a_cap,edgeWeight,4));
    a_cap(:,4) = a_cap(:,4) - update4;
    
    %disp(tvi);
    if(tvi(1)<=length(c1(:,1)))
        xt(1) = c1(tvi(1),1);
        yt(1) = c1(tvi(1),2);
    end
    if(tvi(2)<=length(c2(:,1)))
        xt(2) = c2(tvi(2),1);
        yt(2) = c2(tvi(2),2);
    end
    if(tvi(3)<=length(c3(:,1)))
        xt(3) = c3(tvi(3),1);
        yt(3) = c3(tvi(3),2);
    end
    if(tvi(4)<=length(c4(:,1)))
        xt(4) = c4(tvi(4),1);
        yt(4) = c4(tvi(4),2);
    end
    %xt = [c1(tvi(1),1) c2(tvi(2),1) c3(tvi(3),1) c4(tvi(4),1)]';    %Target co-ordinate
    %yt = [c1(tvi(1),2) c2(tvi(2),2) c3(tvi(3),2) c4(tvi(4),2)]';
    x=x+delt*ux(x,xt);    %Position update
    y=y+delt*uy(y,yt);
    t = t+delt;
    X1 = [X1 x];
    Y1 = [Y1 y];
    T = [T t];
    
   % disp((x-xt)^2 + (y-yt)^2);
    if((((x(1)-xt(1))^2 + (y(1)-yt(1))^2)<0.000005) )
        tvi(1) = tvi(1)+1;
        disp('1');
        disp(x(1));
        %disp(c1(tvi(1),1));
    end
    if((((x(2)-xt(2))^2 + (y(2)-yt(2))^2)<0.000005) )
        tvi(2) = tvi(2)+1;
        disp('2');
        disp((2));
        %disp(c2(tvi(2),1));
    end
    if((((x(3)-xt(3))^2 + (y(3)-yt(3))^2)<0.000005) )
        tvi(3) = tvi(3)+1;
        disp('3');
        disp(x(3));
        %disp(c3(tvi(3),1));
    end
    if((((x(4)-xt(4))^2 + (y(4)-yt(4))^2)<0.000005) )
        disp('karan');
        tvi(4) = tvi(4)+1;
        %disp(c4(tvi(4),1));
        disp(x(4));
    end
       
    if(t >= 0 && (tvi(1) >= (length(c1(:,1))+1)) && (tvi(2) >= (length(c2(:,1))+1)) && (tvi(3) >= (length(c3(:,1))+1)) && (tvi(4) >= (length(c4(:,1)))+1))
        disp(a_cap);
        disp(update1);
        disp(t);
        loop = 0;
        %plot(X1,Y1)
        %plot(T,X1)
    end    
end