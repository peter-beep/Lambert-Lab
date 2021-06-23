y = [ ];
y(1)=4;
syms x t
% 500 #2 
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 50) .* exp(- 1./(50) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 50) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4 .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 50) .* exp(-1./(50).* (4.9-5)^2) + ( 0.8./( x.^3 .* 50) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 50 ) .* exp((-x./0.4 .* (4.9-5).^2)) ));

figure(6)

func_test_3 = @(x, T) (1.6./( sqrt(pi) .* x .* T) .* exp(-1./(T) .* (4.9-5)^2) - (0.8./( x.^2 .* T) .* 2./sqrt(pi) .* int( exp(-t.^2) , 0 , x./0.4 .*(4.9-5))));
func_test_4 = @(x, T) (-1.6./(sqrt(pi).* x.^2 .* T) .* exp(-1./(T) .* (4.9-5)^2) + (1.6./( x.^3 .* T) .* 2./sqrt(pi) .* int(exp(-t.^2) ,0, x./0.4 .* (4.9-5))) - (0.8./(x.^2 .*T .* exp(-x./0.4 .* (4.9-5).^2))));


for yi = 90:15:850
for i = 1:20
     y(i+1) = y(i) - func_test_3(y(i), yi)/func_test_4(y(i), yi);
    disp(y(i+1))
    plot(i,y(i+1),'.')
    hold on;
end 
end

for i = 1:20
    y(i+1) = y(i) - func_test(y(i))/func_test_2(y(i));
    disp(y(i+1))
    plot(i,y(i+1),'.')
    hold on;
end 

yz = [ ];
% 75
yz(1)=23;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 150) .* exp(- 1./(150) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 150) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4 .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 150) .* exp(-1./(150) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 150) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 150 ) .* exp((-x./0.4 .* (4.9-5)).^2) ));


for i = 1:20
    yz(i+1) = yz(i) - func_test(yz(i))/func_test_2(yz(i));
    disp(yz(i+1))
    plot(i,yz(i+1),'.')
    hold on;
end 

yx = [ ];
% 100 
yx(1)=-20.0000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 100) .* exp(- 1./(100) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 100) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4 .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 100) .* exp(-1./(100) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 100) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 100 ) .* exp((-x./0.4 .* (4.9-5)).^2) ));


for i = 1:20
    yx(i+1) = yx(i) - func_test(yx(i))/func_test_2(yx(i));
    disp(yx(i+1))
    plot(i,yx(i+1),'.')
    hold on;
end 

figure(10) 

yzs = [ ];
yzs(1) = -2.4564575525689;
syms x t T
func_test_3 = @(x, T) (1.6./( sqrt(pi) .* x .* T) .* exp(-1./(T) .* (4.9-5)^2) - (0.8./( x.^2 .* T) .* 2./sqrt(pi) .* int( exp(-t.^2) , 0 , x./0.4 .*(4.9-5))));
func_test_4 = @(x, T) (-1.6./(sqrt(pi).* x.^2 .* T) .* exp(-1./(T) .* (4.9-5)^2) + (1.6./( x.^3 .* T) .* 2./sqrt(pi) .* int(exp(-t.^2) ,0, x./0.4 .* (4.9-5))) - (0.8./(x.^2 .*T .* exp(-x./0.4 .* (4.9-5).^2))));


for y = 150:10:750

for i = 1: 20
    
    yzs(i+1) = yzs(i) - func_test_3(yzs(i),y)/func_test_4(yzs(i),y);
    plot(i,yzs(i+1),'.')
    hold on;
      
end 

end

yd= [ ]; 
% 350
yd(1)=-3.000000000000000000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 350) .* exp(- 1./(350) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 350) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 350) .* exp(-1./(350) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 350) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 350 ) .* exp((-x./0.4 .* (4.9-5).^2)) ));


for i = 1:20
    yd(i+1) = yd(i) - func_test(yd(i))/func_test_2(yd(i));
    disp(yd(i+1))
    plot(i,yd(i+1),'.')
    hold on;
end 

ye = [ ];
% 400
ye(1)=-0.000000000000000000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 220) .* exp(- 1./(220) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 220) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 220) .* exp(-1./(220) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 220) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 220 ) .* exp((-x./0.4 .* (4.9-5).^2)) ));


for i = 1:20
    ye(i+1) = ye(i) - func_test(ye(i))/func_test_2(ye(i));
    disp(ye(i+1))
    plot(i,ye(i+1),'.')
    hold on;
end 

yf = [ ];
% 450
yf(1)=-7.000000000000000000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 450) .* exp(- 1./(450) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 450) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 450) .* exp(-1./(450) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 450) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 450 ) .* exp((-x./0.4 .* (4.9-5).^2)) ));


for i = 1:20
    yf(i+1) = yf(i) - func_test(yf(i))/func_test_2(yf(i));
    disp(yf(i+1))
    plot(i,yf(i+1),'.')
    hold on;
end 

figure(11)
yg = [ ];
% 450 #2
yg(1)=-8.0000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 550) .* exp(- 1./(550) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 550) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4 .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 550) .* exp(-1./(550) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 550) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 550 ) .* exp((-x./0.4 .* (4.9-5)).^2) ));

for i = 1:20
    yg(i+1) = yg(i) - func_test(yg(i))/func_test_2(yg(i));
    disp(yg(i+1))
    plot(i,yg(i+1),'.')
    hold on;
end 

yh = [ ];
% 450 #3
yh(1)=-5.0000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 550) .* exp(- 1./(550) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 550) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4 .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 550) .* exp(-1./(550) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 550) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 550 ) .* exp((-x./0.4 .* (4.9-5)).^2) ));


for i = 1:20
    yh(i+1) = yh(i) - func_test(yh(i))/func_test_2(yh(i));
    disp(yh(i+1))
    plot(i,yh(i+1),'.')
    hold on;
end 

yi = [ ];
% 300 #2
yi(1)=-15.0000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 400) .* exp(- 1./(400) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 400) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4 .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 400) .* exp(-1./(400) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 400) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 400 ) .* exp((-x./0.4 .* (4.9-5)).^2) ));

for i = 1:20
    yi(i+1) = yi(i) - func_test(yi(i))/func_test_2(yi(i));
    disp(yi(i+1))
    plot(i,yi(i+1),'.')
    hold on;
end 

figure(7)
ya = [ ];
% 150
ya(1)=-1.0000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 250) .* exp(- 1./(250) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 250) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4 .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 250) .* exp(-1./(250) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 250) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 250) .* exp((-x./0.4 .* (4.9-5)).^2) ));


for i = 1:20
    ya(i+1) = ya(i) - func_test(ya(i))/func_test_2(ya(i));
    disp(ya(i+1))
    plot(i,ya(i+1),'.')
    hold on;
end 

yb = [ ];
% 200
yb(1)=-10.000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 200) .* exp(- 1./(200) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 200) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 200) .* exp(-1./(200) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 200) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 200 ) .* exp((-x./0.4 .* (4.9-5)).^2) ));


for i = 1:20
    yb(i+1) = yb(i) - func_test(yb(i))/func_test_2(yb(i));
    disp(yb(i+1))
    plot(i,yb(i+1),'.')
    hold on;
end 

yc = [ ];
% 250
yc(1)=-5.000000000000000000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 175) .* exp(- 1./(175) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 175) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 175) .* exp(-1./(175) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 175) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 175 ) .* exp((-x./0.4 .* (4.9-5).^2)) ));


for i = 1:20
    yc(i+1) = yc(i) - func_test(yc(i))/func_test_2(yc(i));
    disp(yc(i+1))
    plot(i,yc(i+1),'.')
    hold on;
end 

ycf = [ ];
% 250
ycf(1)=-8.000000000000000000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 220) .* exp(- 1./(220) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 220) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 220) .* exp(-1./(220) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 220) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 220 ) .* exp((-x./0.4 .* (4.9-5).^2)) ));

for i = 1:20
    ycf(i+1) = ycf(i) - func_test(ycf(i))/func_test_2(ycf(i));
    disp(ycf(i+1))
    plot(i,ycf(i+1),'.')
    hold on;
end 

ycfg = [ ];
% 250
ycfg(1)=-8.000000000000000000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 260) .* exp(- 1./(260) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 260) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 260) .* exp(-1./(260) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 260) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 260 ) .* exp((-x./0.4 .* (4.9-5).^2)) ));

for i = 1:20
    ycfg(i+1) = ycfg(i) - func_test(ycfg(i))/func_test_2(ycfg(i));
    disp(ycfg(i+1))
    plot(i,ycfg(i+1),'.')
    hold on;
end 

ycfgi = [ ];
% 250
ycfgi(1)=-3.000000000000000000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 215) .* exp(- 1./(215) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 215) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 215) .* exp(-1./(215) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 215) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 215 ) .* exp((-x./0.4 .* (4.9-5).^2)) ));

for i = 1:20
    ycfgi(i+1) = ycfgi(i) - func_test(ycfgi(i))/func_test_2(ycfgi(i));
    disp(ycfgi(i+1))
    plot(i,ycfgi(i+1),'.')
    hold on;
end 

ycfgia = [ ];
% 250
ycfgia(1)=-23.000000000000000000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 215) .* exp(- 1./(215) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 215) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 215) .* exp(-1./(215) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 215) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 215 ) .* exp((-x./0.4 .* (4.9-5).^2)) ));

for i = 1:20
    ycfgia(i+1) = ycfgia(i) - func_test(ycfgia(i))/func_test_2(ycfgia(i));
    disp(ycfgia(i+1))
    plot(i,ycfgia(i+1),'.')
    hold on;
end 


ycfgi = [ ];
% 250
ycfgi(1)=-3.000000000000000000000000000000000000000000003;
syms x t
func_test = @(x) ( 1.6./(sqrt(pi).* x .* 50) .* exp(- 1./(50) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* 50) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2 = @(x) ( -1.6./(sqrt(pi).* x.^2 .* 215) .* exp(-1./(50) .* (4.9-5)^2) + ( 0.8./( x.^3 .* 50) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* 50 ) .* exp((-x./0.4 .* (4.9-5).^2)) ));

for i = 1:20
    ycfgi(i+1) = ycfgi(i) - func_test(ycfgi(i))/func_test_2(ycfgi(i));
    disp(ycfgi(i+1))
    plot(i,ycfgi(i+1),'.')
    hold on;
end 


ycfgiba = [ ];
% 250
ycfgiba(1)=-3.000000000000000000000000000000000000000000003;
syms x t T
func_test_3_1 = @(x,T) ( 1.6./(sqrt(pi).* x .* T) .* exp(- 1./(T) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* T) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2_3 = @(x,T) ( -1.6./(sqrt(pi).* x.^2 .* T) .* exp(-1./(T) .* (4.9-5)^2) + ( 0.8./( x.^3 .* T) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* T ) .* exp((-x./0.4 .* (4.9-5).^2)) ));
