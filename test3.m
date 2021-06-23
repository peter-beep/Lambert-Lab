syms x 
syms t
func_test_3_1 = @(x,T) ( 1.6./(sqrt(pi).* x .* T) .* exp(- 1./(T) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* T) .* 2./sqrt(pi) .* int( exp(-t.^2),0, x./0.4  .*(4.9-5))));
func_test_2_3 = @(x,T) ( -1.6./(sqrt(pi).* x.^2 .* T) .* exp(-1./(T) .* (4.9-5)^2) + ( 0.8./( x.^3 .* T) .* 2./sqrt(pi) .* int(exp(-t.^2) , 0, x./0.4 .* (4.9-5))) - ( 0.8./(x.^2 .* T ) .* exp((-x./0.4 .* (4.9-5).^2)) ));

% ./func_test_2_3(ycfgiba(i),oi)

% initialize height parameters for each base of the target sequence
ycgfiba = [ ];
ycfgiba(1)=-5.000000000000000000000000000000000000000000003;

figure(23)
for oi = 20000 : 30000 : (20000+(30000*100))
    temp = 0;
for i = 1:20
    if (1 <= i) && (i<= 8)
    ycfgiba(i+1) = ycfgiba(i) - func_test_3_1(ycfgiba(i),oi)/func_test_2_3(ycfgiba(i),oi);
    ycfgiba_2 = ycfgiba(i+1);
    disp(double(ycfgiba_2))
    P1 = [i-1, temp, oi*(i-1)];
    P2 = [i,ycfgiba_2,oi];
    pts = [P1;P2];
    plot3(pts(:,1),pts(:,2),pts(:,3),'-')
    temp = ycfgiba_2;
    hold on;
    elseif (9<=i) && (i<=17)
    ycfgiba(i+1) = ycfgiba(i) - func_test_3_1(ycfgiba(i)/10,oi)/func_test_2_3(ycfgiba(i)/10,oi);
    ycfgiba_2 = ycfgiba(i+1);
    disp(double(ycfgiba_2))
    P1 = [i-1, temp, oi*(i-1)];
    P2 = [i,ycfgiba_2,oi];
    pts = [P1;P2];
    plot3(pts(:,1),pts(:,2),pts(:,3),'-')
    temp = ycfgiba_2;
    hold on;
    elseif (18<=i) && (i<=20)
    ycfgiba(i+1) = ycfgiba(i) - func_test_3_1(ycfgiba(i)/100,oi)/func_test_2_3(ycfgiba(i)/100,oi);
    ycfgiba_2 = ycfgiba(i+1);
    disp(double(ycfgiba_2))
    P1 = [i-1, temp, oi*(i-1)];
    P2 = [i,double(ycfgiba(i+1)),oi];
    pts = [P1;P2];
    plot3(pts(:,1),pts(:,2),pts(:,3),'-')
    temp = ycfgiba(i+1);
    hold on;   
    end 
end 

end 


ycfgib = [ ];
% 250
ycfgib(1)=-15.000000000000000000000000000000000000000000003;
syms x t T


figure(24)
temp_2=0;
for oa = 2000 : 3500 : (2000 + (3500*80))
for i = 1:20
    if (1 <= i) && (i<= 12)
    ycfgib(i+1) = ycfgib(i) - func_test_3_1(ycfgib(i),oa)/func_test_2_3(ycfgib(i),oa);
    ycfgib_2 = ycfgib(i+1);
    disp(ycfgib(i+1))
    P1 = [i-1, temp, oi*(i-1)];
    P2 = [i,double(ycfgib_2),oi];
    pts = [P1;P2];
    plot3(pts(:,1),pts(:,2),pts(:,3), '-')
    temp_2 = ycgfib(i+1);
    hold on;
    elseif (9<=i) && (i<=17)
    ycfgib(i+1) = ycfgib(i) - func_test_3_1(ycfgib(i)/500,oa)/func_test_2_3(ycfgib(i)/500,oa);
    ycfgib_2 = ycfgib(i+1);
    disp(double(ycfgib_2))
    P1 = [i-1, temp, oi*(i-1)];
    P2 = [i,double(ycfgib_2),oi];
    pts = [P1;P2];
    plot3(pts(:,1),pts(:,2),pts(:,3), '-')
    temp_2 = ycgfib(i+1);
    hold on;  
    elseif (13<=i) && (i<=20)
    ycfgib(i+1) = ycfgib(i) - func_test_3_1(ycfgib(i)/5000,oa)/func_test_2_3(ycfgib(i)/5000,oa);
    ycfgib_2 = ycfgiba(i+1);
    disp(ycfgib(i+1))
    P1 = [i-1, temp, oi*(i-1)];
    P2 = [i,double(ycfgib_2),oi];
    pts = [P1;P2];
    plot3(pts(:,1),pts(:,2),pts(:,3), '-')
    temp_2 = ycgfib(i+1);
    hold on;
    end 
end 

end 



