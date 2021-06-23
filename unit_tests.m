function unit_tests()

critical_vec = [1.0500e+21 1.1e+21 1.12e+21 1.13e+21 1.14e+21 1.15e+21 1.16e+21 1.17e+21 1.18e+21 1.19e+21 1.20e+21 1.21e+21 1.22e+21 1.23e+21 1.24e+21 1.25e+21 1.26e+21 1.27e+21 1.28e+21 1.29e+21 1.2955555e+21 1.3e+21 1.31e+21 1.32e+21 1.33e+21 1.34e+21 1.35e+21 1.36e+21 1.37e+21 1.38e+21 1.39e+21 1.4e+21 1.41e+21 1.42e+21 1.43e+21 1.44e+21 1.45e+21 1.46e+21 1.47e+21 1.48e+21 1.49e+21 1.5e+21 1.51e+21 1.52e+21 1.53e+21 1.54e+21 1.55e+21 1.56e+21 1.57e+21 1.58e+21 1.59e+21 1.6e+21 1.61e+21 1.62e+21 1.63e+21 1.64e+21 1.65e+21 1.66e+21 1.67e+21 1.68e+21 1.69e+21 1.7e+21 1.71e+21 1.72e+21 1.73e+21 1.74e+21 1.75e+21 1.76e+21 1.77e+21 1.78e+21 1.79e+21 1.8e+21 1.81e+21 1.82e+21 1.83e+21 1.84e+21 1.85e+21 1.86e+21 1.87e+21 1.88e+21 1.89e+21 1.9e+21 1.91e+21 1.92e+21 1.93e+21]
for i = 1 : length(critical_vec)
S=solve(1.6./(sqrt(pi).* x .* (1.3806e-23) .* critical_vec(i)) .* exp(- 1./((1.3806e-23) .* critical_vec(i)) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* (1.3806e-23).* critical_vec(i)) .* 2./sqrt(pi) .* int( exp(-t.^2),0, (x.*(1.3806e-23).* critical_vec(i))./0.4  .*(4.9-5)))==0.00000001,x)
S=double(S);
plot(critical_vec(i), S,'.')
hold on;
xlabel('Barrier height')
ylabel('Temperature')
end
critical_vec = [1.0500e+21 1.1e+21 1.12e+21 1.13e+21 1.14e+21 1.15e+21 1.16e+21 1.17e+21 1.18e+21 1.19e+21 1.20e+21 1.21e+21 1.22e+21 1.23e+21 1.24e+21 1.25e+21 1.26e+21 1.27e+21 1.28e+21 1.29e+21 1.2955555e+21 1.3e+21 1.31e+21 1.32e+21 1.33e+21 1.34e+21 1.35e+21 1.36e+21 1.37e+21 1.38e+21 1.39e+21 1.4e+21 1.41e+21 1.42e+21 1.43e+21 1.44e+21 1.45e+21 1.46e+21 1.47e+21 1.48e+21 1.49e+21 1.5e+21 1.51e+21 1.52e+21 1.53e+21 1.54e+21 1.55e+21 1.56e+21 1.57e+21 1.58e+21 1.59e+21 1.6e+21 1.61e+21 1.62e+21 1.63e+21 1.64e+21 1.65e+21 1.66e+21 1.67e+21 1.68e+21 1.69e+21 1.7e+21 1.71e+21 1.72e+21 1.73e+21 1.74e+21 1.75e+21 1.76e+21 1.77e+21 1.78e+21 1.79e+21 1.8e+21 1.81e+21 1.82e+21 1.83e+21 1.84e+21 1.85e+21 1.86e+21 1.87e+21 1.88e+21 1.89e+21 1.9e+21 1.91e+21 1.92e+21 1.93e+21 1.95e+21 1.96e+21 1.97e+21 1.98e+21 1.99e+21 1.0e+22 1.1e+22 1.2e+22 1.3e+22 1.4e+22 5.0e+23 5.1e+23 5.2e+23 5.3e+23 5.4e+23 5.5e+23 5.6e+23 5.7e+23 5.8e+23 6.5e+24 6.6e+24 5.5e+25 5.6e+25]
for i = 1 : length(critical_vec)
S=solve(1.6./(sqrt(pi).* x .* (1.3806e-23) .* critical_vec(i)) .* exp(- 1./((1.3806e-23) .* critical_vec(i)) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* (1.3806e-23).* critical_vec(i)) .* 2./sqrt(pi) .* int( exp(-t.^2),0, (x.*(1.3806e-23).* critical_vec(i))./0.4  .*(4.9-5)))==0.00000001,x)
S=double(S);
plot(critical_vec(i), S,'.')
hold on;
xlabel('Barrier height')
ylabel('Temperature')
end
critical_vec = [1.0500e+21 1.1e+21 1.12e+21 1.13e+21 1.14e+21 1.15e+21 1.16e+21 1.17e+21 1.18e+21 1.19e+21 1.20e+21 1.21e+21 1.22e+21 1.23e+21 1.24e+21 1.25e+21 1.26e+21 1.27e+21 1.28e+21 1.29e+21 1.2955555e+21 1.3e+21 1.31e+21 1.32e+21 1.33e+21 1.34e+21 1.35e+21 1.36e+21 1.37e+21 1.38e+21 1.39e+21 1.4e+21 1.41e+21 1.42e+21 1.43e+21 1.44e+21 1.45e+21 1.46e+21 1.47e+21 1.48e+21 1.49e+21 1.5e+21 1.51e+21 1.52e+21 1.53e+21 1.54e+21 1.55e+21 1.56e+21 1.57e+21 1.58e+21 1.59e+21 1.6e+21 1.61e+21 1.62e+21 1.63e+21 1.64e+21 1.65e+21 1.66e+21 1.67e+21 1.68e+21 1.69e+21 1.7e+21 1.71e+21 1.72e+21 1.73e+21 1.74e+21 1.75e+21 1.76e+21 1.77e+21 1.78e+21 1.79e+21 1.8e+21 1.81e+21 1.82e+21 1.83e+21 1.84e+21 1.85e+21 1.86e+21 1.87e+21 1.88e+21 1.89e+21 1.9e+21 1.91e+21 1.92e+21 1.93e+21 1.95e+21 1.96e+21 1.97e+21 1.98e+21 1.99e+21 1.0e+22 1.1e+22 1.2e+22 1.3e+22 1.4e+22 1.5e+22 1.6e+22 1.7e+22 1.8e+22 1.9e+22 2.1e+22 3.5e+22 3.7e+22 3.8e+22 3.9e+22 4.0e+22 4.8e+22 9.8e+23 1.0e+23 1.2e+23 1.3e+23 1.4e+23 1.5e+23 1.6e+23 1.7e+23 1.8e+23 1.9e+23 2.0e+23 2.1e+23 2.2e+23 2.4e+23 2.5e+23 2.6e+23 2.7e+23 2.8e+23 3.0e+23 3.1e+23 3.2e+23 3.3e+23 3.4e+23 3.6e+23 3.8e+23 5.0e+23 5.1e+23 5.2e+23 5.3e+23 5.4e+23 5.5e+23 5.6e+23 5.7e+23 5.8e+23 1e+24 1.1e+24 1.2e+24 1.3e+24 1.4e+24 1.5e+24 1.6e+24 1.7e+24 1.8e+24 1.9e+24 2.2e+24 2.3e+24 2.5e+24 2.6e+24 2.7e+24 2.8e+24 2.9e+24 3.0e+24 3.1e+24 3.2e+24 3.3e+24 3.5e+24 3.6e+24 3.7e+24 3.8e+24 3.9e+24 4.3e+24 4.4e+24 4.5e+24 4.6e+24 4.8e+24 5.0e+24 5.1e+24 5.2e+24 5.3e+24 5.4e+24 5.5e+24 5.6e+24 5.8e+24 5.9e+24 6.0e+24 6.0005e+24 6.1e+24 6.2e+24 6.3e+24 6.4e+24 6.5e+24 6.6e+24 5.5e+25 5.6e+25 1.0e+26 1.5e+26 2.0e+26 2.5e+26 3.0e+26 3.5e+26 4.0e+26 4.5e+26 5.0e+26 5.5e+26 1.0e+27 1.5e+27 2.0e+27 2.5e+27 3.0e+27 3.5e+27 4.0e+27 4.5e+27 5.0e+27]
for i = 1 : length(critical_vec)
S=solve(1.6./(sqrt(pi).* x .* (1.3806e-23) .* critical_vec(i)) .* exp(- 1./((1.3806e-23) .* critical_vec(i)) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* (1.3806e-23).* critical_vec(i)) .* 2./sqrt(pi) .* int( exp(-t.^2),0, (x.*(1.3806e-23).* critical_vec(i))./0.4  .*(4.9-5)))==0.00000001,x)
S=double(S);
plot(critical_vec(i), S,'.')
hold on;
xlabel('Barrier height')
ylabel('Temperature')
end
critical_vec_2 = [0.5e+22 0.51e+22 0.52e+22 0.53e+22 0.57e+22 0.7e+22 0.73e+22 0.8e+22 0.81e+22 0.82e+22 0.83e+22 0.84e+22 0.85e+22 0.9e+22 0.91e+22 0.92e+22 0.93e+22 0.94e+22 2.3e+22 2.4e+22 2.5e+22 2.6e+22]
for i = 1:length(critical_vec_2)
S=solve(1.6./(sqrt(pi).* x .* (1.3806e-23) .* critical_vec_2(i)) .* exp(- 1./((1.3806e-23) .* critical_vec_2(i)) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* (1.3806e-23).* critical_vec_2(i)) .* 2./sqrt(pi) .* int( exp(-t.^2),0, (x.*(1.3806e-23).* critical_vec_2(i))./0.4  .*(4.9-5)))==0.00000001,x)
S=double(S);plot(critical_vec_2(i),S,'.')
hold on;
end
critical_vec_3 = [0.2e+22 0.15e+22 0.18e+22 0.21e+22 0.2222111e+22 0.1947573e+22 0.5e+22 0.45e+22 0.345e+22]
for i =1:length(critical_vec_3)
S=solve(1.6./(sqrt(pi).* x .* (1.3806e-23) .* critical_vec_3(i)) .* exp(- 1./((1.3806e-23) .* critical_vec_3(i)) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* (1.3806e-23).* critical_vec_3(i)) .* 2./sqrt(pi) .* int( exp(-t.^2),0, (x.*(1.3806e-23).* critical_vec_3(i))./0.4  .*(4.9-5)))==0.00000001,x)
S=double(S);
plot(critical_vec_3(i),S,'.')
hold on;
end
critical_vec_4 = [1.51e+22 1.5000048602e+22 1.5784e+22 1.5362e+22 1.5154e+22 1.50486848e+22]
for i =1:length(critical_vec_4)
S=solve(1.6./(sqrt(pi).* x .* (1.3806e-23) .* critical_vec_4(i)) .* exp(- 1./((1.3806e-23) .* critical_vec_4(i)) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* (1.3806e-23).* critical_vec_4(i)) .* 2./sqrt(pi) .* int( exp(-t.^2),0, (x.*(1.3806e-23).* critical_vec_4(i))./0.4  .*(4.9-5)))==0.00000001,x)
S=double(S);
plot(critical_vec_4(i),S,'.')
hold on;
end
figure(2)
for i = 140:180
S=solve(1.6./(sqrt(pi).* x .* (1.3806e-23) .* vec_heights(i)) .* exp(- 1./((1.3806e-23) .* vec_heights(i)) .* (4.9-5)^2) - ( 0.8 ./(x.^2 .* (1.3806e-23).* vec_heights(i)) .* 2./sqrt(pi) .* int( exp(-t.^2),0, (x.*(1.3806e-23).* vec_heights(i))./0.4  .*(4.9-5)))==0.00000001,x)
S=double(S);plot(vec_heights(i),S,'.')
hold on;
xlabel('Barrier height')
ylabel('Temperature')
end

for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-A)/sigma_temp));
temp_sum = temp_sum + vpaintegral(poly_exp,0,A);
end
A=6;
height = 1.08e19
temp_sum=0;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-A)/sigma_temp));
temp_sum = temp_sum + vpaintegral(poly_exp,0,A);
end
bt_temp=600;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-A)/sigma_temp));
temp_sum = temp_sum + vpaintegral(poly_exp,0,A);
end
mu_vec=[0.999999 1.999999 2.999999 3.999999 4.999999 5.999999 6.999999];
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-A)/sigma_temp));
temp_sum = temp_sum + vpaintegral(poly_exp,0,A);
end
mu_vec(1)
mu_vec(2)
mu_vec(0)
sigma_temp = 0.4/height;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-A)/sigma_temp));
temp_sum = temp_sum + vpaintegral(poly_exp,0,A);
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-A)/sigma_temp));
temp_sum = temp_sum + ntegral(poly_exp,0,A);
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-A)/sigma_temp));
temp_sum = temp_sum + integral(poly_exp,0,A);
end
syms u
syms x
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-A)/sigma_temp));
temp_sum = temp_sum + integral(poly_exp,0,A);
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp));
temp_sum = temp_sum + integral(poly_exp,0,A);
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp));
temp_sum = temp_sum + integral(double(poly_exp),0,A);
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral(constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp)) , 0 , A);
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral(matlabFunction(str2sym(constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp)))) , 0 , A);
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral(str2sym(constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp))) , 0 , A);
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp))) , 0 , A);
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp))) , 0 , A);
disp(temp_sum)
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp))) , 0 , A);
plot(temp_sum,i)
hold on;
end
temp_sum=0;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp))) , 0 , A);
plot(temp_sum,i)
hold on;
end
temp_sum=0;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp))) , 0 , A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec= [1.08e+13 1.08e+15 1.08e+16 1.08e+17 1.08e+18 1.08e+23]
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/sigma_temp));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp))) , 0 , A);
end
sigma_vec=0.4/height_vec
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1))))) , 0 , A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1)))) , 0 , A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1))))),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1)))),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(constant_value .* (u/2) .* exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1))))),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
% poly exponential handle
temp_sum = temp_sum + integral( matlabFunction(),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
str_temp = constant_value .* (u/2) .* exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1))))
temp_sum = temp_sum + integral( matlabFunction(sym2str(str_temp)),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
str_temp = constant_value .* (u/2) .* exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1))))
temp_sum = temp_sum + integral( matlabFunction(str_temp),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
str_temp = constant_value .* (u/2) .* exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(1))
temp_sum = temp_sum + integral( matlabFunction(str_temp),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
str_temp = constant_value .* (u/2) .* exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(1)))
temp_sum = temp_sum + integral( matlabFunction(str_temp),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
str_temp = constant_value .* (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(1)))
temp_sum = temp_sum + integral( matlabFunction(str_temp),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
0.4/height_vec(1)
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1); sigma_temp=0.4/height_vec(i+1);
str_temp = constant_value .* (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(sigma_temp))
temp_sum = temp_sum + integral( matlabFunction(str_temp),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
sigma_temp=0.4/height_vec(i+1);
str_temp = constant_value .* (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp))
temp_sum = temp_sum + integral( matlabFunction(str_temp),0,A);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
sigma_temp=0.4/height_vec(i+1);
str_temp = (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/sigma_temp))
temp_sum = temp_sum + (constant_value .* integral( matlabFunction(str_temp),0,A));
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
sigma_temp=0.4/height_vec(i+1);
temp_sum = temp_sum + (constant_value .* integral( matlabFunction((u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/heigh_vec(i+1))))),0,A));
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
sigma_temp=0.4/height_vec(i+1);
temp_sum = temp_sum + (constant_value .* integral( matlabFunction((u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1))))),0,A));
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
sigma_temp=0.4/height_vec(i+1);
temp_sum = temp_sum + (constant_value .* vpaintegral( matlabFunction((u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1))))),0,A));
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
A=6;
temp_sum=0; syms u; syms x;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
sigma_temp=0.4/height_vec(i+1);
temp_sum = temp_sum + (constant_value .* vpaintegral( matlabFunction((u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1))))),0,A));
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  0 : A
exponential_handle = @(x) exp(0.4/(height_vec(i+1) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-x)/(0.4/height_vec(i+1))));
constant_value = exponential_handle(i+1);
sigma_temp=0.4/height_vec(i+1);
temp_sum = temp_sum + (constant_value .* vpaintegral( (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i+1)-A)/(0.4/height_vec(i+1)))),0,A));
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  1 : A
exponential_handle = @(x) exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/(0.4/height_vec(i))));
constant_value = exponential_handle(i);
sigma_temp=0.4/height_vec(i);
temp_sum = temp_sum + (constant_value .* vpaintegral( (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-A)/(0.4/height_vec(i)))),0,A));
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  1 : A
exponential_handle = @(x) exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/(0.4/height_vec(i))));
constant_value = exponential_handle(i);
sigma_temp=0.4/height_vec(i);
temp_sum = temp_sum + (constant_value .* double(vpaintegral( (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-A)/(0.4/height_vec(i)))),0,A)));
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  1 : A
exponential_handle = @(x) exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/(0.4/height_vec(i))));
constant_value = exponential_handle(i);
sigma_temp=0.4/height_vec(i);
x=double(vpaintegral( (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-A)/(0.4/height_vec(i)))),0,A));
disp(x)
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  1 : A
exponential_handle = @(x) exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/(0.4/height_vec(i))));
constant_value = exponential_handle(i);
x=vpaintegral( (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-A)/(0.4/height_vec(i)))),0,A);
disp(x)
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  1 : A
exponential_handle = @(x) exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/(0.4/height_vec(i))));
constant_value = exponential_handle(i);
x=integral( (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-A)/(0.4/height_vec(i)))),0,A);
disp(x)
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  1 : A
exponential_handle = @(x) exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/(0.4/height_vec(i))));
constant_value = exponential_handle(i);
x=integral( matlabFunction((u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-A)/(0.4/height_vec(i))))),0,A);
disp(x)
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  1 : A
exponential_handle = @(x) exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/(0.4/height_vec(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(1 .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-A)/(0.4/height_vec(i)))),0,A);
disp(x)
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  1 : A
exponential_handle = @(x) exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/(0.4/height_vec(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-A)/(0.4/height_vec(i)))),0,A);
disp(x)
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec)
exponential_handle = @(x) exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/(0.4/height_vec(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-A)/(0.4/height_vec(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_2= [1.0800e+18 9.5e+18 1.0800e+19 9.5e+19 1.0800e+20 9.5e+20 1.0800e+21 9.5e+21 1.080e+22 9.5e+22 1.080e+23 9.5e+23];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_2)
exponential_handle = @(x) exp(0.4/(height_vec_2(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-x)/(0.4/height_vec_2(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_2(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec(i)-A)/(0.4/height_vec_2(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
mu_vec_2 = [0.999 1.999 2.999 3.999 4.999 5.999 6.999 7.999 8.999 9.999 10.999 11.999 12.999];
mu_vec_2 = [0.999 1.999 2.999 3.999 4.999 5.999 6.999 7.999 8.999 9.999 10.999 11.999 ];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_2)
exponential_handle = @(x) exp(0.4/(height_vec_2(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_2(i)-x)/(0.4/height_vec_2(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_2(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_2(i)-A)/(0.4/height_vec_2(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_3 = [1.0800e+17 3.5e+17 9.5e+17 0.888e+18 3.5e+18 6.5e+18 1.0800e+21 3.5e+21 6.5e+21 1.5e+24 3.5e+21 1.5e+21 7.8e+20 3.1e+17 5.5e+16 1.5e+16 1e+16 1.5e+15 3.5e+14 6.5e+13]
mu_vec_3 = [0.999 1.999 2.999 3.999 4.999 5.999 6.999 7.999 8.999 9.999 10.999 11.999 12.999 13.999 14.999 15.999 16.999 17.999 18.999 19.999];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_3)
exponential_handle = @(x) exp(0.4/(height_vec_3(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_3(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_3(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_3(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_4 = [1.0800e+17 3.5e+17 9.5e+17 0.888e+18 3.5e+18 6.5e+18 1.0800e+21 3.5e+21 6.5e+21 1.5e+24 3.5e+21 1.5e+21 7.8e+20 3.1e+19 5.5e+19 1.5e+19 1e+19 1.5e+19 3.5e+19 6.5e+19];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_4)
exponential_handle = @(x) exp(0.4/(height_vec_4(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_4(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_4(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_4(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_5 = [1.0800e+17 3.5e+17 9.5e+17 0.888e+18 3.5e+18 6.5e+18 1.0800e+21 3.5e+21 6.5e+21 1.5e+24 3.5e+21 1.5e+21 7.8e+20 3.1e+19 5.5e+18 1.5e+18 1e+18 1.5e+18 3.5e+18 6.5e+18];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_5 = [1.0800e+17 3.5e+17 9.5e+17 0.888e+18 3.5e+18 6.5e+18 1.0800e+21 3.5e+21 6.5e+21 1.5e+24 3.5e+21 1.5e+21 7.8e+20 3.1e+19 5.5e+16 1.5e+16 1e+16 1.5e+16 3.5e+16 6.5e+16];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_5 = [1.0800e+17 3.5e+17 9.5e+17 0.888e+18 3.5e+18 6.5e+18 1.0800e+21 3.5e+21 6.5e+21 1.5e+24 3.5e+21 1.5e+21 7.8e+20 3.1e+19 5.5e+12 1.5e+12 1e+12 1.5e+12 3.5e+12 6.5e+12];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_5 = [1.0800e+17 3.5e+17 9.5e+17 0.888e+18 3.5e+18 6.5e+18 1.0800e+21 3.5e+21 6.5e+21 1.5e+24 3.5e+21 1.5e+21 7.8e+20 3.1e+19 5.5e+19 1.5e+19 1e+19 1.5e+19 3.5e+19 6.5e+19];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_5 = [1.0800e+17 3.5e+17 9.5e+17 0.888e+18 3.5e+18 6.5e+18 1.0800e+21 3.5e+21 6.5e+21 1.5e+24 3.5e+21 1.5e+21 7.8e+20 3.1e+19 5.5e+27 1.5e+27 1e+27 1.5e+27 3.5e+27 6.5e+27];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_5 = [1.0800e+17 3.5e+17 9.5e+17 0.888e+18 3.5e+18 6.5e+18 1.0800e+22 3.5e+22 6.5e+23 1.5e+24 3.5e+23 1.5e+23 7.8e+23 3.1e+23 5.5e+27 1.5e+27 1e+27 1.5e+27 3.5e+27 6.5e+27];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_5 = [1.0800e+17 3.5e+17 9.5e+17 0.888e+18 3.5e+18 6.5e+18 1.0800e+24 3.5e+33 6.5e+33 1.5e+34 3.5e+30 1.5e+34 7.8e+35 3.1e+35 5.5e+27 1.5e+27 1e+27 1.5e+27 3.5e+27 6.5e+27];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
height_vec_6 = [1.0800e+17 3.5e+17 9.5e+17 0.888e+18 3.5e+18 6.5e+18 1.0800e+16 3.5e+16 6.5e+16 1.5e+16 3.5e+22 1.5e+23 7.8e+23 3.1e+24 5.5e+21 1.5e+27 1e+27 1.5e+27 3.5e+27 6.5e+27];
temp_sum=0; syms u; syms x;
for i =  1 : length(height_vec_6)
exponential_handle = @(x) exp(0.4/(height_vec_6(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_6(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_6(i) .* 1.08e-23 .* bt_temp) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_6(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [600 600 600 600 600 600 600 3000 4500 4700 10000 15000 8000 4000 3000 2700 2200 1800 1700 900]; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [600 600 600 600 600 600 600 3000 4500 4700 10000 15000 8000 4000 3000 2700 2200 1800 1700 900]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [600 600 600 600 600 600 600 3000 4500 47000 100000 150000 80000 40000 3000 2700 2200 1800 1700 900]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [6e+5 6e+5 6e+5 6e+5 6e+5 6e+5 6e+3 3e+7 5e+8 6e+9 1e+15 3e+18 8e+20 4e+25 3e+28 2e+15 2e+12 1e+10 7e+7 9e+5]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [6e+5 6e+5 6e+5 6e+5 6e+5 6e+5 6e+15 3e+17 5e+18 6e+19 1e+23 3e+24 8e+35 4e+45 3e+38 2e+15 2e+12 1e+10 7e+7 9e+5]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [6e+5 6e+5 6e+5 6e+5 6e+5 6e+5 6e+15 3e+17 5e+18 6e+19 1e+23 3e+24 8e+25 4e+23 3e+20 2e+19 2e+18 1e+16 7e+7 9e+5]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; temp_vec= [6e+5 6e+5 6e+5 6e+5 6e+5 6e+5 6e+15 3e+17 5e+18 6e+45 1e+46 3e+24 8e+25 4e+23 3e+20 2e+19 2e+18 1e+16 7e+7 9e+5]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; temp_vec= [6e+5 6e+5 6e+5 6e+5 6e+5 6e+5 6e+15 3e+17 5e+18 6e+9 1e+9 3e+9 8e+9 4e+23 3e+20 2e+19 2e+18 1e+16 7e+7 9e+5]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; temp_vec= [6e-5 6e-5 6e-5 6e-5 6e-5 6e-5 6e-3 3e-3 5e-3 6e-3 1e-2 3e-2 8e-3 4e-4 3e-8 2e-7 2e-8 1e-9 7e-10 9e-11]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_sum=0; temp_vec= [6e-5 6e-5 6e-5 6e-5 6e-5 6e-5 6e-7 3e-7 5e-8 6e-10 1e-12 3e-13 8e-22 4e-14 3e-8 2e-7 2e-8 1e-9 7e-10 9e-11]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
temp_sum=0; temp_vec= [6e-5 6e-5 6e-5 6e-5 6e-5 6e-5 6e-7 3e-7 5e-8 6e-15 1e-17 3e-17 8e-18 4e-24 3e-19 2e-7 2e-8 1e-9 7e-10 9e-11]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [6e-1 6e-1 6e-1 6e-1 6e-2 6e-2 6e-5 3e-5 5e-6 6e-8 1e-7 3e-6 8e-5 4e-4 3e-3 2e-3 2e-3 1e-3 7e-2 9e-2]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [10 10 10 10 10 10 6e-5 3e-5 5e-6 6e-8 1e-7 3e-6 8e-5 4e-4 3e-3 2e-3 2e-3 1e-3 7e-2 9e-2]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [1000 1000 100000 100000 10000000 100000000000 6e-5 3e-5 5e-6 6e-8 1e-7 3e-6 8e-5 4e-4 3e-3 2e-3 2e-3 1e-3 7e-2 9e-2]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [1000 1000 100000 100000 10000000 100000000000 6e-10 3e-11 5e-12 6e-10 1e-11 3e-7 8e-9 4e-4 3e-3 2e-3 2e-3 1e-3 7e-2 9e-2]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end
temp_vec= [1000 1000 100000 100000 10000000 100000000000 6e-3 3e-4 5e-4 6e-4 1e-4 3e-7 8e-9 4e-4 3e-3 2e-3 2e-3 1e-3 7e-2 9e-2]; temp_sum=0; for i =  1 : length(height_vec_5)
exponential_handle = @(x) exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-x)/(0.4/height_vec_5(i))));
constant_value = exponential_handle(i);
x=integral(@(u) (u/2) .* exp(0.4/(height_vec_5(i) .* 1.08e-23 .* temp_vec(i)) .* erf((mu_vec_3(i)-A)/(0.4/height_vec_5(i)))),0,A);
temp_sum = temp_sum + (constant_value .* x);
disp(temp_sum)
plot(temp_sum,i)
hold on;
end


end 