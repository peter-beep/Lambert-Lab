function et_fluctuation(n, mu_1, mu_2, sigma_1,T,delta)
for sigma_2=0.0000000000000005:0.000000000000000000001:(0.0000000000000005+(0.000000000000000000001*1000))
height_1= 0.4./sigma_1;
height_2= 0.4./sigma_2;
syms x 
x_1 = @(x) (mu_2-x)./sigma_2;
x_2 = @(x) (mu_1-x)./sigma_1;
temp_fac = 0.4./T;
fluctuation_exp_1 = @(x) exp( temp_fac .* ((0.5.* erf(x_1(x)))-erf(x_2(x))));
syms u 
val_exp_1_0 = sym(fluctuation_exp_1(0));
val_exp_1_1 = sym(fluctuation_exp_1(x./2));
val_exp_2 = sym(fluctuation_exp_1(x));

trap_exp_1 = (x./2) .* symsum(val_exp_1_0 , val_exp_1_1, val_exp_2) .* (exp((0.4./(height_2.* T)).* erf(x_1(x))));
trap_exp_2 = (x./2) .* symsum(val_exp_1_0 , val_exp_1_1 , val_exp_2) .* exp((0.4./(height_1*T)).*erf(x_2(x)));
final_handle_1 = convertCharsToStrings(trap_exp_1);
final_handle_1 =  matlabFunction(final_handle_1);
final_handle_2 = convertCharsToStrings(trap_exp_2);
final_handle_2 = matlabFunction(final_handle_2);
approx_val_1= integral(final_handle_1,0,n);
approx_val_2 = integral(final_handle_2,0,n-1);
final_res_1 = approx_val_1-approx_val_2;

syms x 
syms h
second_handle_1 = @(x,h) (exp(temp_fac .*((0.5 .* (erf(x_1(x))- erf(x_1(0))) -(erf(x_2(x))- erf(x_2(0))))))-1).*(h-2) .* exp((0.4./(height_2.*T)).* erf(x_1(x)));
result_second_handle_1 = convertCharsToStrings(second_handle_1(x, height_2));
result_2_1= vpaintegral(result_second_handle_1,0,(n-1)+ delta);

syms x
second_handle_2 = @(x,h) (exp(temp_fac .* ((0.5 .* (erf(x_1(x)) - erf(x_1(0))) - (erf(x_2(x))-erf(x_2(0)))))) -1) .*(h-2) .* exp((0.4./(height_1.* T)).* erf(x_2(x)));
result_second_handle_2 = convertCharsToStrings(second_handle_2(x, height_1));
result_2_2 = vpaintegral(result_second_handle_2,0,n-1);

final_res_2 = result_2_1 + result_2_2;

disp(final_res_1)
disp(result_2_1)
disp(result_2_2)
exit_time_fluctuation= double((2.*n)+1 + final_res_1 + final_res_2);
disp(exit_time_fluctuation)
figure(3)
plot([(n-1)+delta exit_time_fluctuation ],[n-1 1])
hold on;
end
end