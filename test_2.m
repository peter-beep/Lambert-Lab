figure(7)

syms t
func_test_int = @(x) integral(matlabFunction(exp(-t.^2)) ,0,x);

for A=1:50
    disp(func_test_int(A))
    plot(A,func_test_int(A));
    hold on;
    
end 