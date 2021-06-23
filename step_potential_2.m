function step_potential_2()
format longEng


mu_2 = [0.0000000004 0.00000000001 0.000000000000004 0.00000000000003 0.000000000000000000000001 0.000000000000000000000000000001 0.0000000000000000000000000000011111 0.0000000000000000000000000000000000000000000000120 0.0000000000000000000000000000000000000000000000020 0.000000000070 0.0090 0.0040 10 15 20 30 20 1 0.1 0.1];
temp_exit=0;

    for A = 1 : length(mu_2)
    % construct gaussians 
    mu = 0.99999999999+ (A-1);
    disp('mu')

    syms y(x);
    % initialize different potential from product of all distributions
    height = mu_2(A);
    sigma_temp = 0.4./height;
    func_exp_1 = @(x, mu, sigma) exp(-((x-(mu))./sigma).^2);
    func_next = func_exp_1(x , mu + A-1, sigma_temp);
    
    % construct product of variances
    func_variance = @(sigma) sqrt(2 .* pi .* sigma);
    temp_next = func_variance(sigma_temp);
    func_boltzmann_temp = @(temp) temp;
    bt_temp = func_boltzmann_temp(600); 
    disp(bt_temp)
    
    % store result
    syms y(x);
    ode = diff(y,x,2) == - 1./(bt_temp) + (func_next) .* diff(y,x,1)./(bt_temp);
    Dy = diff(y);
    cond_1 = y(A)==0;
    cond_2 = Dy(0)==0;
    conds_temp = [cond_1 cond_2];
    ySol = dsolve(ode,conds_temp);

    pretty(ySol)
    y = sym(ySol);
    y = char(y);
    temp_sol=y;
    disp('sol')
    % process string contents of solution
    temp_pos = strfind(y,'int');
    ySol = dsolve(ode,conds_temp);
    
    pretty(ySol)
    y = sym(ySol);
    % temp_sol = y;
    % array_temp = zeros(1,7 * length(test_vec));
   
        temp_index=1;
            % array_temp(temp_index) = 
            disp('begin')
            disp('exponential term')
            temp_str_zero = temp_sol(1: temp_pos(temp_index)-2);
            disp('fist piece')
            disp(temp_str_zero)
            
            % first part of exponent 
            

            temp_str_one = temp_sol(temp_pos(temp_index) :temp_pos(temp_index+1)-3);
            disp('first integral term')
            disp(temp_str_one)
            x_1 = strfind(temp_str_one, '(');
            x_2 = strfind(temp_str_one,',');
            % temp_str_one = temp_sol(temp_pos(temp_index)+1 :temp_pos(temp_index+1)-3);
            % array_temp(temp_index+1) = 
            % temp_str_one = str2num(temp_sol(temp_pos(temp_index) :temp_pos(temp_index+1)-3));
            temp_str_one = temp_sol(temp_pos(temp_index) :temp_pos(temp_index+1)-3);
            % disp(temp_str_one)
            % loop over temp_str to numerically integrate part of solution
            x_1 = strfind(temp_str_one,'(');
            x_2 = strfind(temp_str_one,',');
            temp_str_one = temp_str_one(x_1(1)+1: x_2(1)-1);
            disp('string output')
            % disp(class(temp_str_one))
            disp(temp_str_one)
            disp(length(temp_str_one))
            % place all contents of string in first entry
            temp_str_one=convertCharsToStrings(temp_str_one);
            disp(class(temp_str_one))
            % disp(temp_str_one)
            % temp_str_one=str2double(temp_str_one);
            syms u;
            % temp_str_one = append("@(u)", temp_str_one);
            disp(temp_str_one)
            % temp_str_one = str2sym(char(temp_str_one));
            syms u;
            syms x;
            func_temp = matlabFunction(str2sym(temp_str_one));

            disp('FIRST, SECOND AND THIRD TERMS')
            
            temp_exp_sum=0;
           % define function handle array 
           syms x
           for i =  0 : A
             exponential_handle = @(x) exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-x)/sigma_temp));
             constant_value = exponential_handle(i+1);
             % poly exponential handle
             poly_exp = @(x) constant_value .* (u/2) .* exp(0.4/(height .* 1.08e-23 .* bt_temp) .* erf((mu-A)/sigma_temp));
             temp_sum = temp_sum + vpaintegral(poly_exp,0,A);
           end 
           

            final= @(x) -(A^2./bt_temp - x^2./bt_temp - integral_term_four_five);
            disp(final(A-1))
            tmpr = final(A-1);
            temp_exit = tmpr+temp_exit;
            
            disp('end')

            % exp((1329227995784915872903807060280344576*pi^(1/2)*erf(25000*u - 12500))/14546975136154065625) .* exp((1329227995784915872903807060280344576*pi^(1/2)*erf(25000*v - 12500))/14546975136154065625),  1, 2,0,1)


            figure(2)
            subplot(2,1,1)
            title('Relationship between exit time and base pair of the target sequence')
            xlabel('Base pair')
            ylabel('Exit time')
            plot(A,double(temp_exit),'.')
            hold on;
            
            % plot inverse of exit times
            subplot(2,1,2)
            title('Inverse of the exit time')
            xlabel('Base pair')
            ylabel('Reaction time')
            plot(A,-1./double(temp_exit),'.')
            hold on;

            
            bt_temp = bt_temp+10;
    end   

end  
    
    % figure(temp);
    % plot(Cell, linspace(1,10,10))
    % hold on;
    
    % hold off;
 
 


