function step_potential(U_initial, U_string)
format longEng
sol_array = [ ];
% solve IVP for first barrier
syms y(x)
Dy = diff(y);
ode = diff(y,x,2) == -1+(U_initial)*diff(y,x,1);
cond1 = y(0) == 0;
cond2 = Dy(0) == 1;
conds = [cond1 cond2];
ySol_1(x) = dsolve(ode,conds);
ySol = matlabFunction(ySol_1(x));
sol_array(1) = ySol(1);

for I = 1 : length(U_string)
    temp = matlabFunction(ySol_1(x));
   
    % formulate new IVP from exit time of previous solution
    ode_temp = diff(y,x,2) == -1+(U_string(I))*diff(y,x,1);
    cond_1 = y(0)==temp(1);
    cond_2 = Dy(0)==1;
    conds_temp = [cond_1 cond_2];
    ySol = dsolve(ode_temp,conds_temp);
    temp_temp = matlabFunction(ySol);
    disp(temp_temp(1))
    sol_array(I+1) = temp_temp(1);
end

% construct figures for different barrier heights

time_vec = zeros(1,1000);

test_vec = linspace(1,500,500);
time_vec(1) = 0;

J=50;
matrix_test = zeros([J 500]);

% generate test cases for potential
% for J = 1 : 4
    % vec_temp = zeros([1 500]);
    % for index = 1 : 500
        % if index <= double(500/(2*J))
            % vec_temp(index) = 100000 - (0.001 * index * log(index)) - (1.0000005 * index^2) + (30 .* log(index^3) .* index^3) + (0.09 .* log(index^4) .* index^4) - (0.007 .* log(index^5) .* index^5) + (50 .* log(index^6) .* index^6);
        % else 
            % vec_temp(index)=5;
        % end 
    % end
    % matrix_test(J,:) = vec_temp;
    
% end 

for J = 1 : 4
    for AJ = 1 : 250
        for height = 1000 : 10000
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
    end 
end 

    J=2;
    for AJ = 251 : 500
        for height = 8800 : 1300
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
    end
    
    for AJ = 501 : 900
        for height = 100 : 1100
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
    end
    
    J=3;
     for AJ = 251 : 600
        for height = 15000 : 20000
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
     end
     
     disp('third row')
     
     for AJ = 601 : 900
        for height = 1500 : 2000
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
     end
     

    J=4;
    for AJ = 251 : 800
        for height = 21000 : 21500
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
    end
     
    for AJ = 801 : 900
        for height = 2100 : 2150
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
     end

for next_temp = 1 : 46
    x = matrix_test(1,:);
    x = x(randperm(length(x)));
    matrix_test(4+next_temp,:)=x;
end 

for B = 1: 50
    figure(B)
    for A = 1 : 900
    % construct gaussians 
    mu = 0.5;
    disp('mu')
    xmin = 500;
    xmax = 3; 
    
    syms y(x);
    % initialize different potential from product of all distributions
    sigma_temp = matrix_test(B,A);
    func_exp_1 = @(x, mu, sigma) exp(-((x-(mu))./sigma).^2);
    func_next = func_exp_1(x , mu + A-1, sigma_temp);
    
    % construct product of variances
    func_variance = @(sigma) sqrt(2 .* pi .* sigma);
    temp_next = func_variance(matrix_test(B,A));
    func_boltzmann_temp = @(temp) (1.380649e-23).* temp;
    bt_temp = func_boltzmann_temp(1000); 
    disp(bt_temp)
    
    % store result
    syms y(x);
    ode = diff(y,x,2) == - 1 + (func_next) .* diff(y,x,1)./(temp_next.* bt_temp);

    cond_1 = y(0)==0;
    cond_2 = Dy(0)==1;
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
            % iterate through output displayed above
            x_pos_temp_1 = strfind(temp_str_zero, '(');
            x_pos_temp_2 = strfind(temp_str_zero, '*');
            x_pos_temp_3 = strfind(temp_str_zero, '/');
            x_pos_temp_4 = strfind(temp_str_zero, ')');
            
            % first part of exponent 
            temp_str_zero_1 = temp_str_zero(x_pos_temp_1(2)+1:x_pos_temp_2(1)-1);
            disp(temp_str_zero_1)
            % second part of exponent
            func_part_two = @(x) sqrt( x .* pi);
            func_part_three = @(x) exp(x);
            func_part_three_temp = 1./func_part_three(1./3);
            func_part_two_temp = func_part_two(1);
            final_result_funcs_two_three = func_part_three_temp .* func_part_two_temp;
            disp('two and three crossterm result')
            disp(final_result_funcs_two_three)
            
            % error function part of the exponent
            
            x_pos_erf = strfind(temp_str_zero, 'f');
            erf_argument = temp_str_zero(x_pos_erf(1)+1:x_pos_temp_3(length(x_pos_temp_3))-2);
            disp('erf-argument')
            disp(erf_argument)
            erf_argument = str2sym(erf_argument);
            disp(erf_argument)
            term_erf = erf(erf_argument);
            term_erf = double(term_erf);
            disp('numerical erf value')
            disp(term_erf)
            % add in conditions to reflect extra '/' character in erf
            
            % if isequal(length(x_pos_temp_3),3)
                % temp_str_zero_3 = temp_str_zero(x_pos_temp_3(3)+1: x_pos_temp_4(length(x_pos_temp_4))-1);
            % else
                % temp_str_zero_3 = temp_str_zero(x_pos_temp_3(2)+1: x_pos_temp_4(length(x_pos_temp_4))-1);   
            % end 
            
            disp('test')
            if isequal(length(x_pos_temp_3),5)
            temp_str_zero_3 = temp_str_zero(x_pos_temp_3(5)+1 : x_pos_temp_4(length(x_pos_temp_4))-1);
            temp_str_zero_3 = convertCharsToStrings(temp_str_zero_3);
            temp_str_zero_3 = double(temp_str_zero_3);
            disp(temp_str_zero_3)
            elseif ~ isequal(length(x_pos_temp_3),5)
            temp_str_zero_3 = temp_str_zero(x_pos_temp_3(length(x_pos_temp_3))+1 : x_pos_temp_4(length(x_pos_temp_4))-1);
            temp_str_zero_3 = convertCharsToStrings(temp_str_zero_3);
            temp_str_zero_3 = double(temp_str_zero_3);
            disp(temp_str_zero_3)
            end 
            % convert each temp variable above to function handles
            disp('temp_str_zero_1')
            disp(temp_str_zero_1)
            temp_str_zero_1 = convertCharsToStrings(temp_str_zero_1);
            disp('temp_str_zero_1')
            disp(temp_str_zero_1)
            temp_str_zero_1 = double(temp_str_zero_1);
            disp('temp_str_zero_1')
            disp(temp_str_zero_1)
            exponent_power = - temp_str_zero_1 .* final_result_funcs_two_three .* term_erf ./ temp_str_zero_3;
            disp('final exponent power result')
            disp(exponent_power)
            temp_str_zero_3 = convertCharsToStrings(temp_str_zero_3);
            
            exponent = exp(exponent_power);
            disp('exponent term')
            disp(exponent)
            
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
            disp(class(temp_str_one))
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
            func_temp = str2sym(temp_str_one);
            integral_term_one = vpaintegral(func_temp, A , A+1);
            disp('printing first term')
            disp(integral_term_one)
            disp('COMPLETE FIRST TERM')
            result = exponent * integral_term_one;
            disp(result)
            
            % repeat
            % temp_str_two = temp_sol(temp_pos(temp_index+1):temp_pos(temp_index+2)-3);
            % disp(temp_str_two)
            % array_temp(temp_index+2) = str2num(temp_pos(temp_index+1)+1:temp_pos(temp_index+2)-3);
            % x_1 = strfind(temp_str_two,'(');
            % x_2 = strfind(temp_str_two,',');
            % temp_str_two = temp_str_two(x_1(1)+1: x_2(1)-1);
            % temp_str_two=convertCharsToStrings(temp_str_two);
            % disp(class(temp_str_two))
            % func_temp = matlabFunction(str2sym(temp_str_two));
            
            % disp('second integral term')
            syms u;
            syms x;
            % integral_term_two = integral(func_temp, 0, A, 'RelTol',0);
            % integral_term_two = (integral_term_two)^2;
            % disp(integral_term_two)
            % disp('FIRST AND SECOND TERMS')
            % result = result - integral_term_two;
            % disp(result)
            % array_temp(temp_index+2) = integral_term_two;

            % repeat, inner integral for fourth & fifth integral terms
            temp_str_four = temp_sol(temp_pos(temp_index+3)+1:temp_pos(temp_index+4)-3);
            % array_temp(temp_index+4) = str2num(temp_pos(temp_index+4)+1:temp_pos(temp_index+5)-3);

            temp_str_four=convertCharsToStrings(temp_str_four);
            % temp_str_four = matlabFunction(str2sym(temp_str_four));
                        
            disp('double integral term')
            syms u;
            % construct trap approx
            func_temp = matlabFunction(func_temp); 
            trap_func_1 = sym(func_temp(0));
            disp(trap_func_1)
            trap_func_2 = sym(func_temp(u./2));
            disp(trap_func_2)
            trap_func_3 = sym(func_temp(u));
            disp(trap_func_3)
            trap_final = (u./2) .* symsum(trap_func_1 , trap_func_2,trap_func_3);
            trap_final = trap_final .* str2sym(temp_str_one);
            final_handle = convertCharsToStrings(trap_final);
            disp(trap_final)
            integral_term_four_five = vpaintegral(final_handle,A, A+1);

            disp(integral_term_four_five)
            % array_temp(temp_index+4) = integral_term_four_five;
            
            % add 3 
            disp('FIRST, SECOND AND THIRD TERMS')
            final= 1 + integral_term_four_five;
            disp(final)
            disp('end')

            subplot(5,1,1)
            title('Relationship between exit time and base pair of the target sequence')
            xlabel('Base pair')
            ylabel('Exit time')
            plot(A,double(final),'.')
            hold on;
            
            % plot inverse of exit times
            subplot(5,1,2)
            title('Plotting the inverse of the exit time')
            xlabel('Base pair')
            ylabel('Exit time')
            plot(A,1./double(final),'.')
            hold on;
            
            subplot(5,1,3)
            title('Relationship between barrier height and magnitude of exit time')
            xlabel('Base pair')
            ylabel('Fluctuation height')
            plot(A,matrix_test(B,A),'.')
            hold on;
            
            subplot(5,1,4)
            title('Relationship between exit time and base pair of the target sequence')
            xlabel('Base pair')
            ylabel('Exit time')
            plot(matrix_test(B,A),A,'.')
            hold on;
            
            subplot(5,1,5)
            title('Three dimensional solution curve')
            xlabel('Base pair')
            ylabel('Exit time')
            zlabel('Temperature')
            plot3(A,double(final),bt_temp,'.')
            hold on;
            
            bt_temp = bt_temp+1.5;
    end   
end  
end  
    
    % figure(temp);
    % plot(Cell, linspace(1,10,10))
    % hold on;
    
    % hold off;
 
 


