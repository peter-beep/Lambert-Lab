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

for temp = 1: 100: 501
    figure(temp)
    for  A = 1 : length(test_vec) 

    % construct gaussians 
    mu = 0.5;
    disp('mu')
    sigma = 0.4./(A.*temp);
    temp_mu = mu * A;
    % temp_mu = sym(temp_mu);
    syms y(x);
    func_exp_1 = @(x, mu) exp(-(x-(mu)).^2);
    % store result
    syms y(x);
    ode = diff(y,x,2) == - 1 - (func_exp_1(x,temp_mu)) .* diff(y,x,1)./sqrt(2 .* pi .* sigma);
    if isequal(A,1)
    cond_1 = y(0)==0;
    cond_2 = Dy(0)==1;
    conds_temp = [cond_1 cond_2];
    ySol = dsolve(ode,conds_temp);
    else 
    cond_1 = y(0)==1;
    cond_2 = Dy(0)==(1);
    conds_temp = [cond_1 cond_2];
    ySol = dsolve(ode,conds_temp);    
    end
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
            % check for length of x_pos_temp_3
            
            % first part of exponent 
            temp_str_zero_1 = temp_str_zero(x_pos_temp_1(2)+1:x_pos_temp_2(1)-1);
            disp(temp_str_zero_1)
            % second part of exponent
            temp_str_zero_2 = erf(A./2);
            disp(temp_str_zero_2)
            % third part of exponent
            % add in conditions to reflect extra '/' character in erf
            if isequal(length(x_pos_temp_3),3)
                temp_str_zero_3 = temp_str_zero(x_pos_temp_3(3)+1: x_pos_temp_4(length(x_pos_temp_4))-1);
            else
                temp_str_zero_3 = temp_str_zero(x_pos_temp_3(2)+1: x_pos_temp_4(length(x_pos_temp_4))-1);   
            end 
            
            disp('test')
            disp(temp_str_zero_3)
            % convert each temp variable above to function handles
            temp_str_zero_1 = convertCharsToStrings(temp_str_zero_1);
            temp_str_zero_3 = convertCharsToStrings(temp_str_zero_3);
            
            func_temp_1 = str2double(temp_str_zero_1);
            func_temp_3 = str2double(temp_str_zero_3);
            
            % exp power
            power_temp = (func_temp_1/func_temp_3) * temp_str_zero_2 * - sqrt(pi);
            
            % exponentiate product of three terms
            first_term = exp(power_temp);
            disp('first term, error function integral')
            disp(first_term)
            
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
            func_temp = matlabFunction(str2sym(temp_str_one));
            integral_term_one = integral(func_temp, 0 , 1, 'RelTol',0);
            disp('printing first term')
            disp(integral_term_one)
            disp('COMPLETE FIRST TERM')
            result = first_term * integral_term_one;
            disp(result)
            
            % repeat
            temp_str_two = temp_sol(temp_pos(temp_index+1):temp_pos(temp_index+2)-3);
            disp(temp_str_two)
            % array_temp(temp_index+2) = str2num(temp_pos(temp_index+1)+1:temp_pos(temp_index+2)-3);
            x_1 = strfind(temp_str_two,'(');
            x_2 = strfind(temp_str_two,',');
            temp_str_two = temp_str_two(x_1(1)+1: x_2(1)-1);
            temp_str_two=convertCharsToStrings(temp_str_two);
            disp(class(temp_str_two))
            func_temp = matlabFunction(str2sym(temp_str_two));
            
            disp('second integral term')
            syms u;
            syms x;
            integral_term_two = integral(func_temp, 0, 1, 'RelTol',0);
            integral_term_two = (integral_term_two)^2;
            disp(integral_term_two)
            disp('FIRST AND SECOND TERMS')
            result = result - integral_term_two;
            disp(result)
            % array_temp(temp_index+2) = integral_term_two;

            % repeat, inner integral for fourth & fifth integral terms
            temp_str_four = temp_sol(temp_pos(temp_index+3)+1:temp_pos(temp_index+4)-3);
            % array_temp(temp_index+4) = str2num(temp_pos(temp_index+4)+1:temp_pos(temp_index+5)-3);

            temp_str_four=convertCharsToStrings(temp_str_four);
            % temp_str_four = matlabFunction(str2sym(temp_str_four));
                        
            disp('double integral term')
            syms u;
            % construct trap approx
            syms u;
            trap_func_1 = sym(func_temp(0));
            disp(trap_func_1)
            trap_func_2 = sym(func_temp(u./2));
            disp(trap_func_2)
            trap_func_3 = sym(func_temp(u));
            disp(trap_func_3)
            trap_final = (u./2) .* symsum(trap_func_1 , trap_func_2,trap_func_3);
            trap_final = trap_final .* str2sym(temp_str_one);
            trap_final = matlabFunction(trap_final);
            final_handle = trap_final;
            disp(trap_final)
            integral_term_four_five = integral(final_handle,0, 1, 'RelTol',0);

            disp(integral_term_four_five)
            % array_temp(temp_index+4) = integral_term_four_five;
            
            % add 3 
            disp('FIRST, SECOND AND THIRD TERMS')
            final= 3 + integral_term_four_five;
            disp(final)
            disp('end')
            
            time_vec(A+1) = double(final);
            disp(time_vec(A+1))
            disp(A)
            subplot(2,1,1)
            plot(A,double(final),'.')
            hold on;
            
            % plot inverse of exit times
            subplot(2,1,2)
            plot(A,1./double(final),'.')
            hold on;
    end 
    
    % figure(temp);
    % plot(Cell, linspace(1,10,10))
    % hold on;
    
    hold off;

end   

end 


