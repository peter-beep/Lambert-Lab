function [Table, T_binary,t_matrix_1, t_matrix_2, temp_Ising] = rw(N)
% Read in data from .csv files for sequences
% and for binding energies
Table = readtable('First6-Data.csv');
% Table = table2cell(Table);
% Table = table2array(Table);
% x = size(Table{:,1});
% for A = 1 : x(1)
%     x = 
% end 

T_binary = readtable('BinaryMismatchData.csv');
% T_binary = T_binary
Table_4 = readtable('Last14Data.csv');
Length_4 = size(Table_4(:,1));

binding_position=N;
beta = 0.5;

% store positions of sequences for which initial_t_prob 
% gives a negative value from either the disagreeing
% or agreeing base case

test_x = size(T_binary(:,1));
bad_positions = zeros(test_x(1),20);

% initialize t prob to a constant value
initial_t_prob = exp(-3/2);

% or, initialize t prob to a value based on agreement 
% from the PAM sequences 
Table_3 = readtable('PAMData.csv');
PAM_length = size(Table_3(:,1));
Table_3_PAM_sequences = zeros(PAM_length(1),6);

% for A = 1:6
% for JJ = 1 : PAM_length(1)

%     TEST_3 = Table_3{JJ,1};
%     TEST_3 = cell2mat(TEST_3);
%     
%     % Table_3_PAM_sequences(JJ,:) = TEST_3;
% end 
% end 

% initialize empty matrix for Ising spins
Test_x = size(Table(:,1));
master_matrix = zeros(Test_x(1),22);

% transition probability matrix, Hamiltonian 1
TTest_x = size(Table(:,1));
t_matrix_1 = zeros(TTest_x(1),20);

% transition probability matrix, Hamiltonian 2
TTTest_x = size(Table(:,1));
t_matrix_2 = zeros(TTTest_x(1), 20);
t_matrix_3 = zeros(TTTest_x(1), 20);
t_matrix_4 = zeros(TTTest_x(1), 20);

test_sequences = { };
PAM_sequences = { };
fold_changes = { };
occupancy = { };
PAM_occupancy = { };

% assigin Ising spins, first through empty vector
% outside of the for loop
    temp_Ising = zeros(TTTest_x(1),20);

for I = 1 : length(size(Table(:,1)))
    TEST = Table{I,1};
    TEST = cell2mat(TEST);
    % sequence_tDNA = x(22:40);
    % obtain individual base pairs
    for J = 1 : 20
        temp1 = TEST(J:J);
        temp2 = TEST(J+20:J+20);
        
        if (eq(temp1, 'A'))
            if (eq(temp2, 'T'))
                temp_Ising(I,J)=1;
                
            elseif ~(eq(temp2, 'T'))
                temp_Ising(I,J)=0;
                
            end 
        end
        
        if (eq(temp1, 'T'))
            if (eq(temp2, 'A'))
                temp_Ising(I,J)=1;
            elseif ~(eq(temp2, 'A'))
                temp_Ising(I,J)=0;
                
            end 
                
        end  
     
        if (eq(temp1, 'G'))
            if (eq(temp2, 'C'))
                temp_Ising(I,J)=1;
            elseif ~(eq(temp2, 'C'))
                temp_Ising(I,J)=0;
                
            end 
        end 
        
        if (eq(temp1, 'C'))
            if (isequal(temp2, 'G'))
                temp_Ising(I,J)=1;
             elseif ~(isequal(temp2, 'G'))
                temp_Ising(I,J)=0;
            end 
        end      
    end 
    % assign temp_Ising vec to master matrix
    
    % master_matrix(I) = temp_Ising;
end  
    
% for A = 1 : 20
%     x=0;
%     temp_pos = [ ];
%    while (isequal(master_matrix(A),1) && isequal(master_matrix(A+1),1))
%       temp_pos(A) = x+1;
%    end
%    
%    master_matrix(A,21) = min(temp_pos);
%    master_matrix(A,22) = max(temp_pos);
% end 

% Apply transition probability formulas depending on
% the values of Ising spins

t_matrix_1(:,1) = initial_t_prob;
t_matrix_2(:,2)= initial_t_prob;

% for now, look at just the 4000 sequences
% from the binary mismatch file
FINAL = size(Table(:,1));

for IJ = 1 : PAM_length(1)
    for I = 1 : FINAL(1)
    for JI = 1 : Length_4
    % iteratively define x for each I
    % supress information for I index in each loop
    
%     TEST = Table{I,1};
%     TEST = cell2mat(TEST);
    
     % x= temp_Ising(I,:);
     
     % generate random PAM sequences
     x = randseq(6);
    
    % define parameters lambda of mismatch
    for J = 2 : 20
       
        % define a parameter lambda for the exponential random variable
        
        % for a mismatch in the first 6 positions, we 
        % assign the highest exponential penalty
        if ((J>=2) && (J<=6))
            lambda_mismatch = exp(-5/2);
            lambda_c = 0.30;
            lambda_p = 0.25;
        end 
        
        % for the middle 6 positions, we assign a
        % smaller exponential penalty
        if ((J>=7) && (J <=13))
            lambda_mismatch = exp(-8/2);
            lambda_c = 0.25;
            lambda_p = 0;
        end 
        
        % while for the remaining 6 positions, we assign
        % the smallest exponential penalty for mismatching
        % bases
        if ((J>=14) && (J<=20))
            lambda_mismatch = exp(-10/2);
            lambda_c = 0.10;
            lambda_p = 0;
        end 

        % define couplings, case 1
        if eq(temp_Ising(I,J),1)
            couplings_ij1 = exp(N-J);
        elseif eq(temp_Ising(I,J),0)
            couplings_ij1 = 1-exp(N-J);
        end 
        
        % define couplins, case 2
         if eq(temp_Ising(I,J),1)
            couplings_ij2 = N-J;
        elseif eq(temp_Ising(I,J),0)
            couplings_ij2 = 1-(N-J);
        end 
        
        % Hamiltonian 1
        hamiltonian_i1 = exp(couplings_ij1 * temp_Ising(I,J)) * exp(lambda_mismatch * (N-J));
        % Hamiltonian 2
        hamiltonian_i2 = exp(couplings_ij2 * temp_Ising(I,J)) * exp(lambda_mismatch * (N-J));
        
        % project each row of temp_Ising onto one less dimension, for easier
        % indexing in the for loops below
        x= temp_Ising(I,:);
        
        % compute transition probability formula for agreeing
        % and disagreeing bases
        if eq(x(J),1)
            % for pos = J : N
            func = @(X) (hamiltonian_i1.* exp(X))./(1 + (lambda_p .* exp(- beta * Table_3{IJ,8})) + (lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));
            new_t_prob = initial_t_prob - integral(func,J,N,'ArrayValued', true);
            
            func_truncated = @(X) (hamiltonian_i1.* exp(X))./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));
            new_t_prob_truncated = initial_t_prob - integral(func_truncated,J,N,'ArrayValued', true);
            % update matrix entry fo each iteration of while loop
            if (new_t_prob <0 || new_t_prob_truncated<0)
                disp('error')
                % bad_positions(I,J) = x(J)
            end 
            t_matrix_1(I,J) = new_t_prob;
            t_matrix_3(I,J) = new_t_prob_truncated;
            % end 
        elseif eq(x(J),0)
            % alternative transition probability formula
            new_t_prob = initial_t_prob - ((hamiltonian_i1)/((lambda_c .* exp(-beta * Table{I,8}))  ));
            % update matrix entry for each iteration of while loop
            t_matrix_1(I,J) = new_t_prob;
        end 

        % compute transition probabilities for second Hamiltonian 
        if eq(x(J),1)
            % for pos = J : N
            fun = @(X) (hamiltonian_i2.* exp(X))/(1 + (lambda_p .* exp(-beta * Table_3{IJ,8}))+(lambda_c .*exp(- beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(-lambda_mismatch .* (N-X))));
            new_t_prob = initial_t_prob - integral(fun,J,N,'ArrayValued', true);
            
            fun_truncated = @(X) (hamiltonian_i2.* exp(X))/(1 + (lambda_c .*exp(- beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(-lambda_mismatch .* (N-X))));
            new_t_prob_truncated = initial_t_prob - integral(fun_truncated,J,N,'ArrayValued', true);
            % update matrix entry fo each iteration of while loop
            if (new_t_prob <0 || new_t_prob_truncated <0)
                disp('error')
                % bad_positions(I,J) = x(J)
            end 
            t_matrix_2(I,J) = new_t_prob;
            t_matrix_4(I,J) = new_t_prob_truncated;
            % end 
        elseif eq(x(J),0)
            % alternative transition probability formula
            new_t_prob = initial_t_prob - ((hamiltonian_i2)/str2double(Table{I,8}));
            t_matrix_2(I,J) = new_t_prob;
            t_matrix_4(I,J) = new_t_prob_truncated;
        end  
    end 
   test_sequences{I} = x;
   PAM_sequences{I} = Table_3{IJ};
   % after saving the generated sequence, also store FC and occupancies
   XY = 1/ (1 + (lambda_c *exp(- beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})));
   fold_changes{I} = XY;
   occupancy{I} = 1 - XY;
   XYX = (lambda_p * exp(-beta * Table_3{IJ,8}))/(1+ (lambda_c .*exp(- beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})));
   PAM_occupancy{I} = XYX;
    end   
    end
    test_sequence{I+1} = 'END OF ONE ITERATION';
end

func_array_values = { };

for I = 1 : FINAL(1)
    for IJ = 1 : Length_4(1)
        func_FC = @(lambda_c) (1 + (lambda_c * exp(- beta .* Table{I,8})) + (lambda_c * exp(-beta * Table_4{IJ,8})))^{-1};
        func_array_values{I,JI} = func_FC(X);
    end 
end 

% From store function values, calculate the greatest change in values of
% the Hamiltonian

length = size(func_array_values);
finite_diff_func_array_values = { };

for I = 1 : FINAL(1)
for II = 1 : length(1)
    finite_diff_func_array_values{II} = func_array_values{I,II+1} - func_array_values{I,II};
end 
end

max_difference = max(finite_diff_func_array_values);

% plot different FC values, by inserting linear segments of growth, and
% then smooth neighborhoods of the corners at which such segments are 
% inserted, by changing the values of func_FC in some open neighborhood
% surrounding the point.

for I = 1 : FINAL(1)
    for IJ = 1 : Length_4(1)
    func_FC = @(lambda_c) (1 + (lambda_c * exp(- beta .* Table{I,8})) + (lambda_c * exp(-beta * Table_4{IJ,8})))^{-1};
    fplot(@(X) func_FC,[0 6],'b')
    hold on
    fplot(@(X) func_FC,[6 20],'r')
    hold off
    grid on
    end 
end 

% plot occupancy times for the Cas protein, from the test_sequences array
h = 0.5;
length_2 = size(test_sequences);
for JK = 1 : length_2(1)
    if test_sequences(JK) > h
        
    end 
end 

% Put all of the plots together, and update for each position of the target
% sequence

% figure(1)
%     subplot(2,1,1)
%     plot(vector_2,vector_1, 'r+' , vx , vy , 'b-')
%     xlim([0 200])
%     ylim([0 200])
%     
%     subplot(2,1,2)

end 
