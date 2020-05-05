function demo()

% READ IN TABLE DATA
Table = readtable('First6-Data.csv');
Table_3 = readtable('PAMData.csv');
Table_4 = readtable('Last14Data.csv');

% define positions for terms in partition function
I = 30;
JI = 45;
IJ=20;
beta = 0.005;
N = 16;
N_1 = 13;

% define possible cases of 3 mismatches
number_1 = 1;
number_2 = 2;
number_3 = 3;

J=3;
lambda_mismatch = exp(-5/2);
lambda_c = 0.025;
lambda_p = 0.015;
couplings_ij1 = exp(N-J);
couplings_ij2 = N-J;
hamiltonian_i1 = exp(couplings_ij1 * 1 * exp(lambda_mismatch * (N-J)));
hamiltonian_i2 = exp(couplings_ij2 * 1) * exp(lambda_mismatch * (N-J));

% generate test functions, given specific instances of binding energies
func = @(X) (hamiltonian_i1.* exp(X))./(1 + (lambda_p .* exp(- beta * Table_3{IJ,8})) + (lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));
func_truncated = @(X) (hamiltonian_i1.* exp(X))./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));
func_2 = @(X) (hamiltonian_i2.* exp(X))./(1 + (lambda_p .* exp(- beta * Table_3{IJ,8})) + (lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));
func_truncated_2 = @(X) (hamiltonian_i2.* exp(X))./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));
func_truncated_ii1 = @(X) exp(hamiltonian_ii)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate FIRST CASE for base pair mismatch, with position of 
% binding, N = 16, with a SINGLE base pair mismatch past the 
% first 6 positions

normalization_1_bm = 1/(((N-number_1-1+2))^3);
normalization_1_bmis = 1/((number_1+2)^3);

% construct the Hamiltonian. 0.5 comes from the
% T_ij factor for the single mismatch, (N-8)/N = 0.5 (!!!)

Hamiltonian_C1 = exp(-((normalization_1_bm * (((16-1)+(16-2)+(16-3)+(16-4)+(16-5)+(16-6)+(16-7)+(16-9)+(16-10)+(16-11)+(16-12)+(16-13)+(16-14)+(16-15))))+(normalization_1_bmis *(((16-8)-1)*(0.5)))));

% generate a test function with specific normalizations

func_C1 = @(X) (Hamiltonian_C1)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate SECOND CASE for base pair mismatch, with 
% same position of binding, N = 16, but with a base 
% pair mismatch within the first 6 positions

normalization_2_bm = 1/(((N-number_1-1+2))^3);
normalization_2_bmis = 1/((number_1+2)^3);

% construct the Hamiltonian. The T_ij factor for
% the single mismatch comes from the fact that
% (16-5)/16 = 11/16 approx 0.69

Hamiltonian_C2 = exp(-((normalization_2_bm * (((16-1)+(16-2)+(16-3)+(16-4)+(16-6)+(16-7)+(16-8)+(16-9)+(16-10)+(16-11)+(16-12)+(16-13)+(16-14)+(16-15))))+(normalization_2_bmis .*(((16-5)-1)* 0.69))));

% generate a test function with specific normalizations

func_C2 = @(X) (Hamiltonian_C2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate THIRD CASE for base pair mismatch, with
% same position of binding, N = 16, but with TWO base
% pair mismatches, BOTH in the first 6 positions

normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-16|/16 = 0.875, while the T_ij factor for the
% SECOND mismatch, at 5, comes from the fact that
% |5-16|/16 = 0.6875.

Hamiltonian_C3 = exp(-((normalization_3_bm * (((16-1)+(16-3)+(16-4)+(16-6)+(16-7)+(16-8)+(16-9)+(16-10)+(16-11)+(16-12)+(16-13)+(16-14)+(16-15))))+(normalization_3_bmis .*((((16-2)-1)* 0.875)+ (((16-5)-1)* 0.6875)))));

% generate a test function with specific normalizations

func_C3 = @(X) (Hamiltonian_C3)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate FOURTH CASE for base pair mismatch, with
% same position of binding, N = 16, but with TWO base
% pair mismatches, within the first 6 positions, that
% are subsequently right next to each other

normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-16|/16 = 0.75, while the T_ij factor for the
% SECOND mismatch, at 5, comes from the fact that
% |5-16|/16 = 0.6875.

Hamiltonian_C4 = exp(-((normalization_3_bm * (((16-1)+(16-3)+(16-4)+(16-6)+(16-7)+(16-8)+(16-9)+(16-10)+(16-11)+(16-12)+(16-13)+(16-14)+(16-15))))+(normalization_3_bmis .*((((16-4)-1)* 0.75)+ (((16-5)-1)* 0.6875)))));

% generate a test function with specific normalizations

func_C4 = @(X) (Hamiltonian_C4)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate FIFTH CASE for base pair mismatch, with 
% same position of binding, N = 16, but with TWO
% base pair mismatches, past the first 6 base pairs

normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 8,comes from the fact that
% |8-16|/16 = 0.50, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |10-16|/16 = 0.375.

Hamiltonian_C5 = exp(-((normalization_3_bm * (((16-1)+(16-3)+(16-4)+(16-6)+(16-7)+(16-8)+(16-9)+(16-10)+(16-11)+(16-12)+(16-13)+(16-14)+(16-15))))+(normalization_3_bmis .*((((16-8)-1)* 0.5)+ (((16-10)-1)* 0.375)))));

% generate a test function with specific normalizations

func_C5 = @(X) (Hamiltonian_C5)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));



% simulate SIXTH CASE for base pair mismatch, with 
% same position of binding, N = 16, but with TWO
% base pair mismatches, past the first 6 base pairs

normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 14,comes from the fact that
% |14-16|/16 = 0.125, while the T_ij factor for the
% SECOND mismatch, at 15, comes from the fact that
% |15-16|/16 = 0.0625.

Hamiltonian_C6 = exp(-((normalization_3_bm * (((16-1)+(16-2)+(16-3)+(16-4)+(16-6)+(16-7)+(16-8)+(16-9)+(16-10)+(16-11)+(16-12)+(16-13)+(16-8)+(16-10))))+(normalization_3_bmis .*((((16-14)-1)* 0.125)+ (((16-15)-1)* 0.0625)))));

% generate a test function with specific normalizations

func_C6 = @(X) (Hamiltonian_C6)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));



% REPEAT SIMILAR CASES AS ABOVE, BUT NOW, PUSH THE POSITION
% OF BINDING FURTHER UP THE SEQUENCE TO N = 13


% simulate FIRST CASE for base pair mismatch, with position of 
% binding, N = 13, with a SINGLE base pair mismatch past the 
% first 6 positions

normalization_1_bm = 1/(((N-number_1-1+2))^3);
normalization_1_bmis = 1/(((number_1+2))^3);

% construct the Hamiltonian. 0.5 comes from the
% T_ij factor for the single mismatch, |13-8|/13 = 0.3846 (!!!)

Hamiltonian_C1_2 = exp(-((normalization_1_bm * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_1_bmis *((((13-8)-1)*(0.3846))))));

% generate a test function with specific normalizations

func_C1_2 = @(X) (Hamiltonian_C1_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate SECOND CASE for base pair mismatch, with 
% same position of binding, N = 13, but with a base 
% pair mismatch within the first 6 positions

normalization_2_bm = 1/(((N-number_1-1+2))^3);
normalization_2_bmis = 1/((number_1+2)^3);

% construct the Hamiltonian. The T_ij factor for
% the single mismatch comes from the fact that
% |5-13|/13 = 11/16 approx 0.615

Hamiltonian_C2_2 = exp(-((normalization_2_bm * (((13-1)+(13-2)+(13-3)+(13-4)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_2_bmis .*((((13-5)-1)* 0.615)))));

% generate a test function with specific normalizations

func_C2_2 = @(X) (Hamiltonian_C2_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate THIRD CASE for base pair mismatch, with
% same position of binding, N = 13, but with TWO base
% pair mismatches, BOTH in the first 6 positions

normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-13|/13 = 0.846, while the T_ij factor for the
% SECOND mismatch, at 5, comes from the fact that
% |5-13|/13 = 0.61538.

Hamiltonian_C3_2 = exp(-((normalization_3_bm * (((13-1)+(13-3)+(13-4)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_3_bmis .*(((13-2)-1)* 0.846)+ (((13-5)-1)* 0.61538))));

% generate a test function with specific normalizations

func_C3_2 = @(X) (Hamiltonian_C3_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate FOURTH CASE for base pair mismatch, with
% same position of binding, N = 16, but with TWO base
% pair mismatches, within the first 6 positions, that
% are subsequently right next to each other

normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-13|/13 = 0.692, while the T_ij factor for the
% SECOND mismatch, at 5, comes from the fact that
% |5-13|/13 = 0.6154.

Hamiltonian_C4_2 = exp(-((normalization_3_bm * (((13-1)+(13-2)+(13-3)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_3_bmis .*((((16-4)-1)* 0.692))+ (((16-5)-1)* 0.6154))));

% generate a test function with specific normalizations

func_C4_2 = @(X) (Hamiltonian_C4_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate FIFTH CASE for base pair mismatch, with 
% same position of binding, N = 16, but with TWO
% base pair mismatches, past the first 6 base pairs

normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 8,comes from the fact that
% |8-13|/13 = 0.3846, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |10-13|/13 = 0.375.

Hamiltonian_C5_2 = exp(-((normalization_3_bm * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_3_bmis .*((((13-8)-1)* 0.3846))+ (((13-10)-1)* 0.2307))));

% generate a test function with specific normalizations

func_C5_2 = @(X) (Hamiltonian_C5_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate SIXTH CASE for base pair mismatch, with 
% same position of binding, N = 13, but with TWO
% base pair mismatches, past the first 6 base pairs

normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 14,comes from the fact that
% |10-13|/13 = 0.2307, while the T_ij factor for the
% SECOND mismatch, at 15, comes from the fact that
% |11-13|/13 = 0.1538.

Hamiltonian_C6_2 = exp(-((normalization_3_bm * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-12))))+(normalization_3_bmis .*((((13-10)-1)* 0.2307))+ (((13-11)-1)* 0.1538))));

% generate a test function with specific normalizations

func_C6_2 = @(X) (Hamiltonian_C6_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));



% SIMULATE MORE TARGET SEQUENCE CASES BY
% CHANGING THE POSITION OF BINDING TO N=17


% simulate FIRST CASE for base pair mismatch, with position of 
% binding, N = 17, with a SINGLE base pair mismatch past the 
% first 6 positions

normalization_1_bm = 1/(((N-number_1-1+2))^3);
normalization_1_bmis = 1/(((number_1+2))^3);

% construct the Hamiltonian. 0.5 comes from the
% T_ij factor for the single mismatch, |8-17|/17 = 0.529 (!!!)

Hamiltonian_C3_1 = exp(-((normalization_1_bm * (((17-1)+(17-2)+(17-3)+(17-4)+(17-5)+(17-6)+(17-7)+(17-9)+(17-10)+(17-11)+(17-12)+(17-13)+(17-14)+(17-15)+(17-16))))+(normalization_1_bmis *(((17-8)-1)*(0.529)))));

% generate a test function with specific normalizations

func_C3_1 = @(X) (Hamiltonian_C3_1)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate SECOND CASE (FROM THE THIRD CASE ABOVE) 
% for base pair mismatch, with same position of 
% binding, N = 17, but with TWO base
% pair mismatches, BOTH in the first 6 positions

normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-17|/17 = 0.8823, while the T_ij factor for the
% SECOND mismatch, at 5, comes from the fact that
% |5-17|/17 = 0.706.

Hamiltonian_C3_2 = exp(-((normalization_3_bm * (((17-1)+(17-3)+(17-4)+(17-6)+(17-7)+(17-8)+(17-9)+(17-10)+(17-11)+(17-12)+(17-13)+(17-14)+(17-15))))+(normalization_3_bmis .*(((17-2)-1)* 0.8823)+ ((17-5)-1)* 0.706)));

% generate a test function with specific normalizations

func_C3_2 = @(X) (Hamiltonian_C3_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));



% simulate THIRD CASE (FROM FIFTH CASE ABOVE) for base pair mismatch, with 
% same position of binding, N = 17, but with TWO
% base pair mismatches, past the first 6 base pairs

normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 8,comes from the fact that
% |8-17|/17 = 0.529, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |10-17|/17 = 0.4118.

Hamiltonian_C3_3 = exp(-((normalization_3_bm * (((17-1)+(17-2)+(17-3)+(17-4)+(17-5)+(17-6)+(17-7)+(17-9)+(17-10)+(17-11)+(17-12))))+(normalization_3_bmis .*((((17-8)-1)* 0.529)+ (((17-10)-1)* 0.4118)))));

% generate a test function with specific normalizations

func_C3_3 = @(X) (Hamiltonian_C3_3)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));



% WITH 15 SEQUENCES SO FAR, GENERATE 30 MORE. FOR 
% THE FIRST 15 ADDITIONAL SEQUENCES BELOW, USE
% THE 3 BASE PAIR MISMATCH CASES TO SHOW THAT
% BINDING CANNOT OCCUR.

% simulate FIRST case of 3 base pair
% mismatch, in which we have 3 base pair
% mismatches (2,3,4), consecutively all within
% the first 6 positions, with N = 13

normalization_triple_bm = 1/(10+2)^3;
normalization_triple_bmis = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 8,comes from the fact that
% |2-13|/13 = 0.846, while the T_ij factor for the
% SECOND mismatch, at 3, comes from the fact that
% |3-13|/13 = 0.769, and the THIRD mismatch factor
% comes from |4-13|/13 = 0.6923.

Hamiltonian_triple = exp(-((normalization_triple_bm * (((13-1)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_triple_bmis .*((((13-2)-1)* 0.846)+ (((13-3)-1)* 0.769)+(((13-4)-1)* 0.6923)))));

% generate a test function with specific normalizations

func_triple = @(X) (Hamiltonian_triple)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate SECOND case of 3 base pair mismatch,
% in which we have 3 base pair mismatches
% as in the first case, BUT we translate 
% the sequence of mismatches down the sequence
% by one unit, again with N = 13

normalization_triple_bm_1 = 1/(11+2)^3;
normalization_triple_bmis_1 = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 8,comes from the fact that
% |4-13|/13 = 0.6923, while the T_ij factor for the
% SECOND mismatch, at 5, comes from the fact that
% |5-13|/13 = 0.615, and the THIRD T_ij factor comes
% from the fact that |6-13|/13 = 0.5384.

Hamiltonian_triple_2_A = exp(-((normalization_triple_bm_1 * (((13-1)+(13-2)+(13-3)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_triple_bmis_1 .*((((13-4)-1)* 0.6923)+ (((13-5)-1)* 0.615)+(((13-6)-1)* 0.5384)))));

% generate a test function with specific normalizations

func_triple_2_A = @(X) (Hamiltonian_triple_2_A)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate THIRD case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 13

normalization_double_1 = 1/(10+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-13|/13 = 0.846, while the T_ij factor for the
% SECOND mismatch, at 4, comes from the fact that
% |4-13|/13 = 0.692

Hamiltonian_double = exp(-((normalization_double_1 * (((13-1)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2 .*((((13-2)-1)* 0.846)+ (((13-4)-1)* 0.692)))));

% generate a test function with specific normalizations

func_double_2 = @(X) (Hamiltonian_double)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate FOURTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 13

normalization_double_1 = 1/(10+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-13|/13 = 0.692, while the T_ij factor for the
% SECOND mismatch, at 6, comes from the fact that
% |6-13|/13 = 0.538

Hamiltonian_double_C2 = exp(-((normalization_double_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2 .*((((13-4)-1)* 0.692)+ (((13-6)-1)* 0.538)))));

% generate a test function with specific normalizations

func_double_C2 = @(X) (Hamiltonian_double_C2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate FIFTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 13

normalization_double_1 = 1/(10+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-13|/13 = 0.846, while the T_ij factor for the
% SECOND mismatch, at 7, comes from the fact that
% |7-13|/13 = 0.462

Hamiltonian_double_C2_2 = exp(-((normalization_double_1 * (((13-1)+(13-3)+(13-4)+(13-5)+(13-6)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2 .*((((13-2)-1)* 0.846)+ (((13-7)-1)* 0.462)))));

% generate a test function with specific normalizations

func_double_C2_2 = @(X) (Hamiltonian_double_C2_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate SIXTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 13

normalization_double_1 = 1/(10+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-13|/13 = 0.692, while the T_ij factor for the
% SECOND mismatch, at 7, comes from the fact that
% |7-13|/13 = 0.462

Hamiltonian_double_C2_2A = exp(-((normalization_double_1 * (((13-1)+(13-3)+(13-2)+(13-5)+(13-6)+(13-8)+(13-9)+(13-10)+(13-12)+(13-11))))+(normalization_double_2 .*((((13-4)-1)* 0.692)+ (((13-7)-1)* 0.462)))));

% generate a test function with specific normalizations

func_double_C2_2A = @(X) (Hamiltonian_double_C2_2A)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate SEVENTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 13

normalization_double_1 = 1/(10+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-13|/13 = 0.846, while the T_ij factor for the
% SECOND mismatch, at 8, comes from the fact that
% |8-13|/13 = 0.385

Hamiltonian_double_C2_2B = exp(-((normalization_double_1 * (((13-1)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-9)+(13-10)+(13-12)+(13-11))))+(normalization_double_2 .*((((13-2)-1)* 0.846)+ (((13-8)-1)* 0.385)))));

% generate a test function with specific normalizations

func_double_C2_2B = @(X) (Hamiltonian_double_C2_2B)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate EIGHTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 13

normalization_double_1 = 1/(10+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-13|/13 = 0.846, while the T_ij factor for the
% SECOND mismatch, at 9, comes from the fact that
% |9-13|/13 = 0.308

Hamiltonian_double_C2_2_2 = exp(-((normalization_double_1 * (((13-1)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-8)+(13-10)+(13-12)+(13-11))))+(normalization_double_2 .*((((13-2)-1)* 0.846)+ (((13-9)-1)* 0.308)))));

% generate a test function with specific normalizations

func_double_C2_2_2 = @(X) (Hamiltonian_double_C2_2_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate NINTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 13

normalization_double_1 = 1/(10+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-13|/13 = 0.692, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |10-13|/13 = 0.231

Hamiltonian_double_C2_2_2_2 = exp(-((normalization_double_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-12)+(13-11))))+(normalization_double_2 .*((((13-4)-1)* 0.692)+ (((13-10)-1)* 0.231)))));

% generate a test function with specific normalizations

func_double_C2_2_2_2 = @(X) (Hamiltonian_double_C2_2_2_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate TENTH case of 3 base pair mismatch,
% in which we have ONE base pair mismatch
% within the first 6 position, while we
% also have 2 additional mismatches past
% the first 6 positions, for a 
% position of binding, N = 13

normalization_double_1 = 1/(9+2)^3;
normalization_double_2 = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-13|/13 = 0.692, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |10-13|/13 = 0.231, while finally, the THIRD
% T_ij factor comes from the fact that 
% 0.15 = |11-13|/13 

Hamiltonian_triple_A = exp(-((normalization_double_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-12))))+(normalization_double_2 .*((((13-4)-1)* 0.6923)+ (((13-10)-1)* 0.231)+(((13-11)-1)* 0.15)))));

% generate a test function with specific normalizations

func_triple_A = @(X) (Hamiltonian_triple_A)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate ELEVENTH case of 3 base pair mismatch,
% in which we have ONE base pair mismatch
% within the first 6 position, while we
% also have 2 additional mismatches past
% the first 6 positions, for a 
% position of binding, N = 13

normalization_double_1 = 1/(9+2)^3;
normalization_double_2 = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-13|/13 = 0.692, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |11-13|/13 = 0.153, while finally, the THIRD
% T_ij factor comes from the fact that 
% 0.077 = |12-13|/13 

Hamiltonian_triple_A_2 = exp(-((normalization_double_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10))))+(normalization_double_2 .*((((13-4)-1)* 0.6923)+ (((13-11)-1)* 0.153)+(((13-12)-1)* 0.077)))));

% generate a test function with specific normalizations

func_triple_A_2 = @(X) (Hamiltonian_triple_A_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate TWELFTH case of 4 base pair mismatch,
% in which we have TWO base pair mismatchES
% within the first 6 position, while we
% also have 2 additional mismatches past
% the first 6 positions, for a 
% position of binding, N = 13

normalization_double_1 = 1/(8+2)^3;
normalization_double_2 = 1/(2+4)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2, comes from the fact that
% the quantity T_ij is equal to |2-13|/13 = 0.846. While 
% for the SECOND mismatch, at 4,comes from the fact that
% |4-13|/13 = 0.692, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |11-13|/13 = 0.153, while finally, the THIRD
% T_ij factor comes from the fact that 
% 0.077 = |12-13|/13, and FINALLY, the 
% FOURTH mismatch factor comes from
% the fact that 1/13 = 0.0769

Hamiltonian_quad = exp(-((normalization_double_1 * (((13-1)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10))))+(normalization_double_2 .*((((13-2)-1)* 0.846)+(((13-4)-1)* 0.6923)+ (((13-11)-1)* 0.153)+(((13-12)-1)* 0.0769)))));

% generate a test function with specific normalizations

func_quad = @(X) (Hamiltonian_quad)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate THIRTEENTH case of 3 base pair mismatch,
% in which we have 3 base pair mismatches, 
% each of which are past the first 6 positions
% with position of binding, N = 13

normalization_double_1 = 1/(9+2)^3;
normalization_double_2 = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 10, comes from the fact that
% |10-13|/13 = 0.231, while finally, the SECOND
% T_ij factor comes from the fact that 
% 0.153 = |11-13|/13, and FINALLY, the 
% THIRD mismatch factor comes from
% the fact that 1/13 = 0.0769

Hamiltonian_triple_II = exp(-((normalization_double_1 * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9))))+(normalization_double_2 .*((((13-10)-1)* 0.231)+(((13-11)-1)* 0.153)+ (((13-12)-1)* 0.0769)))));

% generate a test function with specific normalizations

func_triple_II = @(X) (Hamiltonian_triple_II)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate FOURTEENTH case of 2 base pair mismatch,
% in which we have 3 base pair mismatches, 
% each of which are past the first 6 positions
% with position of binding, N = 13

normalization_double_1 = 1/(10+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 11, comes from the fact that
% |10-13|/13 = 0.153, while finally, the SECOND
% T_ij factor comes from the fact that 
% 0.0769 = |12-13|/13, at 12.

Hamiltonian_double_II = exp(-((normalization_double_1 * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10))))+(normalization_double_2 .*((((13-11)-1)* 0.153)+(((13-12)-1)* 0.0769)))));

% generate a test function with specific normalizations

func_double_II = @(X) (Hamiltonian_double_II)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate FIFTEENTH case of 2 base pair mismatch,
% in which we have 3 base pair mismatches, 
% each of which are past the first 6 positions
% with position of binding, N = 13

normalization_double_1 = 1/(10+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 9, comes from the fact that
% |9-13|/13 = 0.308, while the SECOND
% T_ij factor comes from the fact that 
% 0.231 = |10-13|/13, at 10, and finally,
% the THIRD mismatch at 11 comes from the 
% fact that |11-13|/13 = 0.154.

Hamiltonian_TRIPLE_II = exp(-((normalization_double_1 * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-8))))+(normalization_double_2 .*((((13-9)-1)* 0.308)+(((13-10)-1)* 0.231) + (((13-11)-1)*0.154)))));

% generate a test function with specific normalizations

func_TRIPLE_II = @(X) (Hamiltonian_TRIPLE_II)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% NOW, WE WILL GENERATE THE LAST 15 SEQUENCES, FOR THE FINAL THIRD PLOT, 
% BY TAKING ALL OF THE MUTATIONS IN THE PREVIOUS 15 AND 
% CHANGING THE POSITION OF BINDING. ALL OF THE BASE
% MISMATCH PATTERNS ARE THE SAME. 

% THE BINDING POSITION WILL BE
% CHANGED TO N = 15.

% simulate FIRST case of 3 base pair
% mismatch, in which we have 3 base pair
% mismatches (2,3,4), consecutively all within
% the first 6 positions, with N = 15

normalization_triple_bm = 1/(11+2)^3;
normalization_triple_bmis = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 8,comes from the fact that
% |2-15|/15 = 0.867, while the T_ij factor for the
% SECOND mismatch, at 3, comes from the fact that
% |3-15|/15 = 0.8, and the THIRD mismatch factor
% comes from |4-15|/15 = 0.733

Hamiltonian_triple_2_B = exp(-((normalization_triple_bm * (((15-1)+(15-5)+(15-6)+(15-7)+(15-8)+(15-9)+(15-10)+(15-11)+(15-12))))+(normalization_triple_bmis .*((((15-2)-1)* 0.867)+ (((15-3)-1)* 0.800)+(((15-4)-1)* 0.733)))));

% generate a test function with specific normalizations

func_triple_2_B = @(X) (Hamiltonian_triple_2_B)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate SECOND case of 3 base pair mismatch,
% in which we have 3 base pair mismatches
% as in the first case, BUT we translate 
% the sequence of mismatches down the sequence
% by one unit, again with N = 15

normalization_triple_bm_1 = 1/(11+2)^3;
normalization_triple_bmis_1 = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 8,comes from the fact that
% |4-15|/15 = 0.733, while the T_ij factor for the
% SECOND mismatch, at 5, comes from the fact that
% |5-15|/15 = 2/3, and the THIRD T_ij factor comes
% from the fact that |6-15|/15 = 0.6

Hamiltonian_triple_2_2A = exp(-((normalization_triple_bm_1 * (((15-1)+(15-2)+(15-3)+(15-7)+(15-8)+(15-9)+(15-10)+(15-11)+(15-12))))+(normalization_triple_bmis_1 .*((((15-4)-1)* 0.733)+ (((15-5)-1)* 0.667)+(((15-6)-1)* 0.600)))));

% generate a test function with specific normalizations

func_triple_2_2A = @(X) (Hamiltonian_triple_2_2A)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate THIRD case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 15

normalization_double_1 = 1/(12+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-15|/15 = 0.867, while the T_ij factor for the
% SECOND mismatch, at 4, comes from the fact that
% |4-15|/15 = 0.733

Hamiltonian_double_2 = exp(-((normalization_double_1 * (((15-1)+(15-3)+(15-5)+(15-6)+(15-7)+(15-8)+(15-9)+(15-10)+(15-11)+(15-12))))+(normalization_double_2 .*((((15-2)-1)* 0.867)+ (((15-4)-1)* 0.733)))));

% generate a test function with specific normalizations

func_double_2_2 = @(X) (Hamiltonian_double_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate FOURTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 15

normalization_double_1 = 1/(12+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-15|/15 = 0.733, while the T_ij factor for the
% SECOND mismatch, at 6, comes from the fact that
% |6-15|/15 = 0.6

Hamiltonian_double_C2_2 = exp(-((normalization_double_1 * (((15-1)+(15-2)+(15-3)+(15-5)+(15-7)+(15-8)+(15-9)+(15-10)+(15-11)+(15-12))))+(normalization_double_2 .*((((15-4)-1)* 0.733)+ (((15-6)-1)* 0.600)))));

% generate a test function with specific normalizations

func_double_C2_2 = @(X) (Hamiltonian_double_C2_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate FIFTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 15

normalization_double_1 = 1/(12+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-15|/15 = 0.867, while the T_ij factor for the
% SECOND mismatch, at 7, comes from the fact that
% |7-15|/15 = 0.533

Hamiltonian_double_C2_2_2A = exp(-((normalization_double_1 * (((15-1)+(15-3)+(15-4)+(15-5)+(15-6)+(15-8)+(15-9)+(15-10)+(15-11)+(15-12))))+(normalization_double_2 .*((((15-2)-1)* 0.867)+ (((15-7)-1)* 0.533)))));

% generate a test function with specific normalizations

func_double_C2_2_2A = @(X) (Hamiltonian_double_C2_2_2A)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% simulate SIXTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 15

normalization_double_1 = 1/(12+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-15|/15 = 0.733, while the T_ij factor for the
% SECOND mismatch, at 7, comes from the fact that
% |7-15|/15 = 0.533

Hamiltonian_double_C2_2A_2A = exp(-((normalization_double_1 * (((15-1)+(15-3)+(15-2)+(15-5)+(15-6)+(15-8)+(15-9)+(15-10)+(15-12)+(15-11))))+(normalization_double_2 .*((((15-4)-1)* 0.733)+ (((15-7)-1)* 0.533)))));

% generate a test function with specific normalizations

func_double_C2_2A_2A = @(X) (Hamiltonian_double_C2_2A_2A)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate SEVENTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 15

normalization_double_1 = 1/(12+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-15|/15 = 0.867, while the T_ij factor for the
% SECOND mismatch, at 8, comes from the fact that
% |8-15|/15 = 0.467

Hamiltonian_double_C2_2B_2 = exp(-((normalization_double_1 * (((15-1)+(15-3)+(15-4)+(15-5)+(15-6)+(15-7)+(15-9)+(15-10)+(15-12)+(15-11))))+(normalization_double_2 .*((((15-2)-1)* 0.867)+ (((15-8)-1)* 0.467)))));

% generate a test function with specific normalizations

func_double_C2_2B_2 = @(X) (Hamiltonian_double_C2_2B_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate EIGHTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 15

normalization_double_1 = 1/(12+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2,comes from the fact that
% |2-15|/15 = 0.867, while the T_ij factor for the
% SECOND mismatch, at 9, comes from the fact that
% |9-15|/15 = 0.4

Hamiltonian_double_C2_2_2_2 = exp(-((normalization_double_1 * (((15-1)+(15-3)+(15-4)+(15-5)+(15-6)+(15-7)+(15-8)+(15-10)+(15-12)+(15-11))))+(normalization_double_2 .*((((15-2)-1)* 0.867)+ (((15-9)-1)* 0.400)))));

% generate a test function with specific normalizations

func_double_C2_2_2_2_2 = @(X) (Hamiltonian_double_C2_2_2_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate NINTH case of 2 base pair mismatch,
% in which we have 2 base pair mismatches, at
% position of binding N = 15

normalization_double_1 = 1/(12+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-15|/15 = 0.733, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |10-15|/15 = 1/3

Hamiltonian_double_C2_2_2_2_2A = exp(-((normalization_double_1 * (((15-1)+(15-2)+(15-3)+(15-5)+(15-6)+(15-7)+(15-8)+(15-9)+(15-12)+(15-11))))+(normalization_double_2 .*((((15-4)-1)* 0.733)+ (((15-10)-1)* 0.333)))));

% generate a test function with specific normalizations

func_double_C2_2_2_2_2A = @(X) (Hamiltonian_double_C2_2_2_2_2A)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate TENTH case of 3 base pair mismatch,
% in which we have ONE base pair mismatch
% within the first 6 position, while we
% also have 2 additional mismatches past
% the first 6 positions, for a 
% position of binding, N = 15

normalization_double_1 = 1/(11+2)^3;
normalization_double_2 = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-15|/15 = 0.733, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |10-15|/15 = 1/3, while finally, the THIRD
% T_ij factor comes from the fact that 
% 0.267 = |11-15|/15 

Hamiltonian_triple_A_2 = exp(-((normalization_double_1 * (((15-1)+(15-2)+(15-3)+(15-5)+(15-6)+(15-7)+(15-8)+(15-9)+(15-12))))+(normalization_double_2 .*((((15-4)-1)* 0.733)+ (((15-10)-1)* 0.333)+(((15-11)-1)* 0.267)))));

% generate a test function with specific normalizations

func_triple_A_2 = @(X) (Hamiltonian_triple_A_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate ELEVENTH case of 3 base pair mismatch,
% in which we have ONE base pair mismatch
% within the first 6 position, while we
% also have 2 additional mismatches past
% the first 6 positions, for a 
% position of binding, N = 15

normalization_double_1 = 1/(11+2)^3;
normalization_double_2 = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-15|/15 = 0.733, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |11-15|/15 = 0.267, while finally, the THIRD
% T_ij factor comes from the fact that 
% 1/5 = |12-15|/15 

Hamiltonian_triple_A_2_2 = exp(-((normalization_double_1 * (((15-1)+(15-2)+(15-3)+(15-5)+(15-6)+(15-7)+(15-8)+(15-9)+(15-10))))+(normalization_double_2 .*((((15-4)-1)* 0.733)+ (((15-11)-1)* 0.267)+(((15-12)-1)* 0.200)))));

% generate a test function with specific normalizations

func_triple_A_2_2 = @(X) (Hamiltonian_triple_A_2_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate TWELFTH case of 4 base pair mismatch,
% in which we have TWO base pair mismatchES
% within the first 6 position, while we
% also have 2 additional mismatches past
% the first 6 positions, for a 
% position of binding, N = 15

normalization_double_1 = 1/(10+2)^3;
normalization_double_2 = 1/(2+4)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 2, comes from the fact that
% the quantity T_ij is equal to |2-15|/15 = 0.867. While 
% for the SECOND mismatch, at 4,comes from the fact that
% |4-15|/15 = 0.733, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |11-15|/15 = 0.267, while finally, the THIRD
% T_ij factor comes from the fact that 
% 1/5 = |12-15|/15, and FINALLY, the 
% FOURTH mismatch factor comes from
% the fact that 3/15 = 1/5

Hamiltonian_quad_2 = exp(-((normalization_double_1 * (((15-1)+(15-3)+(15-5)+(15-6)+(15-7)+(15-8)+(15-9)+(15-10))))+(normalization_double_2 .*((((15-2)-1)* 0.867)+(((15-4)-1)* 0.733)+ (((15-11)-1)* 0.267)+(((15-12)-1)* 0.200)))));

% generate a test function with specific normalizations

func_quad_2 = @(X) (Hamiltonian_quad_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate THIRTEENTH case of 3 base pair mismatch,
% in which we have 3 base pair mismatches, 
% each of which are past the first 6 positions
% with position of binding, N = 15

normalization_double_1 = 1/(11+2)^3;
normalization_double_2 = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 10, comes from the fact that
% |10-15|/15 = 1/3 , while finally, the SECOND
% T_ij factor comes from the fact that 
% 0.267 = |11-15|/15, and FINALLY, the 
% THIRD mismatch factor comes from
% the fact that 3/15 = 1/5

Hamiltonian_triple_II_2 = exp(-((normalization_double_1 * (((15-1)+(15-2)+(15-3)+(15-4)+(15-5)+(15-6)+(15-7)+(15-8)+(15-9))))+(normalization_double_2 .*((((15-10)-1)* 0.333)+(((15-11)-1)* 0.267)+ (((15-12)-1)* 0.200)))));

% generate a test function with specific normalizations

func_triple_II_2 = @(X) (Hamiltonian_triple_II_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate FOURTEENTH case of 2 base pair mismatch,
% in which we have 3 base pair mismatches, 
% each of which are past the first 6 positions
% with position of binding, N = 15

normalization_double_1 = 1/(12+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 11, comes from the fact that
% |10-15|/15 = 1/3, while finally, the SECOND
% T_ij factor comes from the fact that 
%  1/5 = |12-15|/15, at 12.

Hamiltonian_double_II_2 = exp(-((normalization_double_1 * (((15-1)+(15-2)+(15-3)+(15-4)+(15-5)+(15-6)+(15-7)+(15-8)+(15-9)+(15-10))))+(normalization_double_2 .*((((15-11)-1)* 0.333)+(((15-12)-1)* 0.200)))));

% generate a test function with specific normalizations

func_double_II_2 = @(X) (Hamiltonian_double_II_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% simulate FIFTEENTH case of 2 base pair mismatch,
% in which we have 3 base pair mismatches, 
% each of which are past the first 6 positions
% with position of binding, N = 15

normalization_double_1 = 1/(12+2)^3;
normalization_double_2 = 1/(2+2)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 9, comes from the fact that
% |9-15|/15 = 0.400, while the SECOND
% T_ij factor comes from the fact that 
% 1/3 = |10-15|/15, at 10, and finally,
% the THIRD mismatch at 11 comes from the 
% fact that |11-15|/15 = 0.267.

Hamiltonian_TRIPLE_II_2 = exp(-((normalization_double_1 * (((15-1)+(15-2)+(15-3)+(15-4)+(15-5)+(15-6)+(15-7)+(15-8))))+(normalization_double_2 .*((((15-9)-1)* 0.400)+(((15-10)-1)* 0.333) + (((15-11)-1)*0.267)))));

% generate a test function with specific normalizations

func_TRIPLE_II_2 = @(X) (Hamiltonian_TRIPLE_II_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


% PLOTTING

subplot(3,1,1)

sequences_Hamiltonian_FIRSTfifteen = [Hamiltonian_C1; Hamiltonian_C2; Hamiltonian_C3; Hamiltonian_C4; Hamiltonian_C5; Hamiltonian_C6;  Hamiltonian_C1_2; Hamiltonian_C2_2; Hamiltonian_C3_2; Hamiltonian_C4_2; Hamiltonian_C5_2; Hamiltonian_C6_2; Hamiltonian_C3_1; Hamiltonian_C3_2; Hamiltonian_C3_3];
sorted = sort([0.854726811659131 0.754147819067206 0.731027355324388 0.767607861569519 0.893925783791197 0.963186056694457 0.927929013307228 0.838170039325183 0.011593287855590 0.001852884668431 0.296059374682612 0.833068835454632 0.828824704709133 0.011593287855589 0.894347497534619]);
x_3 = linspace(1000,15000,15);
Hamiltonian_C5_2 = exp(-((normalization_3_bm * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_3_bmis .*((((13-8)-1)* 0.3846))+ (((13-10)-1)* 0.2307))));
% x = transpose(x);
y_3 = sequences_Hamiltonian_FIRSTfifteen;
plot(x_3,y_3,'or')
xticks([1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 15000]);
xticklabels({ 'AATTTA,ACAGAG,CaGTCAGAAnATCCACGTA ' , ' AATTTA,ACAGaG,CAGTCAGAAnAATCCACGTA' , ' AATTTA,AcAGaG,CAGTCAGAAnAATCCACGTA ' , ' AATTTA,ACAgaG,CAGTCAGAAnATCCACGTA ' , ' AATTTA,ACAGAG,CaGtCAGAAnTCCACGTA' , '{\color{orange}AATTTA,ACAgAG,CAGTCaanAAATCCACGTA}', '{\color{orange}AATTTA,ACAGAG,CaGTCAnAAAAATCCACGTA}' , '{\color{orange}AATTTA,ACAGaG,CAGTCAnAAAAATCCACGTA}' , ' AATTTA,AcAgaG,CAGtCAnAAAAATCCACGTA' , ' AATTTA,ACAgaG,CAGtCAAAnAAAATCCACGTA' , ' AATTTA,ACAGAG,CaGtCAAAAnAAATCCACGTA', 'AATTTA,AcAgAG,CAgtCAnAAAAATCCACGTA ' , ' AATTTA,ACAGAG,CaGTCAGAAnATCCACGTA ' , '  AATTTA,AcAGaG,CAGTCAGAAnAATCCACGTA ' , '{\color{orange}AATTTA,ACAGAG,CaGtCAGAAnAATCCACGTA}'});
xtickangle(35)
yticks(sorted)
title('CASE 1: Confirming that the Hamiltonian assigns {\color{blue}less energy} to sequences with base mismatches from sequences of Type I, II or III')
xlabel('TARGET SEQUENCE')
ylabel('HAMILTONIAN')
legend({'Hamiltonian value corresponding to a sequence'})
hold on;
refline

subplot(3,1,2)

sequences_Hamiltonian_SECONDfifteen = [Hamiltonian_triple; Hamiltonian_triple_2_A;Hamiltonian_double;Hamiltonian_double_C2;Hamiltonian_double_C2_2;Hamiltonian_double_C2_2A;Hamiltonian_double_C2_2B;Hamiltonian_double_C2_2_2; Hamiltonian_double_C2_2_2_2; Hamiltonian_triple_A ; Hamiltonian_triple_A_2;Hamiltonian_quad; Hamiltonian_triple_II; Hamiltonian_double_II; Hamiltonian_TRIPLE_II];
sorted_values = sort([ 0.822725246606518 0.784810520925704 0.777048652445489 0.841291448057122 0.802988704367723 0.852951933232366 0.825199891201852 0.832698761282255 0.799266926400893 0.907295259867857 0.893125988846724 0.886463133649195 0.942692943836891 0.955239291590360 0.938554501313303]);
x = linspace(1000,14000,15);
% x = transpose(x);
y = sequences_Hamiltonian_SECONDfifteen;
plot(x,y,'or')
xticks([1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000]);
xticklabels({ 'AATTTA,AcagAG,CAGTCAnAAAAATCCACGTA ' , ' AATTTA,ACAgag,CAGTCAnAAAAATCCACGTA' , '{\color{orange}AATTTA,AcAgAG,CAGTCAGnAAAATCCACGTA}' , 'AATTTA,ACAgAg,CAGTCAnAAAAATCCACGTA ' , ' {\color{orange}AATTTA,AcAGAG,cAGTCAnAAAAATCCACGTA}' , ' AATTTA,ACAgAG,cAGTCAnAAAAATCCACGTA', ' {\color{orange}AATTTA,AcAGAG,cAGTCAnAAAAATCCACGTA}' , ' {\color{orange}AATTTA,AcAGAG,CAgTCAnAAAAATCCACGTA}' , ' AATTTA,ACAgAG,CAGtCAnAAAAATCCACGTA' , ' AATTTA,ACAgAG,CAGtcAnAAAAATCCACGTA' , ' AATTTA,ACAgAG,CAGTcanAAAAATCCACGTA', 'AATTTA,AcAgAG,CAGTcanAAAAATCCACGTA ' , ' AATTTA,ACAgAG,CAGtcanAAAAATCCACGTA ' , '  AATTTA,ACAgAG,CAGTcanAAAAATCCACGTA ' , ' AATTTA,ACAgAG,CAgtcAnAAAAATCCACGTA ' });
xtickangle(35)
yticks(sorted_values)
title('CASE 2: Confirming that the Hamiltonian assigns {\color{blue}less energy} to sequences with base mismatches from sequences of Type I, II or III')
xlabel('TARGET SEQUENCE') 
ylabel('HAMILTONIAN') 
legend({'Hamiltonian value corresponding to a sequence'})
hold on;
refline

% sequences_Hamiltonian_THIRDfifteen = [Hamiltonian_triple_2; Hamiltonian_triple_2_2A; Hamiltonian_double_2; Hamiltonian_double_C2_2; Hamiltonian_double_C2_2_2A; Hamiltonian_double_C2_2A_2A; Hamiltonian_double_C2_2B_2; Hamiltonian_double_C2_2_2_2; Hamiltonian_double_C2_2_2_2_2A; Hamiltonian_triple_A_2; Hamiltonian_triple_A_2_2; Hamiltonian_quad_2; Hamiltonian_triple_II_2; Hamiltonian_double_II_2 ; Hamiltonian_TRIPLE_II_2];

subplot(3,1,3)

sequences_Hamiltonian_THIRDfifteen = [Hamiltonian_triple_2_B; Hamiltonian_triple_2_2A; Hamiltonian_double_2; Hamiltonian_double_C2_2; Hamiltonian_double_C2_2_2A; Hamiltonian_double_C2_2A_2A; Hamiltonian_double_C2_2B_2; Hamiltonian_double_C2_2_2_2; Hamiltonian_double_C2_2_2_2_2A; Hamiltonian_triple_A_2; Hamiltonian_triple_A_2_2; Hamiltonian_quad_2; Hamiltonian_triple_II_2; Hamiltonian_double_II_2 ; Hamiltonian_TRIPLE_II_2];
sorted_values_1 = sort([0.784810520925704 0.837079725306906 0.736740607519727 0.802988704367723 0.778506198652469 0.816216256906491 0.789601309899058 0.799266926400893 0.846466326508041 0.893125988846724 0.898991269958390 0.879193736507195 0.940607860875021 0.945085039872948 0.909199833712599]);
x_1 = linspace(1000,14000,15);
% x = transpose(x);
y_1 = sequences_Hamiltonian_THIRDfifteen;
plot(x_1,y_1,'or')
xticks([1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000]);
xticklabels({ 'AATTTA,AcagAG,CAGTCAGAnAAATCCACGTA ' , ' AATTTA,ACAgag,CAGTCAGAnAAATCCACGTA' , '{\color{orange}AATTTA,AcAgAG,CAGTCAGGAnAATCCACGTA} ' , '{\color{orange}AATTTA,ACAgAg,CAGTCAGAnAAATCCACGTA} ' , ' AATTTA,AcAGAG,cAGTCAGAnAAATCCACGTA' , ' AATTTA,ACAgAG,cAGTCAGAnAAATCCACGTA', ' AATTTA,AcAGAG,cAGTCAnAAAAATCCACGTA' , ' AATTTA,AcAGAG,CAgTCAnAAAAATCCACGTA' , ' AATTTA,ACAgAG,CAGtCAnAAAAATCCACGTA' , ' AATTTA,ACAgAG,CAGtcAnAAAAATCCACGTA' , ' AATTTA,ACAgAG,CAGTcaGAnAAATCCACGTA', 'AATTTA,AcAgAG,CAGTcaGAnAAATCCACGTA ' , ' AATTTA,ACAgAG,CAGtcaGAnAAATCCACGTA ' , '  AATTTA,ACAgAG,CAGTcaGAnAAATCCACGTA ' , ' AATTTA,ACAgAG,CAgtcAGAnAAATCCACGTA ' });
xtickangle(35)
yticks(sorted_values_1)
title('CASE 3: Confirming that the Hamiltonian assigns {\color{blue}less energy} to sequences with base mismatches from sequences of Type I, II or III')
xlabel('TARGET SEQUENCE') 
ylabel('HAMILTONIAN') 
legend({'Hamiltonian value corresponding to a sequence'})
hold on;
refline

% SPARE CODE FOR OLD HAMILTONIANS AND TEST FUNCTION

% hamiltonian_test_13 = exp(-((normalization * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+ 0.5*(((13-7)-1)*(1-0.4))));
% hamiltonian_test_16 = exp(-((normalization * (((16-1)+(16-2)+(16-3)+(16-4)+(16-5)+(16-6)+(16-8)+(16-9)+(16-10)+(16-11)+(16-12))))+ 0.5*(((13-7)-1)*(1-0.4))));

% func_truncated_ii_test_2 = @(X) exp(hamiltonian_ii_test)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

end 





