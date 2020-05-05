function [visit_array_11, visit_array_111, visit_array_1111] = demo_2()

% [visit_array_1, visit_array_2, visit_array_3, Hamiltonian_array] = demo_2()

% [visit_array_1, visit_array_2, visit_array_3, visit_array_3_2, sequences_Hamiltonian_3bp_1 , sequences_Hamiltonian_3bp_2] = demo_2()

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
normalization_3_bm = 1/(13+2)^3;
normalization_3_bmis = 1/(2+2)^3;

% define possible cases of 3 mismatches
number_1 = 1;
number_2 = 2;
number_3 = 3;

J=3;
lambda_mismatch = exp(-10/2);
lambda_c = 0.025;
lambda_p = 0.015;
couplings_ij1 = exp(N-J);
couplings_ij2 = N-J;
% hamiltonian_i1 = exp(couplings_ij1 * 1 * exp(lambda_mismatch * (N-J)));
% hamiltonian_i2 = exp(couplings_ij2 * 1) * exp(lambda_mismatch * (N-J));

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

% REPEAT THE SAME SEQUENCES FOR A DIFFERENT POSITION OF BINDING,
% ORIGINALLY FROM THE THIRD SET OF 15 SEQUENCES IN THE DEMO .M FILE

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

Hamiltonian_double_C2_2_2_2X = exp(-((normalization_double_1 * (((15-1)+(15-3)+(15-4)+(15-5)+(15-6)+(15-7)+(15-8)+(15-10)+(15-12)+(15-11))))+(normalization_double_2 .*((((15-2)-1)* 0.867)+ (((15-9)-1)* 0.400)))));

% generate a test function with specific normalizations

func_double_C2_2_2_2_2 = @(X) (Hamiltonian_double_C2_2_2_2X)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

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

Hamiltonian_double_C2_2Y = exp(-((normalization_double_1 * (((15-1)+(15-2)+(15-3)+(15-5)+(15-7)+(15-8)+(15-9)+(15-10)+(15-11)+(15-12))))+(normalization_double_2 .*((((15-4)-1)* 0.733)+ (((15-6)-1)* 0.600)))));

% generate a test function with specific normalizations

func_double_C2_2 = @(X) (Hamiltonian_double_C2_2Y)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));


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


% THE ABOVE SEQUENCES CONCLUDE THE 2 BP GROUP
% CONTINUE WITH 3 BP GROUP, FOR EACH BINDING
% POSITION, ORIGINALLY FROM THE SECOND
% AND THIRD GROUP IN DEMO .M 

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

% CHANGE THE BINDING POSITION (!!!)


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

Hamiltonian_triplee_II_2 = exp(-((normalization_double_1 * (((15-1)+(15-2)+(15-3)+(15-4)+(15-5)+(15-6)+(15-7)+(15-8)+(15-9))))+(normalization_double_2 .*((((15-10)-1)* 0.333)+(((15-11)-1)* 0.267)+ (((15-12)-1)* 0.200)))));

% generate a test function with specific normalizations

func_triple_II_2 = @(X) (Hamiltonian_triplee_II_2)./(1+(lambda_c .* exp(-beta * Table{I,8})) + (lambda_c * exp(-beta * Table_4{JI,8})) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

% PLOTTING THE SEQUENCES SO FAR, SEPARATELY
% AND INDIVIDUALLY

figure(1)
% subplot(2,1,1)
% define sequences_Hamiltonian as in demo .m
sequences_Hamiltonian_2bp_1 = [Hamiltonian_double; Hamiltonian_double_C2_2; Hamiltonian_double_C2_2B; Hamiltonian_double_C2_2_2; Hamiltonian_double_C2; Hamiltonian_double_C2_2A; Hamiltonian_double_C2_2_2_2];  
sequences_Hamiltonian_2bp_2 = [Hamiltonian_double_2; Hamiltonian_double_C2_2_2A; Hamiltonian_double_C2_2B_2; Hamiltonian_double_C2_2_2_2X; Hamiltonian_double_C2_2Y; Hamiltonian_double_C2_2A_2A; Hamiltonian_double_C2_2_2_2_2A];
sorted = sort([0.777048652445489 0.815803162187861 0.825199891201852 0.832698761282255 0.841291448057122 0.852951933232366 0.876417069648107 0.736740607519727 0.778506198652469 0.789601309899058 0.799266926400893 0.802988704367723 0.816216256906491 0.846466326508041]);
x_3_1 = linspace(1000,13000,7);
% Hamiltonian_C5_2_1 = exp(-((normalization_3_bm * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_3_bmis .*((((13-8)-1)* 0.3846))+ (((13-10)-1)* 0.2307))));
% x = transpose(x);
y_3_1 = sequences_Hamiltonian_2bp_1;
y_3_2 = sequences_Hamiltonian_2bp_2;
scatter(x_3_1,y_3_1,'g') 
hold on;
scatter(x_3_1 , y_3_2, 'b')
xticks([1000 3000 5000 7000 9000 11000 13000]);
%xticklabels({ 'AATTTA,ACAGAG,CaGTCAGAAnATCCACGTA ' , ' AATTTA,ACAGaG,CAGTCAGAAnAATCCACGTA' , ' AATTTA,AcAGaG,CAGTCAGAAnAATCCACGTA ' , ' AATTTA,ACAgaG,CAGTCAGAAnATCCACGTA ' , ' AATTTA,ACAGAG,CaGtCAGAAnTCCACGTA' , , '{\color{orange}AATTTA,ACAGAG,CaGTCAnAAAAATCCACGTA}' , '{\color{orange}AATTTA,ACAGaG,CAGTCAnAAAAATCCACGTA}' , ' AATTTA,AcAgaG,CAGtCAnAAAAATCCACGTA' , ' AATTTA,ACAgaG,CAGtCAAAnAAAATCCACGTA' , ' AATTTA,ACAGAG,CaGtCAAAAnAAATCCACGTA', 'AATTTA,AcAgAG,CAgtCAnAAAAATCCACGTA ' , ' AATTTA,ACAGAG,CaGTCAGAAnATCCACGTA ' , '  AATTTA,AcAGaG,CAGTCAGAAnAATCCACGTA ' , '{\color{orange}AATTTA,ACAGAG,CaGtCAGAAnAATCCACGTA}'});
xticklabels({ '{\color{orange}Sequence 1}', ' Sequence 2', 'Sequence 3' , ' Sequence 4 '  , '{\color{orange}Sequence 5}' , ' Sequence 6 ' , 'Sequence 7'})
xtickangle(25)
yticks(sorted)
title('CASE 1: {\color{blue}Two base pair mismatches}')
xlabel('TARGET SEQUENCE')
ylabel('HAMILTONIAN')
legend({'Computed Hamiltonian for Binding Position 13' 'Computed Hamiltonian for Binding Position 15'}, 'Location', 'southoutside')
% hold on;
% plot(x_3_1 , y_3_2, '-o', '*', [0 1 0])
hold on;
refline
% uitable('Data', [AATTTAAcAGaGCAGTCAGAAnAATCCACGTA  AATTTAACAgaGCAGTCAGAAnATCCACGTA], 'ColumnName', {'Orange sequence form'}, 'Position', [20 20 500 150]);

% subplot(2,1,2)

figure(2)
% subplot(2,1,1)
% define sequences_Hamiltonian as in demo .m

sequences_Hamiltonian_3bp_1 = [Hamiltonian_triple; Hamiltonian_triple_2_A; Hamiltonian_triple_A; Hamiltonian_triple_A_2; Hamiltonian_triple_II; Hamiltonian_TRIPLE_II];
sequences_Hamiltonian_3bp_2 = [Hamiltonian_triple_2_B; Hamiltonian_triple_2_2A; Hamiltonian_triple_A_2; Hamiltonian_triple_A_2_2; Hamiltonian_TRIPLE_II_2; Hamiltonian_triplee_II_2];

sorted = sort([0.777048652445489 0.815803162187861 0.825199891201852 0.832698761282255 0.841291448057122 0.852951933232366 0.876417069648107 0.736740607519727 0.778506198652469 0.789601309899058 0.799266926400893 0.802988704367723 0.816216256906491 0.846466326508041]);
x_3_11 = linspace(500,12000,6);
% Hamiltonian_C5_2_1 = exp(-((normalization_3_bm * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_3_bmis .*((((13-8)-1)* 0.3846))+ (((13-10)-1)* 0.2307))));
% x = transpose(x);
y_3_11 = sequences_Hamiltonian_3bp_1;
y_3_22 = sequences_Hamiltonian_3bp_2;
scatter(x_3_11,y_3_11,'g') 
hold on;
scatter(x_3_11 , y_3_22, 'b')
xticks([500 3000 5000 7000 9000 12000]);
%xticklabels({ 'AATTTA,ACAGAG,CaGTCAGAAnATCCACGTA ' , ' AATTTA,ACAGaG,CAGTCAGAAnAATCCACGTA' , ' AATTTA,AcAGaG,CAGTCAGAAnAATCCACGTA ' , ' AATTTA,ACAgaG,CAGTCAGAAnATCCACGTA ' , ' AATTTA,ACAGAG,CaGtCAGAAnTCCACGTA' , , '{\color{orange}AATTTA,ACAGAG,CaGTCAnAAAAATCCACGTA}' , '{\color{orange}AATTTA,ACAGaG,CAGTCAnAAAAATCCACGTA}' , ' AATTTA,AcAgaG,CAGtCAnAAAAATCCACGTA' , ' AATTTA,ACAgaG,CAGtCAAAnAAAATCCACGTA' , ' AATTTA,ACAGAG,CaGtCAAAAnAAATCCACGTA', 'AATTTA,AcAgAG,CAgtCAnAAAAATCCACGTA ' , ' AATTTA,ACAGAG,CaGTCAGAAnATCCACGTA ' , '  AATTTA,AcAGaG,CAGTCAGAAnAATCCACGTA ' , '{\color{orange}AATTTA,ACAGAG,CaGtCAGAAnAATCCACGTA}'});
xticklabels({ '{\color{orange}Sequence 1}', ' Sequence 2', 'Sequence 3' , ' Sequence 4 '  , '{\color{orange}Sequence 5}' , ' Sequence 6 ' })
xtickangle(25)
yticks(sorted)
title('CASE 2: {\color{blue}Three base pair mismatches}')
xlabel('TARGET SEQUENCE')
ylabel('HAMILTONIAN')
legend({'Computed Hamiltonian for Binding Position 13' 'Computed Hamiltonian for Binding Position 15'}, 'Location', 'southoutside')
% hold on;
% plot(x_3_1 , y_3_2, '-o', '*', [0 1 0])
hold on;
refline

visit_positions_1 = zeros(1,20);
visit_positions_2 = zeros(1,20);
visit_positions_3 = zeros(1,20);
visit_positions_4 = zeros(1,20);
visit_positions_5 = zeros(1,20);
visit_positions_6 = zeros(1,20);
visit_positions_7 = zeros(1,20);
visit_positions_8 = zeros(1,20);
visit_positions_9 = zeros(1,20);
visit_positions_10 = zeros(1,20);
visit_positions_11 = zeros(1,20);
visit_positions_12 = zeros(1,20);
visit_positions_13 = zeros(1,20);

visit_positions_11 = zeros(1,22);
visit_positions_22 = zeros(1,22);
visit_positions_33 = zeros(1,22);
visit_positions_44 = zeros(1,22);
visit_positions_55 = zeros(1,22);
visit_positions_66 = zeros(1,22);
visit_positions_77 = zeros(1,22);
visit_positions_88 = zeros(1,22);
visit_positions_99 = zeros(1,22);
visit_positions_100 = zeros(1,22);
visit_positions_111 = zeros(1,22);
visit_positions_122 = zeros(1,22);
visit_positions_133 = zeros(1,22);

visit_positions_1(1)=1;
visit_positions_2(1)=1;
visit_positions_3(1)=1;
visit_positions_4(1)=1;
visit_positions_5(1)=1;
visit_positions_6(1)=1;
visit_positions_7(1)=1;
visit_positions_8(1)=1;
visit_positions_9(1)=1;
visit_positions_10(1)=1;
visit_positions_11(1)=1;
visit_positions_12(1)=1;

% visit_array_1 = {visit_positions_1; visit_positions_2; visit_positions_3; visit_positions_4; visit_positions_5; visit_positions_6; visit_positions_7; visit_positions_8; visit_positions_9; visit_positions_10; visit_positions_11; visit_positions_12};

visit_array_1 = {visit_positions_1 ; visit_positions_2; visit_positions_3; visit_positions_4};

visit_array_2 = {visit_positions_1; visit_positions_2; visit_positions_3; visit_positions_4; visit_positions_5; visit_positions_6; visit_positions_7; visit_positions_8; visit_positions_9; visit_positions_10; visit_positions_11; visit_positions_12};
visit_array_3 = {visit_positions_1; visit_positions_2; visit_positions_3; visit_positions_4; visit_positions_5; visit_positions_6; visit_positions_7; visit_positions_8; visit_positions_9; visit_positions_10; visit_positions_11; visit_positions_12};

% visit_array_3_2 = {visit_positions_11; visit_positions_22; visit_positions_33; visit_positions_44; visit_positions_55; visit_positions_66; visit_positions_77; visit_positions_88; visit_positions_99; visit_positions_100; visit_positions_111; visit_positions_122; visit_positions_133};

visit_array_3_2 = {visit_positions_11; visit_positions_22; visit_positions_33};

visit_array_3_single = {visit_positions_11};

L = size(visit_array_1);
L_2 = size(visit_array_3_2);
% GENERALIZE the code by having arbitrary
% sequences, as shown below, of transition
% probabilities
TEST = linspace(exp(-1/2), exp(-3/2),20);
% try out the random walk for TEST above

for m = 1 : length(sequences_Hamiltonian_3bp_1) + length(sequences_Hamiltonian_3bp_2)
 for J = 1 : L(1)
  for I = 2 : length(TEST)
    for n = 1 : 250
    if rand > TEST(I)
      % Change only x.
      % x_t(n+1) = x_t(n) + 1;
      visit_array_1{J}(I)=visit_array_1{J}(I)+1;
    else
      % Change only y.
      % x_t(n+1) = x_t(n)-1;
      visit_array_1{J}(I-1)=visit_array_1{J}(I-1)+1;
    end
    % distance = sqrt(x_t(n+1)^2 + y_t(n+1)^2);
    end
  end 
 end 
end

TEST_2 = linspace(exp(-1/2), exp(-10/2),20);

for m = 1 : length(sequences_Hamiltonian_3bp_1) + length(sequences_Hamiltonian_3bp_2)
 for J = 1 : L(1)
  for I = 2 : length(TEST)
    for n = 1 : 250
    if rand > TEST_2(I)
      % Change only x.
      % x_t(n+1) = x_t(n) + 1;
      visit_array_2{J}(I)=visit_array_2{J}(I)+1;
    else
      % Change only y.
      % x_t(n+1) = x_t(n)-1;
      visit_array_2{J}(I-1)=visit_array_2{J}(I-1)+1;
    end
    % distance = sqrt(x_t(n+1)^2 + y_t(n+1)^2);
    end
  end 
 end 
end

TEST_3 = [exp(-10/2) ; exp(-8/2) ; exp(-7/2); exp(-6/2); exp(-5/2); exp(-4/2); exp(-3/2); exp(-2/2);  exp(-1); exp(-1)/2; exp(-1)/4; exp(-1)/6; 0.00000000009; exp(-1)/7; exp(-1)/8; exp(-1)/9; exp(-1)/10; exp(-1)/11; exp(-1)/12; exp(-1)/12]; 

% TEST_3_2 = [0;exp(-10/2) ; exp(-8/2) ; exp(-7/2); exp(-6/2); exp(-5/2); exp(-4/2); exp(-3/2); exp(-2/2);  exp(-1); exp(-1)/2; exp(-1)/4; exp(-1)/6; 0;0;0;0;0;0;0]; 
format long
TEST_3_2 = [exp(-1)/3; exp(-1)/4 ; exp(-1)/5 ; exp(-1)/6 ; exp(-1)/7; exp(-1)/8; exp(-1)/9 ; exp(-1)/10; exp(-1)/11 ; exp(-1)/12 ; exp(-1)/13 ; exp(-1)/14 ; exp(-1)/15; exp(-1)/100; exp(-1)/105 ; exp(-1)/110; exp(-1)/120; exp(-1)/125 ; exp(-1)/130; exp(-1)/135; exp(-1)/140 ; 0.000001];

% initialize the TEST_3_NEW array

TEST_3_NEW = zeros(1,20);
TEST_3_NEW(1) = exp(-1/2);
N_bind = 15;
N_mis = 5;
syms sym1 
syms sym2
func_triple_A_11_test_3 = @(N,N_1, X_2, mis, lambda_mismatch_1) exp(-(symprod(1./((N- mis) + 2)^3 *  abs(N-sym1), sym1 , 1 , (N-mis)+1) * symprod(2./(mis + 2)^3 *  (1-abs(N-sym2)), sym2 , 1, mis+1)))./(1+(lambda_c * exp(-beta * Table_4{JI,8}))  + (N_1 * lambda_mismatch_1 * exp(- beta * lambda_mismatch_1 * (-X_2))));

% test func_triple values, for fixed N,N_1
func_value_array = [ ];
format long
figure(3)
subplot(2,1,1)
for A = 1 : 20
    func_value_array(A) = double(func_triple_A_11_test_3(N_bind,(25-A)^3, A , 1 , exp(-(20-A)/2)));
    plot(A,func_value_array(A), ' . ');
    hold on;
    title('{\color{blue}Single base pair mismatch, N_{mis}=25}')
    xlabel('TARGET SEQUENCE POSITION')
    ylabel('APPROXIMATED HAMILTONIAN')
end 

func_value_array_2 = [ ];
format long
subplot(2,1,2)
for A = 1 : 20
    func_value_array(A) = double(func_triple_A_11_test_3(N_bind,(20-A)^3, A , 1 , exp(-(20-A)/2)));
    plot(A,func_value_array(A), ' . ');
    hold on;
    title('{\color{blue}Single base pair mismatch, N_{mis}=20}')
    xlabel('TARGET SEQUENCE POSITION')
    ylabel('APPROXIMATED HAMILTONIAN')
end 

% single base pair mismatch, formula 1
for IJ = 2 : length(TEST_3_NEW)
    if isequal(IJ,N_mis)
        TEST_3_NEW(IJ) = (1 - double(func_triple_A_11_test_3(N_bind, (20 - IJ)^3 ,IJ, 1, exp(-(N_bind-IJ)/2)))) * TEST_3_NEW(IJ-1);
    else
        good_t_1 = double(func_triple_A_11_test_3(N_bind, (20 - IJ)^3 , IJ , 1, exp(-(N_bind-IJ))/2));
        good_t_2 = double(func_triple_A_11_test_3(15, (20 - 2)^3 , 2 , 1, exp(-(15-2))/2));
        TEST_3_NEW(IJ) = abs((1 - (good_t_2+good_t_1))) * TEST_3_NEW(IJ-1);
    end
end

%single base pair mismatch, formula 2
for IJ = 2 : length(TEST_3_NEW)
    if isequal(IJ,N_mis)
        TEST_3_NEW(IJ) = (1 - double(func_triple_A_11_test_3(N_bind, (20 - IJ)^3 ,IJ, 1, exp(-(N_bind-IJ)/2)))) * TEST_3_NEW(IJ-1);
    else
        good_t_1 = double(func_triple_A_11_test_3(N_bind, (20 - IJ)^3 , IJ , 1, exp(-(N_bind-IJ))/2));
        good_t_2 = double(func_triple_A_11_test_3(15, (20 - 2)^3 , 2 , 1, exp(-(15-2))/2));
        TEST_3_NEW(IJ) = abs((1 - (good_t_2+good_t_1))) * TEST_3_NEW(IJ-1);
    end
end 

% double base pair mismatch

visit_array_11 = zeros(1,13);
TESTY = [ 0.998080624523906; 0.996168604406074 ;  0.994263897446976 ;0.992366461769276; 0.743529869249711 ; 0.988593238341331 ; 0.986717368419950 ; 0.984848605431730 ; 0.982986909064973 ; 0.981132239312280 ;  0.979284556467681 ; 0.977443821123807 ; 0.975609994169083 ];
% for m = 1 : length(sequences_Hamiltonian_3bp_1) + length(sequences_Hamiltonian_3bp_2)
 % for J = 1 : L(1)
  for I = 2 : length(TESTY)
    for n = 1 : 250
    if rand >= TESTY(I)
      % Change only x.
      % x_t(n+1) = x_t(n) + 1;
      visit_array_11(I)=visit_array_11(I)+1;
    else
      % Change only y.
      % x_t(n+1) = x_t(n)-1;
      visit_array_11(I-1)=visit_array_11(I-1)+1;
    end
    % distance = sqrt(x_t(n+1)^2 + y_t(n+1)^2);
    end
  end 
 % end 
% end 

visit_array_111 = zeros(1,13);
TESTY_2 = [ 0.998080624523906; 0.996168604406074 ;  0.994263897446976 ;0.992366461769276; 0.743529869249711 ; 0.751047119144556 ; 0.749963941025334 ; 0.748883883607377 ; 0.747806933423691 ; 0.746733077084662 ;  0.745662301277492 ; 0.744594592765657 ];

% for m = 1 : length(sequences_Hamiltonian_3bp_1) + length(sequences_Hamiltonian_3bp_2)
 % for J = 1 : L(1)
  for I = 2 : length(TESTY_2)
    for n = 1 : 250
    if rand >= TESTY_2(I)
      % Change only x.
      % x_t(n+1) = x_t(n) + 1;
      visit_array_111(I)=visit_array_111(I)+1;
    else
      % Change only y.
      % x_t(n+1) = x_t(n)-1;
      visit_array_111(I-1)=visit_array_111(I-1)+1;
    end
    % distance = sqrt(x_t(n+1)^2 + y_t(n+1)^2);
    end
  end 

visit_array_1111 = zeros(1,13);
TESTY_3 = [ 0.440713480549430; 0.996168604406074 ; 0.994263897446976 ; 0.992366461769276; 0.990476255814766 ; 0.988593238341331 ; 0.986717368419950 ; 0.984848605431730 ; 0.982986909064973 ; 0.981132239312280 ;  0.979284556467681 ; 0.977443821123807 ; 0.975609994169083 ];
% for m = 1 : length(sequences_Hamiltonian_3bp_1) + length(sequences_Hamiltonian_3bp_2)
 % for J = 1 : L(1)
  for I = 2 : length(TESTY_3)
    for n = 1 : 250
    if rand >= TESTY_3(I)
      % Change only x.
      % x_t(n+1) = x_t(n) + 1;
      visit_array_1111(I)=visit_array_1111(I)+1;
    else
      % Change only y.
      % x_t(n+1) = x_t(n)-1;
      visit_array_1111(I-1)=visit_array_1111(I-1)+1;
    end
    % distance = sqrt(x_t(n+1)^2 + y_t(n+1)^2);
    end
  end 
  
% try to implement different boundary conditions 
% past the 1st and 20th base pairs

% define initial transition probability

% initial_t_prob = exp(-1);  

% for m = 1 : length(sequences_Hamiltonian_3bp_1) + length(sequences_Hamiltonian_3bp_2)
    
 
     
%  for J = 1 : L_2(1)
%   for I = 2 : length(TEST_3_2)
%     for n = 1 : 250
%         X = rand;
%     if X >= TEST_3_2(I)
%       % Change only x.
%       % x_t(n+1) = x_t(n) + 1;
%       visit_array_3_single{J}(I)=visit_array_3_single{J}(I)+1;
%     else
%         
%         
%         
%       % Change only y.
%       % x_t(n+1) = x_t(n)-1;
% %       if isequal(I,2)
% %          visit_array_3_2{J}{I}=visit_array_3_2{J}(I-1)+1;
% %       else 
% %         visit_array_3_2{J}(I-1)=visit_array_3_2{J}(I-1)+1;
% %       end 
%     end
%     % distance = sqrt(x_t(n+1)^2 + y_t(n+1)^2);
%     end
%   end 
%  end 
% end

% Try out different normalizations to the Hamiltonian
% so that a sharper decrease in transition probabilities occurs

% simulate TENTH case of 3 base pair mismatch,
% in which we have ONE base pair mismatch
% within the first 6 position, while we
% also have 2 additional mismatches past
% the first 6 positions, for a 
% position of binding, N = 13

normalization_double_1_1 = 1/(9+2)^3;
normalization_double_2_1 = 1/(2+3)^3;

% construct the Hamiltonian. The T_ij factor for
% the FIRST mismatch, at 4,comes from the fact that
% |4-13|/13 = 0.692, while the T_ij factor for the
% SECOND mismatch, at 10, comes from the fact that
% |10-13|/13 = 0.231, while finally, the THIRD
% T_ij factor comes from the fact that 
% 0.15 = |11-13|/13 

Hamiltonian_triple_A_1 = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-12))))+(normalization_double_2_1 .*((((13-4)-1)* 0.6923)+ (((13-10)-1)* 0.231)+(((13-11)-1)* 0.15)))));

% generate a test function with specific normalizations

func_triple_A_1 = @(X)(Hamiltonian_triple_A_1)./(1+(exp(X)*(lambda_c .* exp(-beta * Table{I,8}) + lambda_c * exp(-beta * Table_4{JI,8}))) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X))));

TEST_func = linspace(1,20,20);

figure(5)
for II = 1 : length(TEST_func)
    X_i = func_triple_A_1(II);
    plot(TEST_func(II),X_i, '.')
    hold on;
end 

func_triple_A_1_2 = @(X)(Hamiltonian_triple_A_1)./(1+(exp(X)*(lambda_c .* exp(-beta * Table{I,8}) + lambda_c * exp(-beta * Table_4{JI,8}))) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (17-X))));
figure(6)
for II = 1 : length(TEST_func)
    X_i = func_triple_A_1_2(II);
    plot(TEST_func(II),X_i, '.')
    hold on;
end 

partition_func = @(X) 1+(exp(-X)*(lambda_c .* exp(-beta * Table{I,8}) + lambda_c * exp(-beta * Table_4{JI,8}))) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (13-X)));
p_func_array = [ ];

for K = 1 : 12
    X_i = partition_func(K);
    p_func_array(K)=X_i;
end 

Hamiltonian_array = [ ];
Hamiltonian_array(1) = exp(-((normalization_double_1_1 * (((13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2_1 .*(13-1)-1)));
Hamiltonian_array(2) = exp(-((normalization_double_1_1 * (((13-1)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2_1 .*(13-2)-1)));
Hamiltonian_array(3) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-4)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2_1 .*(13-3)-1)));
Hamiltonian_array(4) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2_1 .*(13-4)-1)));
Hamiltonian_array(5) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2_1 .*(13-5)-1)));
Hamiltonian_array(6) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2_1 .*(13-6)-1)));
Hamiltonian_array(7) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-8)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2_1 .*(13-7)-1)));
Hamiltonian_array(8) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-9)+(13-10)+(13-11)+(13-12))))+(normalization_double_2_1 .*(13-8)-1)));
Hamiltonian_array(9) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-10)+(13-11)+(13-12))))+(normalization_double_2_1 .*(13-9)-1)));
Hamiltonian_array(10) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-11)+(13-12))))+(normalization_double_2_1 .*(13-10)-1)));
Hamiltonian_array(11) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-12))))+(normalization_double_2_1 .*(13-11)-1)));
Hamiltonian_array(12) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-11))))+(normalization_double_2_1 .*(13-12)-1)));
% Hamiltonian_array(13) = exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-10)+(13-12))))+(normalization_double_2_1 .*(13-1)-1)));
% Hamiltonian_array(14) =
% Hamiltonian_array(15) =
% Hamiltonian_array(16) =
% Hamiltonian_array(17) =
% Hamiltonian_array(18) =
% Hamiltonian_array(19) =
% Hamiltonian_array(20) =

% energy_array = [ ];
% 
% for A = 1:12
%     energy_array(A) = Hamiltonian_array(A)./p_func_array(A);
% end 
% 
% figure(5)
% plot(energy_array,linspace(1,12,12))

% define match and mismatch vector
% then define each of the functions


% func_triple_A_1_3 = @(X)(Hamiltonian_triple_A_1)./(1+(exp(-X)*(lambda_c .* exp(-beta * Table{I,8}) + lambda_c * exp(-beta * Table_4{JI,8}))) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (13-X))));
% figure(5)
% for II = 1 : length(TEST_func)
%     X_i = func_triple_A_1_3(II);
%     plot(TEST_func(II),X_i, '.')
%     hold on;
% end 

% plot func for different binding positions
% we will take possible binding positions
% from 6 to 17

% figure(6)
% for K = 1 : 5
%     subplot(5,1,K)
% for JJ =  12 : 17
%     Hamiltonian_triple_A_1 = exp(-((normalization_double_1_1 * (((JJ-1)+(JJ-2)+(JJ-3)+(JJ-5)+(JJ-6)+(JJ-7)+(JJ-8)+(JJ-9)+(JJ-12))))+(normalization_double_2_1 .*((((JJ-4)-1)* 0.6923)+ (((JJ-10)-1)* 0.231)+(((JJ-11)-1)* 0.15)))));
%     func_triple_A_11 = @(X)(Hamiltonian_triple_A_1)./(1+(exp(-X)*(lambda_c .* exp(-beta * Table{I,8}) + lambda_c * exp(-beta * Table_4{JI,8}))) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (JJ-X))));
%     % binding_position_vec = linspace(JJ,20, 20-JJ+1);
%     for JK = 1 : 20
%         X_i = func_triple_A_11(JK);
%         plot(TEST_func(JK),X_i, '.')
%         hold on;
%     end
% end 
% end 

% prepare each Hamiltonian individually

% Hamiltonian_array = [ ];
% 
% normalization_double_1_1 = 1/(9+2)^3;
% normalization_double_2_1 = 1/(2+3)^3;
% 
% for JJ = 12 : 17
%      Hamiltonian_triple_A_1 = exp(-((normalization_double_1_1 * (((JJ-1)+(JJ-2)+(JJ-3)+(JJ-5)+(JJ-6)+(JJ-7)+(JJ-8)+(JJ-9)+(JJ-12))))+(normalization_double_2_1 .*((((JJ-4)-1)* 0.6923)+ (((JJ-10)-1)* 0.231)+(((JJ-11)-1)* 0.15)))));
%      Hamiltonian_array(JJ-11) = Hamiltonian_triple_A_1;
% end 

% prepare the functions to be evaluated


% figure(6)
% for K = 1 : 5
%     subplot(5,1,K)
% for JJ = 12:17
% func_triple_A_11 = @(X) exp(-((normalization_double_1_1 * (((JJ-1)+(JJ-2)+(JJ-3)+(JJ-5)+(JJ-6)+(JJ-7)+(JJ-8)+(JJ-9)+(JJ-12))))+(normalization_double_2_1 .*((((JJ-4)-1)* (abs(JJ-4))./JJ)+ (((JJ-10)-1)* (abs(JJ-10))./JJ)+(((JJ-11)-1)* (abs(JJ-11))./JJ))))./(1+(exp(X)*(lambda_c .* exp(-beta * Table{I,8}) + lambda_c * exp(-beta * Table_4{JI,8}))) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (JJ-X)))));
% for I = 1:20
%         X_i = func_triple_A_11(I);
%         plot(TEST_func(I),X_i, '.')
%         % hold on;
% end 
% end 
% end 

% binding position 12
% figure(6)
% for I = 1:20
% func_triple_A_11 = @(X) exp(-((normalization_double_1_1 * (((12-1)+(12-2)+(12-3)+(12-5)+(12-6)+(12-7)+(12-8)+(12-9))))+(normalization_double_2_1 .*((((12-4)-1)* (abs(12-4))./12)+ (((12-10)-1)* (abs(12-10))./12)+(((JJ-11)-1)* (abs(12-11))./12))))./(1+(exp(X)*(lambda_c .* exp(-beta * Table{I,8}) + lambda_c * exp(-beta * Table_4{JI,8}))) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (12-X)))));
% % plot
% end
% 
% % func_triple_A_11 = @(N,X) exp(-((normalization_double_1_1 * (((12-1)+(12-2)+(12-3)+(12-5)+(12-6)+(12-7)+(12-8)+(12-9))))+(normalization_double_2_1 .*((((12-4)-1)* (abs(12-4))./12)+ (((12-10)-1)* (abs(12-10))./12)+(((JJ-11)-1)* (abs(12-11))./12))))./(1+(exp(X)*(lambda_c .* exp(-beta * Table{I,8}) + lambda_c * exp(-beta * Table_4{JI,8}))) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (N-X)))));
% 
% 
% % binding position 13
% figure(7)
% 
% func_triple_A_11 = @(X) exp(-((normalization_double_1_1 * (((13-1)+(13-2)+(13-3)+(13-5)+(13-6)+(13-7)+(13-8)+(13-9)+(13-12))))+(normalization_double_2_1 .*((((13-4)-1)* (abs(13-4))./13)+ (((13-10)-1)* (abs(13-10))./13)+(((13-11)-1)* (abs(13-11))./13))))./(1+(exp(X)*(lambda_c .* exp(-beta * Table{I,8}) + lambda_c * exp(-beta * Table_4{JI,8}))) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (13-X)))));
% 
% 
% % binding position 14
% figure(8)
% 
% func_triple_A_11 = @(X) exp(-((normalization_double_1_1 * (((JJ-1)+(JJ-2)+(JJ-3)+(JJ-5)+(JJ-6)+(JJ-7)+(JJ-8)+(JJ-9)+(JJ-12))))+(normalization_double_2_1 .*((((JJ-4)-1)* (abs(JJ-4))./JJ)+ (((JJ-10)-1)* (abs(JJ-10))./JJ)+(((JJ-11)-1)* (abs(JJ-11))./JJ))))./(1+(exp(X)*(lambda_c .* exp(-beta * Table{I,8}) + lambda_c * exp(-beta * Table_4{JI,8}))) + (lambda_mismatch .* exp(- beta * lambda_mismatch .* (JJ-X)))));
% 


% VISIT_AVG = [ ];
% compute average, the mean field type approximation
% for A = 1 : L(1)
%     temp = [];
%     for AI = 1 : L(1)
%        temp(:,AI) = visit_array{AI}; 
%        VISIT_AVG(A) = sum(temp(:,AI))/length(temp(:,AI));
%     end 
% end

% figure(3)
% Y = linspace(0,22,22);
% for J = 1 : L(1)
%     scatter(visit_array_1{J},Y)
%     hold on;
% end

% figure(4)
% Y_1 = linspace(1,20,20);
% for J = 1 : L(1)
%     scatter(visit_array_2{J},Y_1)
%     hold on;
% end 

% plot(visit_positions_2, [1,20], 'LineWidth', 2);
%   caption = sprintf('Random Walk #%d of %d', m,  length(sequences_Hamiltonian_3bp_1) + length(sequences_Hamiltonian_3bp_2));
%   title(caption, 'FontSize', 25, 'FontWeight', 'bold');
%   hold on
%   drawnow;

end  