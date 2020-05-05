function [visit_array_1,visit_array_2, visit_array_3,visit_array_4]=dist_visits(vec1,vec2,vec3,vec4)
visit_array_1 = zeros(1,13);
visit_array_2 = zeros(1,13);
visit_array_3 = zeros(1,13);
visit_array_4 = zeros(1,13);
% for m = 1 : length(sequences_Hamiltonian_3bp_1) + length(sequences_Hamiltonian_3bp_2)
 % for J = 1 : L(1)
  for I = 2 : length(vec1)
    for n = 1 : 250
    if rand > vec1(I)
      % Change only x.
      % x_t(n+1) = x_t(n) + 1;
      visit_array_1(I)=visit_array_1(I)+1;
    else
      % Change only y.
      % x_t(n+1) = x_t(n)-1;
      visit_array_1(I-1)=visit_array_1(I-1)+1;
    end
    % distance = sqrt(x_t(n+1)^2 + y_t(n+1)^2);
    end
  end 
  
   for I = 2 : length(vec2)
    for n = 1 : 250
    if rand > vec2(I)
      % Change only x.
      % x_t(n+1) = x_t(n) + 1;
      visit_array_2(I)=visit_array_2(I)+1;
    else
      % Change only y.
      % x_t(n+1) = x_t(n)-1;
      visit_array_2(I-1)=visit_array_2(I-1)+1;
    end
    % distance = sqrt(x_t(n+1)^2 + y_t(n+1)^2);
    end
   end 
  
   for I = 2 : length(vec3)
    for n = 1 : 250
    if rand > vec3(I)
      % Change only x.
      % x_t(n+1) = x_t(n) + 1;
      visit_array_3(I)=visit_array_3(I)+1;
    else
      % Change only y.
      % x_t(n+1) = x_t(n)-1;
      visit_array_3(I-1)=visit_array_3(I-1)+1;
    end
    % distance = sqrt(x_t(n+1)^2 + y_t(n+1)^2);
    end
   end 
  
  for I = 2 : length(vec4)
    for n = 1 : 250
    if rand > vec4(I)
      % Change only x.
      % x_t(n+1) = x_t(n) + 1;
      visit_array_4(I)=visit_array_4(I)+1;
    else
      % Change only y.
      % x_t(n+1) = x_t(n)-1;
      visit_array_4(I-1)=visit_array_4(I-1)+1;
    end
    % distance = sqrt(x_t(n+1)^2 + y_t(n+1)^2);
    end
  end 
  
   
  figure(1)
  subplot(5,1,1)
  plot(linspace(1,13,13),visit_array_1,'-o')
  
  subplot(5,1,2)
  plot(linspace(1,13,13),visit_array_2,'-o')
  
  subplot(5,1,3)
  plot(linspace(1,13,13),visit_array_3,'-o')
  
  subplot(5,1,4)
  plot(linspace(1,13,13),visit_array_4,'-o')
  
  subplot(5,1,5)
  plot(linspace(1,13,13),visit_array_1,'-o')
  hold on;
  plot(linspace(1,13,13),visit_array_2,'-o')
  hold on;
  plot(linspace(1,13,13),visit_array_3,'-o')
  
end 