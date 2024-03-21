clc
clear all
close all

% Optimization Problem
% Simplex - Nelder Mead
% Rosenbrok Function - 3 variables - initial point [1 -2 5]
% 6/4/2023
% By Amin Pirmohammadi

% Define Function
F = @(x) (1-x(1))^2 + (2-x(2))^2 + 65*(x(2)-x(1)^2)^2 + 65*(x(3)-x(2)^2)^2;

% Initial Values
alpha_r = .7;
alpha_e = 3;
alpha_q = -0.5;
delta_j = 1;
delta_sh = 0.7;
epsilon = 1e-6;

x_0 = [1; -2; 5];
% x_0 = [1; 2; 3];
% x_0 = [0; 0; 0];
n = length(x_0);             % Number of Dimensions

% Directions
e_i = [1; 0; 0];
e_j = [0; 1; 0];
e_k = [0; 0; 1];

x_1 = x_0 + delta_j*e_i;
x_2 = x_0 + delta_j*e_j;
x_3 = x_0 + delta_j*e_k;

ConverganceTest = 10000;   % Just to define the value of Convergance

Total_Iter = 0;
reflection_Iter = 0;
shrinkage_Iter = 0;
expansion_Iter = 0;
contraction_Iter = 0;

while ConverganceTest > epsilon
Total_Iter = Total_Iter+1;

f0 = F(x_0);
f1 = F(x_1);
f2 = F(x_2);
f3 = F(x_3);

% NumberOfFuncEvaluation = n + 1;

A = [f0 f1 f2 f3]
A_sort = sort(A,'ascend')                            % Sort F(x) values without F(x_w)
Max_f = A_sort(4);
Min_f = A_sort(1);

if Max_f == f0
    x_w = x_0;
    B = [x_1 x_2 x_3];                               % Sort x values without x_w  (B)
    disp('x_w = x_0');
elseif Max_f == f1
    x_w = x_1;
    B = [x_0 x_2 x_3];
    disp('x_w = x_1');
elseif Max_f == f2
    x_w = x_2;
    B = [x_1 x_0 x_3];
    disp('x_w = x_2');
elseif Max_f == f3
    x_w = x_3;
    B = [x_1 x_2 x_0];
    disp('x_w = x_3');
end

if Min_f == f0
    x_b = x_0;
    C = [x_1 x_2 x_3];                               % Sort x values without x_b  (B)
    disp('x_b = x_0');
elseif Min_f == f1
    x_b = x_1;
    C = [x_0 x_2 x_3];
    disp('x_b = x_1');
elseif Min_f == f2
    x_b = x_2;
    C = [x_1 x_0 x_3];
    disp('x_b = x_2');
elseif Min_f == f3
    x_b = x_3;
    C = [x_1 x_2 x_0];
    disp('x_b = x_3');
end

f1 = A_sort(1);
f2 = A_sort(2);
f3 = A_sort(3);                     % Sort f1 < f2 < f3 < f4(MAX_F)
f4 = A_sort(4);

x_c = (1/n)*([sum(B(1,:)); sum(B(2,:)); sum(B(3,:))])      % centroid point of n points

% Convergance Test
ConverganceTest = (sqrt((f1-F(x_c))^2 + (f2-F(x_c))^2 + (f3-F(x_c))^2 + (f4-F(x_c))^2))/(n+1);

% Reflection
x_r = (1 + alpha_r)*x_c - alpha_r*x_w;
f_r = F(x_r);

if f1 <= f_r && f_r < f3
    % x_r is a good reflection
    D = [B x_r];                                 % Accept and replacement
    reflection_Iter = reflection_Iter + 1;
elseif f_r >= f3
    if f_r >= f3 && f_r < f4
        % Calculate x_q by outside contraction
        x_q = (1+alpha_q)*x_r - alpha_q*x_c;
        f_q = F(x_q);
        if f_q < f3
            D = [B x_r];                         % Accept and replacement
            contraction_Iter = contraction_Iter + 1;
        else
            % Shrinkage
            C(:,1) = x_b + delta_sh*(C(:,1) - x_b);
            C(:,2) = x_b + delta_sh*(C(:,2) - x_b);
            C(:,3) = x_b + delta_sh*(C(:,3) - x_b);

            D = [x_b C(:,1) C(:,2) C(:,3)];
            shrinkage_Iter = shrinkage_Iter + 1;
        end
    elseif f_r >= f4
        % Calculate x_q by inside contraction
        x_q = (1+alpha_q)*x_c - alpha_q*x_w;
        f_q = F(x_q);
        if f_q < f3
            D = [B x_r];                         % Accept and replacement
            contraction_Iter = contraction_Iter + 1;
        else
            % Shrinkage
            C(:,1) = x_b + delta_sh*(C(:,1) - x_b);
            C(:,2) = x_b + delta_sh*(C(:,2) - x_b);
            C(:,3) = x_b + delta_sh*(C(:,3) - x_b);

            D = [x_b C(:,1) C(:,2) C(:,3)];
            shrinkage_Iter = shrinkage_Iter + 1;
        end
    end
elseif f_r < f1
    %Expansion
    x_e = (1+alpha_e)*x_c - alpha_e*x_w;
    f_e = F(x_e);
    if f_e < f_r
        % x_e will be used as the new replacement point
        D = [B x_e];
        expansion_Iter = expansion_Iter + 1;
    else
        D = [B x_r];
        reflection_Iter = reflection_Iter + 1;
    end

    each_nIter{Total_Iter} = Total_Iter;
    f_eachIter{Total_Iter} = F(D(:,1));
end

x_0 = D(:,1);
x_1 = D(:,2);
x_2 = D(:,3);
x_3 = D(:,4);

x{Total_Iter} = [x_0, x_1, x_2, x_3];


end

f0 = F(x_0);
f1 = F(x_1);
f2 = F(x_2);
f3 = F(x_3);

x_optimum = D
f_optimum = min([f0,f1,f2,f3])

% Iterations
Total_Iter
reflection_Iter
contraction_Iter
shrinkage_Iter
expansion_Iter

%% Visualize the Process

ShematicNelderMead;

for i = 1:Total_Iter

    NelderX = x{1,i}(1,:);

    NelderY = x{1,i}(2,:);

    NelderZ = x{1,i}(3,:);

    figure(1)

    grid on
    Plot_Nelder_Mead
    
    % Uncomment line below if you want to make plot axis constant.
%     axis([NelderX(1)-1 NelderX(1)+1 NelderY(1)-1 NelderY(1)+1 NelderZ(1)-1 NelderZ(1)+1]) 

 pause(.1)

end

% ---------------------------------Plot Convergance ----------------------------------------
figure(2)
AA = [each_nIter{:}];
BB = [f_eachIter{:}];
plot(AA,BB)
axis([0 max(AA) -10 max(BB)])
title('Convergance (function value in each iteration)')
xlabel('Number of Iterations')
ylabel('f value')
print(gcf,'SimplexConverge.png','-dpng','-r600');
