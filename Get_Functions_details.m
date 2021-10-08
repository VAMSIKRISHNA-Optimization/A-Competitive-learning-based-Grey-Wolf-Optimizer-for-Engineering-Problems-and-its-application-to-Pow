% By VAMSI KRISHNA  for Clb-GWO
% STANDARD BENCHMARK FUNCTIONS ONLY
% 11 UNIMODAL
% 13 MULTIMODAL 

% This function containts full information and implementations of the benchmark 
% functions in Table 1, Table 2, and Table 3 in the paper

% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)

function [lb,ub,dim,fobj] = Get_Functions_details(F)


switch F



%%% UNIMODAL
    case 'F1'
%sphere
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F2'
%step
        fobj = @F2;
        lb=-100;
        ub=100;
          dim=100;
        
    case 'F3'
%quartic
        fobj = @F3;
        lb=-1.28;
        ub=1.28;
           dim=100;
        
    case 'F4'
%sum squares
        fobj = @F4;
        lb=-10;
        ub=10;
           dim=100;
        
    case 'F5'
%Sum of different powers function
        fobj = @F5;
        lb=-1;
        ub=1;
           dim=100;
        
    case 'F6'
% ROTATED HYPER-ELLIPSOID FUNCTION
        fobj = @F6;
        lb=-100;
        ub=100;
           dim=100;
           
    case 'F7'
%Schwefel 2.21 Function
        fobj = @F8;
        lb=-100;
        ub=100;
           dim=100;
        
    case 'F8'
% Dixon-price
        fobj = @F8;
        lb=-10;
        ub=10;
           dim=100;
       
   case 'F9'
%Exponential
        fobj = @F9;
        lb=-1;
        ub=1;
           dim=100;
        
           
    case 'F10'
%POWELL
        fobj = @F10;
        lb=-4;
        ub=5;
        dim=100;
        
    case 'F11'
%BROWN FUNCTION
        fobj = @F11;
        lb=-1;
        ub=4;
        dim=100;  
      
  


%%% MULTI-MODAL
    case 'F12'
%Schwefel 1.2
        fobj = @F12;
        lb=-100;
        ub=100;
        dim=100;
        
    case 'F13'
% Periodic Function
        fobj = @F13;
        lb=-10;
        ub=10;
        dim=100;
        
    case 'F14'
	%ROSENBROCK
       fobj = @F14;
        lb=-5;
        ub=10;
        dim=100;    
        
    case 'F15'
%Rastrigin
        fobj = @F15;
        lb=-5.12;
        ub=5.12;
        dim=100;
        
    case 'F16'
%Ackley
        fobj = @F16;
        lb=-32;
        ub=32;
        dim=100;

    case 'F17'
%Griewank
        fobj = @F17;
        lb=-600;
        ub=600;
        dim=100;

    case 'F18'
%Generalized penalized function 1
        fobj = @F18;
        lb=-50;
        ub=50;
        dim=100;
        
    case 'F19'
%Generalized penalized function 2
        fobj = @F19;
        lb=-50;
        ub=50;
        dim=100;
        
   case 'F20'
        % Salomon
	fobj = @F20;
        lb=-100;
        ub=100;
        dim=100;
             
    case 'F21'
     %Alpine N.1
        fobj = @F21;
        lb=0;
        ub=10;
        dim=100;     

 case 'F22'
     %Qing
        fobj = @F22;
        lb=-500;
        ub=500;
        dim=100;              

  case 'F23'
     %Happy Cat Function
        fobj = @F23;
        lb=-2;
        ub=2;
        dim=100;       

case 'F24'
%Ackley N.4 function
        fobj = @F24;
        lb=-35;
        ub=35;
        dim=100; 
end

end





%%%%%%% UNI-MODAL %%%%%%%%%%%
% F1
function o = F1(x)
%SPHERE
o=sum(x.^2);
end


% F2
function o = F2(x)
%STEP
o = sum(floor(x+.5).^2, 2);
end


% F3
function o = F3(x)
%Quartic Function
 n = size(x, 2);
    
    scores = 0;
    for i = 1:n
        scores = scores + i *(x(:, i) .^ 4);
    end
     
    o  = scores + rand;
end


% F4
%Sum Squares Function
function o = F4(x)
[m, n] = size(x);
   x2 = x .^2;
   I = repmat(1:n, m, 1);
   o = sum( I .* x2, 2);
end


% F5
% SUM OF DIFFERENT POWERS FUNCTION
function o = F5(x)
d = length(x);
sum = 0;

for i = 1:d
    xi = x(i);
    new = (abs(xi))^(i+1);
    sum = sum + new;
end

o = sum;
end


% F6
% ROTATED HYPER-ELLIPSOID FUNCTION
function o = F6(x)
 d= length(x);
outer = 0;
for ii = 1:d
    inner = 0;
    for jj = 1:ii
        xj = x(jj);
        inner = inner + xj^2;
    end
    o= outer + inner;
end
y = outer;
end


% F7
function o = F7(x)
%Schwefel 2.21 Function
o = max(abs(x), [], 2);
end


% F8
% DIXON PRICE
function o = F8(x)
x1 = x(1);
d = length(x);
term1 = (x1-1)^2;
sum = 0;
for ii = 2:d
	xi = x(ii);
	xold = x(ii-1);
	new = ii * (2*xi^2 - xold)^2;
	sum = sum + new;
end
o = term1 + sum;
end
    
% F9
%Exponential
function o = F9(x)
  x2 = x .^2;
   
   o = -exp(-0.5 * sum(x2, 2));
end    


% F10
function o = F10(x)
%POWELL
d = length(x);
sum = 0;

for i = 1:(d/4)
	term1 = (x(4*i-3) + 10*x(4*i-2))^2;
	term2 = 5 * (x(4*i-1) - x(4*i))^2;
	term3 = (x(4*i-2) - 2*x(4*i-1))^4;
	term4 = 10 * (x(4*i-3) - x(4*i))^4;
	sum = sum + term1 + term2 + term3 + term4;
end
o = sum;
end

% F11
%Brown Function
function o = F11(x)
n = size(x, 2);  
    scores = 0;
    
    x = x .^ 2;
    for i = 1:(n-1)
        o = scores + x(:, i) .^ (x(:, i+1) + 1) + x(:, i+1).^(x(:, i) + 1);
    end
end


%%%%%%% MULTI-MODAL FUNCTIONS %%%%%%%%%%%


% F12
% Schwefel Function 1.2
function o = F12(x)
dim = size(x,2);
o = 0;
for i=1:dim
    o = o + sum(x(:, 1:i), 2).^2;
end
end


%  F13
% Periodic Function
function o = F13(x)
sin2x = sin(x) .^ 2;
    sumx2 = sum(x .^2, 2);
    o  = 1 + sum(sin2x, 2) -0.1 * exp(-sumx2);
end

%Rosenbrock
function o = F14(x)
scores = 0;
    n = size(x, 2);
    assert(n >= 1, 'Given input X cannot be empty');
    a = 1;
    b = 100;
    for i = 1 : (n-1)
        o = scores + (b * ((x(:, i+1) - (x(:, i).^2)) .^ 2)) + ((a - x(:, i)) .^ 2);
    end
end


% F15
%Rastrigin
function o = F15(x)
n = size(x, 2);
    A = 10;
    o = (A * n) + (sum(x .^2 - A * cos(2 * pi * x), 2));
end


% F16
%Ackley Function
function o = F16(x)
n = size(x, 2);
    ninverse = 1 / n;
    sum1 = sum(x .^ 2, 2);
    sum2 = sum(cos(2 * pi * x), 2);
    
     o = 20 + exp(1) - (20 * exp(-0.2 * sqrt( ninverse * sum1))) - exp( ninverse * sum2);
end


% F17
%Griewank Function
function o = F17(x)
 n = size(x, 2);
    
    sumcomp = 0;
    prodcomp = 1;
    
    for i = 1:n
        sumcomp = sumcomp + (x(:, i) .^ 2);
        prodcomp = prodcomp .* (cos(x(:, i) / sqrt(i)));
    end
    
    o = (sumcomp / 4000) - prodcomp + 1;
end





% F18
function o = F18(x)
%Generalized penalized function 1
dim=size(x,2);
y = 1+(x+1)./4;
o = pi/dim .* (10*sin(pi.*y(:, 1)).^2 + sum((y(:, 1:end-1)-1).^2 .* (1+10.*sin(pi.*y(:, 2:end)).^2), 2) + (y(:, end)-1).^2) + ...
    sum(Ufun(x, 10, 100, 4), 2);
end


% F19
%Generalized penalized function 2
function o = F19(x)
dim = size(x,2);
o = 0.1 .* (sin(3*pi*x(:, 1)).^2 + sum((x(:, 1:dim-1)-1).^2 .* (1+sin(3.*pi.*x(:, 2:dim)).^2), 2)+...
    ((x(:, dim)-1).^2) .* (1+sin(2*pi*x(:, dim)).^2)) + sum(Ufun(x,5,100,4), 2);
end

% F20
function o = F20(x)
%Salomon Function
x2 = x .^ 2;
    sumx2 = sum(x2, 2);
    sqrtsx2 = sqrt(sumx2);
    
    o = 1 - cos(2 .* pi .* sqrtsx2) + (0.1 * sqrtsx2);
end


% F21
function o = F21(x)
%ALPINE N.1
  o = sum(abs(x .* sin(x) + 0.1 * x), 2);
end


% F22
function o = F22(x)
%QING
 n = size(x, 2);
    x2 = x .^2;
    
    scores = 0;
    for i = 1:n
        o = scores + (x2(:, i) - i) .^ 2;
    end
end

% F23
function o = F23(x)
%Happy Cat unctionn 
if nargin < 2 
        alpha = 0.5;
    end
    
    n = size(x, 2);
    x2 = sum(x .* x, 2);
    o = ((x2 - n).^2).^(alpha) + (0.5*x2 + sum(x,2))/n + 0.5;

end 


% F24
function o = F24(x)
%Ackley N. 4 Function
[m, n] = size(x);
    
    scores = zeros(m, 1); 
   
   for i = 1:m
      for j = 1:(n - 1)
            o = scores + exp(-0.2) .* sqrt( x(i, j) .^ 2 + x(i, j + 1) .^ 2 ) ...
            + 3 * ( cos(2 * x(i, j)) + sin(2 * x(i, j + 1)) );
      end
   end
end



function o=Ufun(x,a,k,m)
o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end
