
% ----------------------------------------
%   Solution of least squares  problem  min_x || Ax - y ||_2
%   using the method of normal equations.
%   Matrix A is constructed as a Vandermonde matrix.
%   
% ----------------------------------------

degree=2; % graden på polynomet 
points=200; %antalet diskreta punkter = ( rader i matrisen A)

x=zeros(1,points);
y=zeros(1,points);
Vandermode_matrix=[];
for i=1:1:points
  x = linspace(-10.0,10.0,points);
  %  exakt funktion
  %y(i)= sin(pi*x(i)/5) %+ x(i)/5;
  % y=rand(1,i)*0.1;
  y(i) = 1./(1+x(i).^2) + 0.1*randn(size(x(i)));
 
end

%skapar vandermode matrisen 
Vandermode_matrix = bsxfun(@power,x(:),0:degree);


% beräknar högeledet 
rhs=Vandermode_matrix'*y';

% beräknar vänsterledet
lhs=Vandermode_matrix'*Vandermode_matrix;

l=zeros(degree+1);

% solution of the normal equation using Cholesky decomposition

if all(eig(lhs) > 0)
    l = chol(lhs, 'lower');
    rhs = l' \ (l \ rhs);

    figure(1)
plot(x,y,'- r', 'linewidth',1)
hold on

% beräknar approximation till det exakta polynomet med coefficienter c

approx = Vandermode_matrix*rhs;
plot(x,approx,'- b', 'linewidth',1)
hold off

str_xlabel = ['polynomgrad: ', num2str(degree)];

% Beräknar det relativa felet 
%  norm(approx. value - true value) / norm(true value)
%e1=norm(y'- approx)/norm(y')

error=norm(y'- approx)/norm(y');
%fprintf('Relativt fel: %.2f%%', error);

legend('exact ',str_xlabel);
title(['Relativt fel: ' num2str(error)])
xlabel('x')
else
    error('Matrisen är positivt definit; en matris A sådan att x^(T) Ax > 0 för alla x');
end







