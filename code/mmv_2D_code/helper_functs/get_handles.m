function [A,AH] = get_handles(F,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function turns a matrix into a function handle to be used in the
% optimization code 
% input: 
%         F         matrix
% outputs: 
%         A         function handle associated with the forward operator
%         AH        function handle associated with the adjoint operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = @(u) A_handle(F,N,u); 
AH = @(u) AH_handle(F,N,u); 

end

function y = A_handle(F,N,u)

u = reshape(u,N,N);
y = F*u;
y = y(:); 

end

function y = AH_handle(F,N,u) 

u = reshape(u,N,N);
% y = (u'.*F)';
y = F'*u;
% y = F\u;
y = y(:); 

end