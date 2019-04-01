function [ P ] = PA_Operator_1D( N,m_spa )

if m_spa == 0
    P = speye(N);
else
    
    C_spa = ones(m_spa+1,1);
    for j = 1:m_spa+1
        for js =1:m_spa+1
            if js ~= j
                C_spa(j) = C_spa(j)/(j-js);
            end
        end
    end
    L_SPA = zeros(N);
    m2 = floor((m_spa+1)/2);
    q_norm = sum(C_spa(1:m2,1));
    C_spa = C_spa./q_norm ;
    for j=1:N
        inds = j-m2:j+m_spa-m2;
        inds(find(inds<1))=inds(find(inds<1))+N;
        inds(find(inds>N))=inds(find(inds>N))-N;
        L_SPA(j,inds) = C_spa;
    end
    
    P = sparse(L_SPA)*2^(1-m_spa);
    
end
