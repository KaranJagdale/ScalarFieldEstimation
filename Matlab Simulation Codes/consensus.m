function consensus = consensus(a_cap,edgeWeight,i)
    consensus = zeros(length(a_cap(:,1)),1);
    for j=1:4
        if j~=i
            consensus = consensus + (a_cap(:,i)-a_cap(:,j))*edgeWeight(i,j);
        end    
    end       
end