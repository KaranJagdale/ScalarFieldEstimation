function edge = adjacentvertex(bot_locx,bot_locy,i,j)
% finds if i and j are adjacent vertices
    [xbord1, ybord1] = compute_voronoi(i, [-5 5 5 -5], [-5 -5 5 5], bot_locx, bot_locy);
    [xbord2, ybord2] = compute_voronoi(j, [-5 5 5 -5], [-5 -5 5 5], bot_locx, bot_locy);
    for i=1:length(xbord1)
        for j = 1:length(xbord2)
            xbord1(i) = round(10.^4.*xbord1(i))./(10.^4);
            xbord2(j) = round(10.^4.*xbord2(j))./(10.^4);
        end    
    end    
    for i=1:length(ybord1)
        for j = 1:length(ybord2)
            ybord1(i) = round(10.^4.*ybord1(i))./(10.^4);
            ybord2(j) = round(10.^4.*ybord2(j))./(10.^4);
        end    
    end    
    xcommon = intersect(xbord1,xbord2);
    %disp('xcommon')
    %disp(xcommon);
    if(i == 1)
        %disp(xbord1);
        %disp(xbord2);
        %disp(xcommon);
    end  
    if (j ==3)
        if(xbord1(3)==xbord2(2))
            %disp('ok');
        end    
    end    
    %disp('xcommon');
    %disp(xcommon);
    n = length(xcommon);
    i = 1;
    count = 0;
    matcher = 0;
    
    while(i<=n)
        index1 = find(xbord1 == xcommon(i));
        index2 = find(xbord2 == xcommon(i));
        %disp(index1);
        %disp(index2);
        for l = 1:length(index1)
            for m = 1:length(index2)
                if (ybord1(index1(l))==ybord2(index2(m)))
                    matcher =1;
                end
            end
        end
        if (matcher == 1)
            count = count +1;
        end    
            
        i = i+1;
    end    
            
    if (count >= 2)
        edge = 1;
    else 
        edge = 0;
    end
end
