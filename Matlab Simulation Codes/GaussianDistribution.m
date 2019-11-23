function localgaussian = GaussianDistribution(c,bot_locx,bot_locy,botNo)
%This function takes value of location of bots gaussian centres and outputs
%the gaussian centres lying in voronoi region of bot
    n1 = length(bot_locx);
    n = length(c);
    product = true;
    localgaussian = [];
    condition = [];
    for i = 1:n
        for j = 1:n1
            if (j ~=botNo)
                %disp(j);
                a = ((c(i,1)-(bot_locx(j)+bot_locx(botNo))./2).*(bot_locx(j)-bot_locx(botNo))./(bot_locy(botNo)-bot_locy(j)))-(c(i,2)-(bot_locy(j)+bot_locy(botNo))./2);
                b = ((bot_locx(botNo)-(bot_locx(j)+bot_locx(botNo))./2).*(bot_locx(j)-bot_locx(botNo))./(bot_locy(botNo)-bot_locy(j)))-(bot_locy(botNo)-(bot_locy(j)+bot_locy(botNo))./2);
                condition(j) = a*b;
            else
                condition(j) = 0;
                
            end
        end
        %disp(condition);
        if (condition(1)>=0 && condition(2)>=0 && condition(3)>=0 && condition(4)>=0)
             localgaussian = [localgaussian; c(i,1) c(i,2)];
        
        end
    end 
     %disp('botNo');
     %disp(botNo);
end