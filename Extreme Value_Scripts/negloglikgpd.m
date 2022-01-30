function f=negloglikgpd(theta,excess),
    xi = theta(1);
    
    beta = theta(2);
	
    cond1 = beta <= 0;
    cond2 = (xi <= 0) & (max(excess) > ( - beta/xi));
    if(cond1 | cond2)
		f=inf;
    else

        y = (1+1/xi)*sum(log(1+xi*excess/beta));
        
        f = length(excess) * log(beta) + y;
    end

