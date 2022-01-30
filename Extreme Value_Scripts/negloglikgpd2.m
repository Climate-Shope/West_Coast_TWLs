function f=negloglikgpd2(theta,excess,T);
    
    nu = theta(1);
    xi = theta(2);
    sigma = theta(3);

    
    cond1 = sigma <= 0;
    cond2 = (xi <= 0) & (max(excess) > ( - sigma/xi));
    if(cond1 | cond2)
		f=inf;
    else

        y = (1+1/xi)*sum(log(1+xi*excess/sigma));
        
        z = -(length(excess)*log(nu)) + nu*T;
        
        f = length(excess)*log(sigma) + y + z;
    end
    
   


