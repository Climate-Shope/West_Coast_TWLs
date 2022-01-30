function c= negloglikgev(theta,data)
	
		y = 1 + (theta(1) * (data - theta(3)))./theta(2);
		if((theta(2) < 0) | (min(y) < 0))
			c= inf;
		else 
			term1 = length(data) * log(theta(2));
			term2 = sum((1 + 1/theta(1)) * log(y));
			term3 = sum(y.^(-1/theta(1)));
			c = term1 + term2 + term3;
        end    
		   
