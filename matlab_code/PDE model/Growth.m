function g = Growth(x,r_max,K) %withinhost pathogen growth
         g = r_max .* (1-exp(x)./K);   
end
