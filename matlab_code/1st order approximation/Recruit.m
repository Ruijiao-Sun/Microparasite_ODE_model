function dndt = Recruit(N,r,gamma)
         dndt = r.*N.*exp(-gamma.*N);
end
