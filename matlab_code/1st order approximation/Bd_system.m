function [S,I,P,Z,load,R0,stability]=Bd_system(S_hat,parms)
         equ=arrayfun(@(n) Bd_equlibrium(parms(n, 1), parms(n, 2), parms(n, 3), parms(n, 4), ...
                                      parms(n, 5), parms(n, 6), parms(n, 7), parms(n, 8), ...
                                      parms(n, 9), parms(n, 10), parms(n, 11), parms(n, 12), parms(n, 13)), ...
                                      1:size(parms, 1), 'UniformOutput', false);
         equ=cell2mat(equ)';
         S=equ(:,1);I=equ(:,2);P=equ(:,3);Z=equ(:,4);load=P./I;
         R0=arrayfun(@(n) Bd_R0(S_hat,I(n),P(n),parms(n, 1), parms(n, 2), parms(n, 3), parms(n, 4), ...
                                      parms(n, 5), parms(n, 6), parms(n, 7), parms(n, 8), ...
                                      parms(n, 9), parms(n, 10), parms(n, 11)), ...
                                      1:size(parms, 1), 'UniformOutput', false);
         R0=cell2mat(R0)';

         stability=arrayfun(@(n) Bd_stability(S(n),I(n),P(n),Z(n),parms(n, 1), parms(n, 2), parms(n, 3), parms(n, 4), ...
                                      parms(n, 5), parms(n, 6), parms(n, 7), parms(n, 8), ...
                                      parms(n, 9), parms(n, 10), parms(n, 11), parms(n, 12), parms(n, 13)), 1:size(parms, 1), 'UniformOutput', false);
         stability=cell2mat(stability)';
end