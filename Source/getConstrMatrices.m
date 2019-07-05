%{ 
This file is part of TRANSCOMPP
Copyright (C) 2019  N. Suhas Jagannathan

TRANSCOMPP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

TRANSCOMPP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TRANSCOMPP.  If not, see <https://www.gnu.org/licenses/>.
%}

function [Aeq, Beq, lb, ub] = getConstrMatrices(numEntities,bd_constr,trans_constr,solveForBD)

s = numEntities^2;
if solveForBD
    s = numEntities^2 + numEntities;
end

Aeq = zeros(numEntities,s);

for i = 1:numEntities
    start_index = (i-1)*numEntities+1;
    Aeq(i,start_index:start_index+numEntities-1)=1;
end
Beq = ones(numEntities,1);



lb = 0.5*eye(numEntities); lb = reshape(lb,1,numEntities^2); 
ub = ones(1,numEntities^2);

for i = 1:size(trans_constr,1)
    r = trans_constr(i,1); c = trans_constr(i,2);
    index = (r-1)*numEntities + c;
    lb(index) = trans_constr(i,3);
    ub(index) = trans_constr(i,4);
end

if solveForBD
   lb_bd = 0.1+zeros(1,numEntities); ub_bd = 10*ones(1,numEntities);
   lb_bd(1) = 1; ub_bd(1) = 1;
   for i = 1:size(bd_constr,1)
       index = bd_constr(i,1);
       lb_bd(index) = bd_constr(i,2);
       ub_bd(index) = bd_constr(i,3);
       
   end    
   lb = [lb lb_bd]; ub = [ub ub_bd];
end

end