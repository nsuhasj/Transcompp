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

function [trans, trans_batch, fval_batch, bd_solved] = solveForTransitionMatrix(numEntities, TimeSeriesData,timepoints,numTries,bd_vector,errorTerm,Aeq,Beq,lb,ub,solveForBD,lts_frac)

if(isempty(find(timepoints==0,1)))
    disp ('Please ensure your time-series data includes FACS efficiency as t0 data');
    return;
end

min_f = 100;
trans_batch = zeros(numTries,numEntities*numEntities);
fval_batch = zeros(numTries,1);
f_tol = 0.001;
bd_solved = [];


for index = 1:numTries
    temp_trans = zeros(numEntities);
    for j = 1:numEntities-1
        sum_trans = sum(temp_trans,2);
        temp_trans(:,j) = (1-sum_trans).*rand(numEntities,1);
    end
    temp_trans(:,numEntities) = 1-sum(temp_trans,2);
    
    x0 = reshape(temp_trans',1,numEntities*numEntities);
    if solveForBD
        bd_min = lb(end-numEntities+1:end);
        bd_max = ub(end-numEntities+1:end);
        r = rand(size(bd_min));
        bd_toAdd = bd_min + r.*(bd_max-bd_min);
        bd_toAdd = reshape(bd_toAdd,1,numEntities);
        x0 = [x0 bd_toAdd];
    end
    options = optimset('Algorithm','interior-point', 'MaxFunEvals',50000,'TolFun',0.00000000005,'TolCon',0.00000000001,'display','off');
    [x,fval] = fmincon(@(x)objSolveforTrans(numEntities,x,TimeSeriesData,timepoints,bd_vector,errorTerm,solveForBD,lts_frac),x0,[],[],Aeq,Beq,lb,ub,[],options);
    fval_batch(index) = fval;
    trans_batch(index,:) = x(1:numEntities^2);
    
    if fval< min_f
        min_f = fval;
        trans = reshape(trans_batch(index,:),numEntities,numEntities); trans = trans';
        if solveForBD
            bd_solved = x(numEntities^2+1:end);
        end
    end
    
    if min_f < f_tol
        break;
    end
end

end
