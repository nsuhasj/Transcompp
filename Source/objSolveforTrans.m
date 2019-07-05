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

function f = objSolveforTrans(numEntities, x,TimeSeriesData,timepoints,bd_vector,errorTerm,solveForBD,lts_frac)

expected_bd_vector_size = [numEntities,1];

if solveForBD == 0
    if ~isequal(length(x),numEntities^2)
        disp ('Error in size of transition matrix seed');
        return;
    end
    
    if ~isequal(size(bd_vector),expected_bd_vector_size) && bd_vector ~= 1
        disp ('Error in Supplied BD vector');
        return;
    end
    trans_mat = reshape(x,numEntities,numEntities); trans_mat = trans_mat';
else
   bd_vector = x((numEntities^2+1):end);  
   trans_mat = reshape(x(1:(numEntities^2)),numEntities,numEntities); trans_mat = trans_mat';
end

bd_mat = diag(bd_vector);
index_t0 = find(timepoints==0,1);
numReplicates = size(TimeSeriesData,1);
f = 0;
f_array = [];

if errorTerm == 3 && lts_frac == 0
    lts_frac = 0.8;
end

for facs_sort = 1:size(TimeSeriesData,3)
    
    for sample = 1:numReplicates
        
        index_t0_ts = (index_t0-1)*numEntities+1;
        pop_t0 = squeeze(TimeSeriesData(sample,index_t0_ts:index_t0_ts+(numEntities-1),facs_sort));
        
        for time = 1: length(timepoints)
            if time == index_t0
                continue;
            end
            
            index_t_ts = (time-1)*numEntities+1;
            frac_pop_t = squeeze(TimeSeriesData(sample,index_t_ts:index_t_ts+(numEntities-1),facs_sort));
            
            if sum(frac_pop_t) == -1*numEntities
                continue;
            end

            scaleMat = bd_mat*trans_mat;
            pred_pop_t = pop_t0 * (scaleMat^timepoints(time));
            pred_frac_pop_t = pred_pop_t/sum(pred_pop_t);
            
            if errorTerm == 1
                f = f + sum((pred_frac_pop_t-frac_pop_t).^2);  %SSE
            elseif errorTerm == 2
                f = f + sum(abs(pred_frac_pop_t-frac_pop_t));  %L1
            elseif errorTerm == 3
                f_array = [f_array (pred_frac_pop_t - frac_pop_t).^2];
            elseif errorTerm == 4
                f = f + customError(pred_frac_pop_t,frac_pop_t);
            else
            end           
        end
    end
end

f_array = sort(f_array);
lts_max_ind = ceil(lts_frac*length(f_array));
if errorTerm == 3
    f = sum(f_array(1:lts_max_ind));   % LTS
end

end