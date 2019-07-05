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

function [transMat,bd_solved] = runSimulation(numEntities,TimeSeriesData,timepoints,numTriesOpt,trans_constr,solveForBD,bd_vector,bd_constr,errorTerm,lts_frac,solveType,numSim,bootstrapSize,scStateDef)
if isempty(bd_vector)
    bd_vector = ones(numEntities,1);
end
[Aeq, Beq, lb, ub] = getConstrMatrices(numEntities,bd_constr,trans_constr,solveForBD);

if solveType == 1
    [transMat, ~, ~,bd_solved] = solveForTransitionMatrix(numEntities,TimeSeriesData,timepoints,numTriesOpt,bd_vector,errorTerm,Aeq,Beq,lb,ub,solveForBD,lts_frac);
    return;
else
    transMat = zeros(numSim,numEntities^2);
    bd_solved = zeros(numSim,numEntities);
    wb = waitbar(0,'Processed 0 samples','Name','RESAMPLING PROGRESS','WindowStyle','modal');
    for sim = 1:numSim
        waitbar((sim/numSim),wb);
        if(mod(sim,10)==0)
            waitbar((sim/numSim),wb,['Processed ',num2str(sim),'/',num2str(numSim),' pseudo samples']);
        end
        if solveType == 2
            synth_ts = getPseudoPopTrajFromSCdata(numEntities,TimeSeriesData,timepoints,bootstrapSize,scStateDef);
        elseif solveType == 3
            synth_ts = getPseudoPopTrajFromFracData(numEntities,TimeSeriesData,timepoints,bootstrapSize);
        end
        
        [best_trans, ~, ~, best_bd] = solveForTransitionMatrix(numEntities,synth_ts,timepoints,numTriesOpt,bd_vector,errorTerm,Aeq,Beq,lb,ub,solveForBD,lts_frac);
        transMat(sim,:) = reshape(best_trans',1,numEntities^2);
        if solveForBD
            bd_solved(sim,:) = reshape(best_bd,1,numEntities);
        end
    end
    delete(wb);
end

end