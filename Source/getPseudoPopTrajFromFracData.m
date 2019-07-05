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

function synth_ts = getPseudoPopTrajFromFracData(numEntities,bootstrapData,timepoints,bootstrapSize)

for e = 1:size(bootstrapData,3)
    tmp_bs = squeeze(bootstrapData(:,:,e));
    tmp_synth_ts = [];
    numBatches = size(tmp_bs,1);
    
    for i = 1:size(timepoints,2)
        tmpFrac = zeros(1,numEntities);
        start_index = (i-1)*numEntities + 1;
        end_index = i*numEntities;
        frac = tmp_bs(:,start_index:end_index);
        
        for fp = 1:bootstrapSize
            b = randi(numBatches,1);
            exp_frac = [0 frac(b,:)];
            frac_range = cumsum(exp_frac);
            p = rand;
            entity = find(histc(p,frac_range));
            tmpFrac(entity) = tmpFrac(entity) + 1;
        end
        
        tmpFrac = tmpFrac/bootstrapSize;
        tmp_synth_ts = [tmp_synth_ts tmpFrac];
        
    end  
    synth_ts(:,:,e) = tmp_synth_ts;
    
end
end