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

function synth_ts = getPseudoPopTrajFromSCdata(numEntities,bootstrapData,timepoints,bootstrapSize,scStateDef)

for e = 1:size(bootstrapData,3)
    tmp_bs = squeeze(bootstrapData(:,:,e));

    numBatches = size(tmp_bs,1);
    tmp_synth_ts = [];
    for i = 1:size(timepoints,2)

        tmpFacsPts = zeros(bootstrapSize,size(tmp_bs{1,1},2));
        
        for fp = 1:bootstrapSize
            b = randi(numBatches,1);
            exp_pts = tmp_bs{b,i};
            rand_index = randi(size(exp_pts,1));
            tmpFacsPts(fp,:) = exp_pts(rand_index,:);
        end
        
       % tmpFrac = getFracFromPts(numEntities,tmpFacsPts);
       tmpFrac = getPopFrac(scStateDef,tmpFacsPts');
       tmp_synth_ts = [tmp_synth_ts tmpFrac];
    end
    synth_ts(:,:,e) = tmp_synth_ts;
end
end