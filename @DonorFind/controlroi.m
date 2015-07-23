function [cmcontrol,controlim] = controlroi(imsegm,labelim,nuclei,rim,prm)
%
% CONTROLROI finds control regions for imaging based on special rules.
%
%   [CMCONTROL,CONTROLIM] = CONTROLROI(IMSEGM,LABELIM,NUCLEI,RIM,PRM) finds
%   control regions for imaging. The rules are custom made. Returning the
%   center of mass relative to the LABELIM for highres imaging. CONTROLIM
%   is an image showing where the control, possible non-chosen DiD targest,
%   and non control regions are positioned.
%   Possible DiD regions are labelled with "1" and taken control regions
%   are labelled as "2"
%
%   NB: This a file for a special project and should be replaced by a
%   suitable one!
% 
%
%   =======================================================================================
%   Copyright (C) 2014  Erlend Hodneland
%   Email: erlend.hodneland@gmail.com 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================
%

msg = ['This is ' upper(mfilename) ' finding control regions'];
disp(msg);

tic
% find possible did cells 
bw = imsegm > prm.detect.leveldid;

% remove noise objects
bw = bwareaopen(bw,10);

% label and find CM
[faser,L] = bwlabeln(bw);

msg = ['Finding ' int2str(L) ' potential DiD objects'];
disp(msg);

cm = zeros(L,2);
for i = 1 : L
    % this region
    reg = faser == i;
        
    clear c;
    [c(:,1),c(:,2)] = find(reg);

    % center of mass in voxels
    cm(i,:) = mean(c,1);
end;
clear faser;
cm = round(cm);
toc;

% distance around cells
dist = prm.detect.smallestdistcells;


% label image around the DiD cells not to take these regions
dim = size(labelim);
controlim = zeros(dim(1),dim(2),'uint8');
for i = 1 : L
    x = max(1,cm(i,1) - dist):min(dim(1),cm(i,1)+dist);
    y = max(1,cm(i,2) - dist):min(dim(2),cm(i,2)+dist);
    controlim(x,y) = 1;    
end;

msg = ['Stopping after finding ' int2str(prm.detect.ncontrol) ' control ROIs'];
disp(msg);

% count number of wells
[faser,nwells] = bwlabeln(labelim > 0);

% save test2
for j = 1 : nwells

    
    msg = ['Well number ' int2str(j)];
    disp(msg);
    
    % stay within this well
    reghere = faser == j;
    % and stay within the possible control pixels
    reghere = reghere .* (controlim == 0);

    clear cwell;
    [cwell(:,1),cwell(:,2)] = find(reghere);    
    ntot = size(cwell,1);
    % randomize order
    v = randperm(ntot);
    cwell = cwell(v,:);
    clear reghere;
    
    ncontrolwell = 0;
    ncontroltotal = 0;
    cmcontrol = zeros(prm.detect.ncontrol,2);    
    for i = 1 : ntot
        cmhere = cwell(i,:);
        x = max(1,cmhere(1)-dist):min(dim(1),cmhere(1)+dist);
        y = max(1,cmhere(2)-dist):min(dim(2),cmhere(2)+dist);
        reg = controlim(x,y);

        % is there an overlap to a DiD cell image?
        if ~isempty(find(reg,1))
            msg = ['Overlap to DID or earlier object'];
            disp(msg);        
            continue;
        end;

        % is it dense enough of cells?     
        imhere = nuclei.bw(x,y);
        [a,n] = bwlabeln(imhere);
        if n < prm.detect.minnumneigh
            msg = ['Too sparse cells'];
            disp(msg);
            continue;
        end;

        % is it far enough from rim?
        d = rim.c - repmat(cmhere,rim.nc,1);  
        d = d.*repmat(prm.lowres.h,rim.nc,1);
        d = d.^2;
        d = sum(d,2);
        d = sqrt(d);
        d = min(d);
        if d < prm.detect.smallestdistrim
            msg = ['Too close too rim'];
            disp(msg);
            continue;
        end;

        % label this region as "taken"
        controlim(x,y) = 2;

        % and store the coordinates
        ncontrolwell = ncontrolwell + 1;
        ncontroltotal = ncontroltotal + 1;
        cmcontrol(ncontroltotal,:) = cmhere;
        msg = ['Finding control ROI ' int2str(ncontrolwell)];
        disp(msg);

        if i == ntot
            msg = ['Could not find enough control regions'];
            disp(msg);
            break;
        end;
        
        if ncontrolwell == prm.detect.ncontrol
            msg = ['Reaching maximum number of control ROIS'];
            disp(msg);
            break;
        end;
        
    end;
end;




