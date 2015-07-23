function [cmdsave] = starthighres(obj,cm,cmlabelim,info,labelim,jobname)
% STARTHIGHRES Starting high resolution imaging based on image coordinates
% collected previously.
%   
%   CMDSAVE = STARTHIGHRES(OBJ,CM,CMLABELIM,PLATE,LABELIM,JOBNAME) makes
%   the list of high resolution jobs and sends the starting command. OBJ is
%   the TPI/IP object, CM is the local coordinates within the filed,
%   CMLABELIM is the coordinates within the super image, showing us which
%   field we are in, PLATE is from SCREEN and contains information about
%   the lowres fields, LABELIM is the superimage labelled with the lowres
%   fields, and JOBNAME is the highres jobname.
%   
%   NB: We need both coordinates within the field (CM) and
%   the coordinates of the label image (CMLABELIM) so we know which lowres
%   field we are within. The coordinates within the field are only local
%   coordinates.
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
%

% Delete list 
send = deletelist;
matrixcompute.sendcmd(obj,send);

% make sure all coordinates are rounded numbers
cm = round(cm);

cmdsave = cell(size(cm,1),1);
n = size(cm,1);
msg = ['Number of highres fields: ' int2str(n)];
disp(msg);
for i = 1 : n

    % which lowres image are we in? This is a lowres scan field numbering
    x = round(cmlabelim(i,1)); 
    y = round(cmlabelim(i,2));
    label = labelim(x,y);

    % these are the local image coordiantes within the image
    cmhere = cm(i,:);

    msg = ['Adding to list coordinates ' num2str(cmhere)];
    disp(msg);

    % add the job to cam list only if we are within a lowres field, just an
    % extra check
    if label == 0
        continue;
    end;
    send = matrixcompute.addtocamlist(info{label}.lowres,jobname,cmhere);
    matrixcompute.sendcmd(obj,send);      
    
    % keep the information to restore which well the images belong to
    cmdsave{i,1} = send;
    
end;

% "prepare" the scan
send = matrixcompute.startcamscan;
sendcmd(obj,send);        

% stop wait for camscan, starts the scan
send = matrixcompute.stopwaitcamscan();
sendcmd(obj,send);


