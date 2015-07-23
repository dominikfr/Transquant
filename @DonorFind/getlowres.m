function [plate] = getlowres(prm,obj)
% GETLOWRES Reading the relative paths of lowres images sent from
%   MatrixScreener, then loading the images and stiching them together.
% 
%   PLATE = GETLOWRES(PRM,OBJ) Takes in the parameter argument PRM from 
%   SCREEN as well as the TCP/IP object OBJ.
%
%   Returning the variable struct array PLATE with info about each of the
%   low res images as well as the super images, the stiched together low
%   res images.
%
%   NB: This software is not made for overlap between the lowres images, this
%   would need an image registraiton algorith to optimize the stiching.
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
% Listen for images until waiting job
%
imcount = 0;
mes = '';
while 1
    
    % listen for image arriving
    mesold = mes;
    while isequal(mesold,mes)    
        [mes,count] = fscanf(obj);    
    end;
    msg = ['Receiving ' mes];
    disp(msg);

    % get image info from the specific code used by MatrixScreener
    lowres = matrixcompute.getimageinfo(mes);    
            
    % waiting job, lets go on!
    if isequal(prm.lowres.jobnamewait,lowres.jobnameread)
        msg = ['Finding terminating job number, returning'];
        disp(msg);
        break;
    end;

    % did we get a highres image? Do nothing until we get a valid lowres
    % image with the correct job name
    if ~isequal(prm.lowres.jobname,lowres.jobnameread)
        disp('Received image not a lowres image, waiting');
        continue;
    end;

    % we received a lowres job
    msg = ['Receiving lowres job'];
    disp(msg);
    
    % count one more
    imcount = imcount + 1;    
    
    % store the read information
    plate.info{imcount}.lowres = lowres;
    
    
end;

msg = ['Received ' int2str(imcount) ' images'];
disp(msg);

save test

%
% Loading the images and placing the together from here and to the end
% 


% Find max and min coordinate to know the extension of the super image
n = numel(plate.info);
X = NaN(n,1);Y = NaN(n,1);
for i = 1 : n
    X(i) = plate.info{i}.lowres.attributen.Y;
    Y(i) = plate.info{i}.lowres.attributen.X;
end;
minx = nanmin(X)+1;
maxx = nanmax(X)+1;
miny = nanmin(Y)+1;
maxy = nanmax(Y)+1;
dim = [maxx-minx+1,maxy-miny+1];
% the offset, where to start, relatively
offsetx = minx-1;
offsety = miny-1;


% Collect images and stich them together
plate.lowresim = cell(prm.lowres.nch,1);
for i = 1 : prm.lowres.nch
    % if the images are not uint8, this is not good
    plate.lowresim{i} = zeros(dim(1)*prm.lowres.dim(1),dim(2)*prm.lowres.dim(2),'uint8');
end;
% single format to save memory
plate.labelim = zeros(dim(1)*prm.lowres.dim(1),dim(2)*prm.lowres.dim(2),'single');

% collect images
msg = ['Collecting images and stiching them together'];
disp(msg);
% Read the images
nimages = numel(plate.info);
for i = 1 : nimages
    
    % this field information
    lowres = plate.info{i}.lowres;
    
    % Find where the channel information is to replace by the right channel
    ind = strfind(lowres.relpath,'--C');
    
    % channels to use        
    im = cell(prm.lowres.nch,1);
    for j = 1 : prm.lowres.nch
        % the channel to use
        % not gonna work for more than 9 channels...who cares?
        lowres.relpath(ind+4) = prm.lowres.ch{j};    
    
        % read the image
        pathload = [prm.lowres.exportpath  lowres.relpath];  
        msg = ['Reading ' pathload];
        disp(msg);
        im{j} = imread(pathload);
        % NB: LAS Af rotates the image!
        im{j} = imrotate(im{j},-90);
    
    end;    
    lowres.dim = size(im{1});

    % put the image into the mesh of the plate    
    x = lowres.attributen.Y-offsetx;
    y = lowres.attributen.X-offsety;

    % put the image into the lowres image
    x = x*lowres.dim(1) + 1;
    y = y*lowres.dim(2) + 1;

    for j = 1 : prm.lowres.nch
        plate.lowresim{j}(x:x+lowres.dim(1)-1,y:y+lowres.dim(2)-1) = uint8(im{j});
    end;
    % the label image to know the field for later use in the highres. This
    % format can not be uint8, too few values!!
    plate.labelim(x:x+lowres.dim(1)-1,y:y+lowres.dim(2)-1) = i*ones(lowres.dim,'single');
        
end;




