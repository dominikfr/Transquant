function [imsegm,imnucleus,labelim,did,cm,nuclei,rim,controlim] = getimageobj(plate,prm)
% GETIMAGEOBJ Finding image objects of interest in the lowres super image 
%   for high resolution imaging.
%
%   [IMSEGM,IMNUCLEUS,LABELIM,DID,CMOUT,NUCLEI,RIM,CONTROLIM] = GETIMAGEOBJ(PLATE,PRM)
%   takes the arguments PLATE and PRM as in SCREEN, and finds objects of
%   interest for imaging. Here, we find DiD labelled objects that are sufficiently far
%   away from other DiD labelled objects, with a certain density of cells around, and
%   sufficiently far from the rim. Also, up to a specified number of
%   control regions are found. 
%
%   The algorithm returns the image used for segmentation IMSEGM, the nucleus
%   image IMNUCLEUS, the label image (to know where we are in the lowres 
%   scan field) LABELIM, the DiD image DID, the center of mass CM of 
%   the objects in pixels (used for highres imaging), the found nuclei
%   NUCLEI, the found rim RIM, and the control image CONTROLIM showing 
%   the fields that were used as control and non-control regions.
%
%   The only mandatory return argument for SCREEN is the CM variable
%   with coordinates for highres imaging. Must be given as nx2 array named 
%   CM.ALL.LABEL and CM.ALL.REL
%   
%   This file is project specific
%   and is therefore considered custom made. Thus, it should be considered
%   as an example of algorithm for finding objects in the lowres image for
%   imaging in the highres image.
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
msg = ['This is ' upper(mfilename) ' finding highres image objects'];
disp(msg);


% take out the image
imsegm = plate.lowresim{prm.segmch};
imnucleus = plate.lowresim{prm.nucleusch};
labelim = plate.labelim;
clear plate:
dim = size(imsegm);

% make local coordinate image
a = 1:prm.lowres.dim(1):dim(1);
b = 1:prm.lowres.dim(2):dim(2);
na = numel(a);
nb = numel(b);
d = prm.lowres.dim/2-1;

% make centered grid
[cx,cy] = ndgrid(-d(1)-0.5:d(1)+0.5,-d(2)-0.5:d(2)+0.5);
cx = imrotate(cx,-90);
cy = imrotate(cy,-90);
cx = single(cx);
cy = single(cy);
cimx = repmat(cx,na,nb);
cimy = repmat(cy,na,nb);

% find nuclei
tic
msg = ['Finding nuclei'];
disp(msg);
nuclei.bw = imnucleus > prm.detect.levelnuclei;
th = round(prm.detect.smallestnucleiarea/prm.lowres.voxelvol);
nuclei.bw = bwareaopen(nuclei.bw,th);
th = round(prm.detect.largestnucleiarea/prm.lowres.voxelvol);
nuclei.bw = logical(nuclei.bw - bwareaopen(nuclei.bw,th));
toc;


% find donor cells
tic
msg = ['Finding donor cells'];
disp(msg);
did.bw = imsegm > prm.detect.leveldid;
% did.bw = im2bw(imsegm,I/2);
se = strel('disk',2);
did.bw = imclose(did.bw,se);
% keep only the medium large objects
th = round(prm.detect.smallestcellarea/prm.lowres.voxelvol);
did.bw = bwareaopen(did.bw,th);
th = round(prm.detect.largestcellarea/prm.lowres.voxelvol);
did.bw = logical(did.bw - bwareaopen(did.bw,th));
toc

% find rim
tic
msg = ['Finding rim'];
disp(msg);
rim.bw = imnucleus > prm.detect.levelrim;
th = round(10000/prm.lowres.voxelvol);
rim.bw = bwareaopen(rim.bw,th);
[rim.c(:,1),rim.c(:,2)] = find(rim.bw);
rim.nc = size(rim.c,1);
% true voxel size
toc




% find the center of mass of these objects
tic
[faser,L] = bwlabeln(did.bw);
msg = ['Initially finding ' int2str(L) ' donorcells'];
disp(msg);
cm = NaN(L,2);
for i = 1 : L
    clear c;
    [c(:,1),c(:,2)] = find(faser == i);
    % center of mass in voxels
    cm(i,:) = mean(c,1);    
end;
toc

tic
% only keep objects that are at a certain distance to all other objects
remove = zeros(L,1);
remhere = zeros(L,1);
for i = 1 : L
    for j = i+1 : L
        d = cm(j,:) - cm(i,:);
        d = d.*prm.lowres.h;
        absd = norm(d);
        if absd < prm.detect.smallestdistcells
            remove(i) = 1;            
            remove(j) = 1;
            remhere(i) = 1;            
            remhere(j) = 1;            
        end;
    end;
end;
msg = ['Cell-to-cell distance removed ' int2str(sum(remhere)) ' objects'];
disp(msg);
toc


tic
% only keep cells that are in a dense collection of cells
d = prm.detect.averagecelldiameter;
remhere = zeros(L,1);
for i = 1 : L
    cmhere = cm(i,:);    
    cmhere = round(cmhere);
    m1 = max(1,cmhere(1)-d);m2 = min(dim(1),cmhere(1)+d);
    n1 = max(1,cmhere(2)-d);n2 = min(dim(2),cmhere(2)+d);
    imhere = nuclei.bw(m1:m2,n1:n2);
    [a,n] = bwlabeln(imhere);
    if n < prm.detect.minnumneigh
        remove(i) = 1;
        remhere(i) = 1;
    end;
end;
msg = ['Density requirement removed ' int2str(sum(remhere)) ' objects'];
disp(msg);
toc

tic
% only keep the objects that have a certain distance to the rim
remhere = zeros(L,1);
for i = 1 : L
    cmhere = cm(i,:);
    d = rim.c - repmat(cmhere,rim.nc,1);  
    d = d.*repmat(prm.lowres.h,rim.nc,1);
    d = d.^2;
    d = sum(d,2);
    d = sqrt(d);   
    d = min(d);
    if d < prm.detect.smallestdistrim
        remove(i) = 1;            
        remhere(i) = 1;
    end;
end;
msg = ['Cell-to-rim distance removed ' int2str(sum(remhere)) ' objects'];
disp(msg);
toc
        
% Remove objects
ind = find(remove);
if ~isempty(ind)
    msg = ['Removing objects '];
    disp(msg);
    
    % remove from CM
    cm(ind,:) = [];    
    for i = 1 : numel(ind)            
        % remove from BW image
        did.bw(faser == ind(i)) = 0;
    end;
end;

msg = ['Removed in total ' int2str(sum(remove)) ' objects'];
disp(msg);

% H = figure('Position',[100 100 1000 800]);
% subplot(2,3,1);imagesc(imsegm);colormap(gray);axis image;axis off;
% subplot(2,3,2);imagesc(imnucleus);colormap(gray);axis image;axis off;
% subplot(2,3,3);imagesc(nuclei.bw);colormap(gray);axis image;axis off;
% subplot(2,3,4);imagesc(rim.bw);colormap(gray);axis image;axis off;
% subplot(2,3,5);imagesc(bw);colormap(gray);axis image;axis off;

% Find the centered coordinates of the hits
cm = round(cm);
n = size(cm,1);
cmout.did.label = cm;
cm2 = zeros(n,2);
for i = 1 : size(cm,1)
    cm2(i,1) = cimx(cm(i,1),cm(i,2));
    cm2(i,2) = cimy(cm(i,1),cm(i,2));
end;
cmout.did.rel = cm2;

%
% Get control images
%
[cmout.control.label,controlim] = matrixcompute.controlroi(imsegm,labelim,nuclei,rim,prm);
cm = cmout.control.label;
n = size(cmout.control.label,1);
cm2 = zeros(n,2);
for i = 1 : n
    cm2(i,1) = cimx(cm(i,1),cm(i,2));
    cm2(i,2) = cimy(cm(i,1),cm(i,2));
end;
cmout.control.rel = cm2;

% all coordinates and labels   
cmout.all.rel = [cmout.did.rel;cmout.control.rel];
cmout.all.label = [cmout.did.label;cmout.control.label];
cm = cmout;

