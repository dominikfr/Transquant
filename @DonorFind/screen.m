% SCREEN Script for automatic screening on a Leica SP5 microscope in 
%   combination with MatrixScreener and LAS AF. The
%   communication goes via the CAM in the MatrixScreener software. 
%   This function has on purpose been made as a script in order to be able to
%   execute separate parts independently. The connection to the microscope
%   computer must be made as a TCP/IP connection via port 8895. This software
%   has been made for two computers: (1) The miscroscope computer (here 
%   called the server) and the (2) analysis computer (here called the
%   client).
%
%   There are three jobs for lowres imaging that need to match between this 
%   file and MatrixScreener: 
%
%   1. Low resolution job defined in PRM.LOWRES.JOBNAME,
%
%   2. Waiting job defined in PRM.LOWRES.JOBNAMEWAIT, 
%
%   3. High resolution job defined in PRM.HIGHRES.JOBNAME. 
%
%   The lowres job is the low resolution imaging for screening the wells for 
%   interesting imaging regions. The waiting job is a job essentially doing 
%   nothing, just letting the microscope enter a waiting status, waiting 
%   for the client to analyse the images. There must be one waiting after 
%   each line of wells. The highres job is the high resolution job, 
%   typically in 3D, that you are finally aiming at.
%
%   The scanning is done one per line in order to restrict the number of 
%   wells acquired simultaenously. After one line the analysis for 
%   interesting regions is performed. Too many wells results often in drift 
%   of microscope and autofocus map. After each line there must be a 
%   waiting job.
%
%   The objects for imaging are defined in GETIMAGEOBJ. This file has been
%   tuned to a specific project and must, depending on the project, 
%   be replaced by a custom made software for detection of ROIs for imaging.
%   This is typically very depending on the biological question.
%   
%
%   Check list before imaging:
%   -------------------------
%   * Map a network drive from the client to the server (physically this is 
%   done on the server) such that acquired images are exported directly to a 
%   shared drive on the client.
% 
%   * Open a connection from the client to the server using the function
%   OBJ = OPENCONNECTION.
%
%   * Set the line numbers for scanning, e.g. PRM.LOWRES.LINE = [2,3].
%
%   * The software is using global thresholds in order to be fast. 
%   Remember to set a suitable threshold.
%
%   * Are the correct channels set that MATLAB should read? Set for instance 
%   by PRM.LOWRES.CH = {'0','1'} for reading channel zero (starting counting
%   at 0) and one.
% 
%   * Is the correct export path set in MatrixScreener. Also set the
%   alternative export path to the same.
%
%   * Counting in matlab starts from one and the channels for analysis are
%   named accordingly. For instance, set a segmentation channel in the 
%   second channel (counted as channel one in PRM.LOWRES.CH) as 
%   PRM.SEGMCH = 2.
% 
%   * Remember to set the correct IP address of the server in the file 
%   OPENCONNECTION
%
%   * Remember to set the correct path on the client where the images can
%   be found. Set in PRM.LOWRES.EXPORTPATH as a full path.
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


% number of lines
prm.lowres.line = [1,2];
prm.lowres.nline = numel(prm.lowres.line);
% the channels to read (NB: one number down from the segmentation and 
% nucleus channel, counting starts at 0)
prm.lowres.ch = {'0','1'};

% number of channels
prm.lowres.nch = numel(prm.lowres.ch);

% voxel size in the lowres job
prm.lowres.h = [1.52,1.52];
% voxel volume in the lowres job
prm.lowres.voxelvol = prod(prm.lowres.h);

% Dimension of lowres images in pixels
prm.lowres.dim = [512,512];

% the lowres job name, job number one above
prm.lowres.jobname = 'lowres';

% the job number indicating that we have reached the end of well, job
% number two above
prm.lowres.jobnamewait = 'wait';

% highres job name, job number three above
prm.highres.jobname = 'jobhighres';

% path where the images are exported to (mapped drive on the client)
prm.lowres.exportpath = 'C:\Users\SP5\Documents\MatrixScreenerImages\';


% Folder to save information from the screen
% Modify this according to your wish
prm.foldersave = fullfile('..',['evaluation-' date]);
[a,b,c] = mkdir(prm.foldersave);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BEGIN Custom made settings for a specific type of recognition of ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Settings for finding the candidates of interest in the lowres job
% 

% define channels
prm.nucleusch = 1;
prm.segmch = 2;

% level for thresholding to find DiD cells in lowres job
prm.detect.leveldid = 250;

% level for thresholding to find nuclei in lowres job. 
prm.detect.levelnuclei = 100;

% level for thresholding to find rim in lowres job (quite stable with this threshold)
prm.detect.levelrim = 130;

% All distance defined in microns

% smallest distance between did cells
prm.detect.smallestdistcells = 200;

% smallest distance to the rim
prm.detect.smallestdistrim = 300;

% smallest and largest possible nuclei area
prm.detect.smallestnucleiarea = 80;
prm.detect.largestnucleiarea = 700;

% smallest and largest cell area (used for DiD cell definition)
prm.detect.smallestcellarea = 300;
prm.detect.largestcellarea = 8000;

% average cell diameter (used for density measurements)
prm.detect.averagecelldiameter = 110;

% Minimum number of neighbors
prm.detect.minnumneigh = 15;

% number of control images
prm.detect.ncontrol = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   END Custom made settings for a specific type of recognition of ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% start lowres scan. This starts the low resolution scan if everything is
% set up correctly
send = startscan;
matrixcompute.sendcmd(obj,send);

% clear old data from previous scans
clear labelim im bw plate lowresim

imagelabel = cell(prm.lowres.nline,1);
cmdsave = cell(prm.lowres.nline,1);
infosave = cell(prm.lowres.nline,1);
for i = 1 : prm.lowres.nline
    
    msg = ['Well number ' int2str(i)];
    disp(msg);
    
    % get lowres paths. NB: reading images in GETIMAGEOBJ
    plate = matrixcompute.getlowres(prm,obj);
    
    % find objects to image
    % NB this function is custom made and must be replaced by any finding
    % objects of interest
    clear cm;
    [imsegm,imnucleus,labelim,did,cm,nuclei,rim,controlim] = matrixcompute.getimageobj(plate,prm);
    
    % reorder the coordinates in random order to avoid any bias
    n = size(cm.all.rel,1);
    v = randperm(n);
    cm.all.rel = cm.all.rel(v,:);
    cm.all.label = cm.all.label(v,:);    
    
    msg = ['Printing images'];
    disp(msg);    
    matrixcompute.printfile(imsegm,'DiDimage',prm,i);
    matrixcompute.printfile(did.bw*256,'donorcells',prm,i);
    matrixcompute.printfile(imnucleus,'nucleiimage',prm,i);
    matrixcompute.printfile(nuclei.bw*256,'nuclei',prm,i);
    matrixcompute.printfile(controlim*128,'controlim',prm,i);
    matrixcompute.printfile(rim.bw*256,'rim',prm,i);
        
    % start highres
    cmdsave{i,1} = matrixcompute.starthighres(obj,cm.all.rel,cm.all.label,...
        plate.info,plate.labelim,prm.highres.jobname);
    infosave{i,1} = plate.info;
end

% get the experimental data and save them
pathsave = fullfile(prm.foldersave,'experiment.mat');
save(pathsave,'prm','cmdsave','infosave');

  