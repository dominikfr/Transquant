function [] = TransQuant(folder,trmode,prm)


trth = prm.(trmode).th;
trch = prm.(trmode).ch;
msg = ['Running transfer with threshold ' num2str(trth) ' with mode ' trmode];
disp(msg);

vis = 0;

% for saving with name
[foldersavefinal,b,c] = fileparts(folder{1});
[filesavefinal,b,c] = fileparts(foldersavefinal);
if ~isempty(b)
    filesavefinal = b;
end;


prm.planetrack = 4;

nfolder = length(folder);
% nfolder = 2;
start = 1;
for i = start : nfolder
    
    % for evaluation
    foldereval1 = [folder{i} '/' 'eval' '/' 'single'];
    [a b c] = mkdir(foldereval1);
    foldereval2 = [folder{i} '/' 'eval' '/' 'multiple'];
    [a b c] = mkdir(foldereval2);

    % for storing
    stats.folder(i).name = folder{i};
    stats.folder(i).im = [];
    
    nim = 1000;
%     nim = 2;
    startim = 1;
    nim = 100;
    for j = startim : nim
        

        % load the data
        pathload = [folder{i} '/' 'stack' int2str(j) '-segm.mat'];                    
        try            
            D = load(pathload);
            msg = ['Loaded ' pathload];
            disp(msg);
        catch         
            msg = ['Could not load ' pathload ', continuing'];
            disp(msg);
            continue;
        end;
        
        
        im = D.im;
        prm.h = D.info.prm.h;
        voxelvol = prod(prm.h);
        dim = size(im);
        cellbw = D.cellbw;
        wat = double(D.wat);
        % DiD in first channel, transfer in last channel. If only one
        % channel, they are both equal
        imtransf = im(:,:,:,trch);


        % load the donor cell definition
        pathload = [folder{i} '/' 'stack' int2str(j) '-donorcell.mat'];                            
        try            
            D = load(pathload);
            msg = ['Loaded ' pathload];
            disp(msg);
        catch         
            msg = ['Could not load ' pathload ', continuing'];
            disp(msg);
            continue;
        end;
        % the union of the donor cells
        donorcell = D.donorcellunion;
        donorcellwat = D.donorcellwat;
        % is it a control cell?
        control = D.control;        
        
        % remove from border
        cellbw = imclearborderth(cellbw,0.05);            
        
        % remove small cells        
        conn = 18;
        % 100 correspondts to 500 voxels
        tv = 100;
        t = tv/voxelvol;        
        cellbw = bwareaopenrange(cellbw,t,conn);
        
        % remove donorcell from the watershed definition (watershed is
        % always equal or smaller than the donorcell definition)
        cellbw(donorcell) = 0;
        wat(donorcell) = 0;
                
        % make cellwat
        cellwat = cellbw .* wat; 
        
        % did segmentation should always be in number one!        
%         for k = 1 : numel(trth)
        bwtransf = gt(imtransf,trth);
        bwtransfini = bwtransf;            
        

%         showall(imtransf,bwtransf,donorcell,donorcellwat)        

        % remove signal connected to the donor cell from the transfer image
        [faser,~] = bwlabeln(bwtransf);
        val = faser(donorcell == 1);
        % it is not a control image
        if ~isempty(val)
            [val,maxval] = mostfrequent(val);
            bwtransf(faser == val) = 0;
            % remove the cells from the binary image to do detection of signal
            bwtransf(donorcell == 1) = 0;                   
        end;
                
                
%         showall(imtransf,bwtransf,donorcell)
        
                
        % cell indicator
        cellval = unique(cellwat(:));
        cellval(cellval == 0) = [];
        ncells = numel(cellval);

%         showall(imsegm,imtransf,donorcell,bwtransf,cellbw)
        
        % 
        % Find transferred signal etc in segmented cells
        %
                
        clear cells;
        cells.volsignal = NaN(ncells,1);
        cells.pos = NaN(ncells,1);
        cells.intsignal = NaN(ncells,1);       
        
        % threshold for a signal to be measured for a cell to be positive
        thpos = 5;
        for k = 1 : ncells
            % this cell
            cellhere = cellwat == cellval(k);
                                       
            % NB if you add a term here also change the sorting below!!            
                        
            % volume of signals inside cell
            reg = logical(bwtransf .* cellhere);            
            cells.volsignal(k,1) = sum(reg(:));
                         
            % sum of signal inside cell   
            cells.intsignal(k,1) = sum(imtransf(reg));
            
            % label cell as positive if it has at least 5 spots
            [faserhere,L] = bwlabeln(reg);
            cells.pos(k,1) = L > thpos;
                                    
        end;%k

        %
        % Donor cell analysis
        %        
        
        % find signal in the donorcell to print to textfile
        reg = logical(donorcell.*bwtransfini);                
        donor.volsignal = sum(reg(:));
        donor.intsignal = sum(imtransf(reg));
        
        % number of donorcells
        [~,ndonorcell] = bwlabeln(donorcell);
        
        % surface area
        perim = bwperim(donorcell);
        donor.surfsignal = sum(perim(:));        
        
%         showall(imtransf,donorcell,bwtransfini,perim)
        
        % find total (transferred?) signal outside donor cell        
        tottr.volsignal = sum(bwtransf(:));
        tottr.intsignal = sum(imtransf(bwtransf));                        
        
        % find row relations
        [cells.rownum] = cell2cellcontact(donorcellwat,cellwat,cellval,ncells);                

        %
        % Search for colocalization. 
        % NB this is only done with trmode = 'gfp' since we than work on two different images
        %
        if isequal(trmode,'gfp')
            radhere = [4 4 1];
            im1 = im(:,:,:,prm.did.ch);
            bw1global = im1 > prm.did.th;
            
            im2 = im(:,:,:,prm.gfp.ch);
            bw2global = im2 > prm.gfp.th;
                        
%             showall(im1,im2,bw1global,bw2global)
            
            bw1 = adaptfiltim(im1,radhere,0.10,prm.h);
            bw1 = bw1 .* bw1global .* cellbw;
            
            bw2 = adaptfiltim(im2,radhere,0.10,prm.h);
            bw2 = bw2 .* bw2global .* cellbw;
            
            
            [coloc.voloverlap1,coloc.voloverlap2,coloc.vol1,coloc.vol2] = ...
                colocalize(bw1,bw2,0.5,1);

            coloc.n1 = numel(coloc.voloverlap1);
            coloc.n2 = numel(coloc.voloverlap2);
            % percentage colocalization
            coloc.perc1 = sum(coloc.voloverlap1 > 0) / coloc.n1;
            coloc.perc2 = sum(coloc.voloverlap2 > 0) / coloc.n2;
            coloc.meanvol1 = mean(coloc.vol1);
            coloc.meanvol2 = mean(coloc.vol2);
            
        else
            coloc.voloverlap1 = NaN;
            coloc.voloverlap2 = NaN;
            coloc.vol1 = NaN;
            coloc.vol2 = NaN;            
            coloc.n1 = NaN;
            coloc.n2 = NaN;
            coloc.perc1 = NaN;
            coloc.perc2 = NaN;            
            coloc.meanvol1 = NaN;
            coloc.meanvol2 = NaN;
        end;

        %
        % Find the longest distance in the didcell
        %        
        maxprojdonorcell = maxprojimage(donorcell);
        ind = find(maxprojdonorcell);
        clear c;
        [c(:,1), c(:,2)] = ind2sub(dim,ind);
        n = numel(ind);
        dist = 0;
        for k = 1 : n
            disthere = bsxfun(@minus,c(k:end,:),c(k,:));
            disthere = sqrt(sum(disthere.^2,2));
            disthere = max(disthere);
            dist = max(dist,disthere);
        end            
        donor.maxdist = dist;
        
%         showall(imsegm,bw1,imtransf,bw2,cellbw)

        
        % label the test image with the row number of the cells
        % for evaluation
        val = [1,1.5];
        testim = zeros(dim(1:3));                
        for k = 1 : ncells
            cellhere = cellwat == cellval(k);
            testim(cellhere) = val(cells.rownum(k));
        end;                       
%         a1 = scale(imtransf)*2;
%         a2 = scale(imtransf)*2;
        perim = bwperim(donorcell);
%         load ball1;se = getball(ball,1,1);
%         perim = imdilate(perim,se);
        testim(donorcell == 1) = 0.75;
        testim(perim == 1) = 0.4;
        testim(bwtransfini == 1) = 2.0;

%         % for the paper: 3-1_2c-manual/XWNeg9_s444246_B1
%         val = [1,1.5];
%         testim = zeros(dim(1:3));                
%         for k = 1 : ncells
%             cellhere = cellwat == cellval(k);
%             testim(cellhere) = 1;
%         end;                       
%         perim = bwperim(donorcell);
%         testim(donorcell == 1) = 0.75;
%         testim(donorcell == 1) = 2;
%         testim(perim == 1) = 0.4;
%         testim(bwtransf == 1) = 2;

        
        
        % print images to disc for easy evaluation
        filesave = [foldereval2 '/' 'im' int2str(j) '-eval-' trmode '-th-' int2str(trth) '.tif'];
        msg = ['Saving ' filesave];
        disp(msg);
        mat2tifdirect(testim,filesave);
        filesave = [foldereval1 '/' 'im' int2str(j) '-eval-oneplane-' trmode '-th-' int2str(trth) '.tif'];
        msg = ['Saving ' filesave];
        disp(msg);
        imwrite(prepimwrite(testim(:,:,prm.planetrack)),filesave,'tiff');

        filesave = [foldereval2 '/' 'im' int2str(j) '-imsegm-' trmode '-th-' int2str(trth) '.tif'];
        msg = ['Saving ' filesave];
        disp(msg);
        mat2tifdirect(imtransf,filesave);
        filesave = [foldereval1 '/' 'im' int2str(j) '-imsegm-oneplane-' trmode '-th-' int2str(trth) '.tif'];
        msg = ['Saving ' filesave];
        disp(msg);
        imwrite(prepimwrite(imtransf(:,:,prm.planetrack)),filesave,'tiff');
        
%         ndonorcell
%         showall(imtransf,bwtransf,imsegm,bwsegm,donorcell,donorcellwat)
%         showall(D.imsegm,D.cellbw,cellwat,D.minima,imsegm,donorcellwat)
        
        if vis == 1
            trch
            ndonorcell
            figure(1);
            subplot(2,2,1);imagesc(imsegm(:,:,prm.planetrack));colormap(gray);axis image;axis off;title('Segmentaiton image');
            subplot(2,2,2);imagesc(imtransf(:,:,prm.planetrack));colormap(gray);axis image;axis off;title('Transfer image');
            subplot(2,2,3);imagesc(cellwat(:,:,prm.planetrack));colormap(gray);axis image;axis off;title('Cell watershed');
            subplot(2,2,4);imagesc(testim(:,:,prm.planetrack));colormap(gray);axis image;axis off;title('Testim');            
            pause
        end;
        

        % NB sort the array according row number!!
        [cells.rownum,indsort] = sort(cells.rownum);        
        cells.volsignal = cells.volsignal(indsort);
        cells.intsignal = cells.intsignal(indsort);
        cells.pos = cells.pos(indsort);
        

        
        % store information in struct to print to textfile        
        stats.folder(i).im(j).num = j;
        stats.folder(i).im(j).donor = donor;
        stats.folder(i).im(j).cells = cells;
       
        
        n = numel(cells.rownum);
        stats.folder(i).im(j).n = n;
        stats.folder(i).im(j).ndonorcell = ndonorcell;
        stats.folder(i).im(j).tottr = tottr;
        stats.folder(i).im(j).ind1 = cells.rownum == 1;
        stats.folder(i).im(j).ind2 = cells.rownum == 2;
        stats.folder(i).im(j).ind = cells.rownum > -Inf;
        stats.folder(i).im(j).coloc = coloc;
        stats.folder(i).im(j).control = control;
        
    end;
    
end;

% rearrange data for printing
try
    data = rearrangedata(stats);
catch
    save testanalysetransfer1
    error('error')
end;

% print summary
base = ['summary' filesavefinal '-' trmode '-' 'th' num2str(trth)];
if isempty(foldersavefinal)
    pathsave = [base  '.txt'];
else
    pathsave = [foldersavefinal '/' base  '.txt'];
end;
msg = ['Printing ' pathsave];
disp(msg);
try
    printsummary(data,pathsave);
catch
    save testanalysetransfer2
    error('error')
end;


if isempty(foldersavefinal)
    pathsave = [base  '.mat'];
else
    pathsave = [foldersavefinal '/' base  '.mat'];
end;
msg = ['Printing ' pathsave];
disp(msg);
save(pathsave,'stats','data');

%------------------------------------------------

function [] = printsummary(data,pathsave)


fid = fopen(pathsave,'wt');
p = 40;

A{1}  =['Summary from ' mfilename];
printcell(fid,A,p);
fprintf(fid,'\n');
description{1,1} = 'int-donor: The total signal inside the donor cell, above threshold';
description{2,1} = 'vol-donor: Volume of donor cell, above threshold';
description{3,1} = 'maxdist-did: Maximum distance within the donor cell';
description{4,1} = 'vol-signal-row1: Number of pixels in the first row cells, above threshold';
description{5,1} = 'int-signal-row1: Total intensity in the first row cells, above threshold';
description{6,1} = 'pos-row1: Number of positive cells (with number of spots, above threshold)';
description{7,1} = 'tot-transf-vol: Total volume transferred signal from donor cell, above threshold';
description{8,1} = 'tot-transf-int: Total intensities transferred signal from donor cell, above threshold';
description{9,1} = 'NB: row3 are ALL cells';
description{10,1} = 'surf-donor: Surface volume of donor cell';
description{11,1} = 'Removing all data points if there are no cells found in either first or second row, in both control or noncontrol';


printcell(fid,description,p);
fprintf(fid,'\n');
% varnametext = {'DiD-cell-int','DiD-cell-vol','DiD-max-dist','1st-row-vol','1st-row-int','1st-row-pos(>5)','1st-row-n-cells','2nd-row-vol','2nd-row-int','2nd-row-pos(>5)','2nd-row-n-cell','1and2-row-vol','1and2-row-pos(>5)','1and2-row-n-cells','tot-transf-vol','tot-transf-int'};


%
% Print cell characteristics
%

msg{1} = upper('Cell row characteristics');
printcell(fid,msg,p);

nrow = 3;
Dc = [];
Dnc = [];
nfolder = numel(data.folder);
for i = 1 : data.nvarnamecells

    varname = data.varnamecells{i};
    c = 0;  
    for j = 1 : nrow
    
        Bc = NaN(data.nimmax,nfolder);
        Bnc = NaN(data.nimmax,nfolder);
        for k = 1 : nfolder                
            c = c + 1;
            % NB can not only nansum since nansum of only nan is 0!!! And
            % then problems when its written 0 instead of NaN
            
            % control                        
            varhere = data.folder(k).control.cells.(varname).row{j};  
            var = nansum(varhere,1);
            n1 = size(varhere,1);
            label = isnan(varhere);
            label = sum(label,1);
            label = label == n1;
            var(label) = NaN;
            varhere = var;            
            varhere = varhere(:); 
            nim = numel(varhere);
            Bc(1:nim,k) = varhere;            
            
            % noncontrol
            varhere = data.folder(k).noncontrol.cells.(varname).row{j};                     
            var = nansum(varhere,1);
            n1 = size(varhere,1);
            label = isnan(varhere);
            label = sum(label,1);
            label = label == n1;
            var(label) = NaN;
            varhere = var;
            varhere = varhere(:);                        
            nim = numel(varhere);
            Bnc(1:nim,k) = varhere; 

        end;          
        
        
        % collect the mean values over the images
        meanval = nanmean(Bc,1);
        Dc = [Dc meanval'];
        meanval = nanmean(Bnc,1);
        Dnc = [Dnc meanval'];

        % print to file
        
        % control
        Bc = mat2celldirect(Bc);
        clear header;
        header{1} = [varname '-row' int2str(j) '-' 'control'];        
        printcell(fid,header,p);
        header = [{'Image'} data.foldername'];
        printcell(fid,header,p);
        left = mat2celldirect((1:size(Bc,1))');
        Bc = [left Bc];
        printcell(fid,Bc,p);
        fprintf(fid,'\n');
        
        % non-control
        Bnc = mat2celldirect(Bnc);
        clear header;
        header{1} = [varname '-row' int2str(j) '-' 'noncontrol'];        
        printcell(fid,header,p);
        header = [{'Image'} data.foldername'];
        printcell(fid,header,p);
        left = mat2celldirect((1:size(Bnc,1))');
        Bnc = [left Bnc];
        printcell(fid,Bnc,p);
        fprintf(fid,'\n');

    end;        
    
end;

%
% Print mean values
%

msg{1} = upper('Mean cell row characteristics');
printcell(fid,msg,p);

c = 0;
clear headernc headerc;
for i = 1 : data.nvarnamecells
    varname = data.varnamecells{i};
    for j = 1 : nrow
        c = c + 1;
        headernc{1,c} = [varname '-row' int2str(j) '-' 'noncontrol'];        
        headerc{1,c} = [varname '-row' int2str(j) '-' 'control'];        
    end;
end;
headernc = [{'Folder'} headernc];
headerc = [{'Folder'} headerc];

% control
printcell(fid,headerc,p);
Dc = mat2celldirect(Dc);
Dc = [data.foldername Dc];
printcell(fid,Dc,p);
fprintf(fid,'\n');

% non control
printcell(fid,headernc,p);
Dnc = mat2celldirect(Dnc);
Dnc = [data.foldername Dnc];
printcell(fid,Dnc,p);
fprintf(fid,'\n');

%
% Print total transfer
%

msg{1} = upper('Total and background transfer');
printcell(fid,msg,p);

for i = 1 : data.nvarnametottr
    varname = data.varnametottr{i};
    Bc = NaN(data.nimmax,nfolder);
    Bnc = NaN(data.nimmax,nfolder);
    Cc = NaN(data.nimmax,nfolder);
    Cnc = NaN(data.nimmax,nfolder);
    
    for j = 1 : nfolder
        var = data.folder(j).control.tottr.(varname);
        var = var(:);
        n = numel(var);
        Bc(1:n,j) = var;
        
        var = data.folder(j).noncontrol.tottr.(varname);
        var = var(:);
        n = numel(var);
        Bnc(1:n,j) = var;
        
        var = data.folder(j).control.bcktr.(varname);
        var = var(:);
        n = numel(var);
        Cc(1:n,j) = var;
        
        var = data.folder(j).noncontrol.bcktr.(varname);
        var = var(:);
        n = numel(var);
        Cnc(1:n,j) = var;
    end;

    n = (1:size(Bc,1))'; 
    clear header;
    header{1} = ['tottransf' '-' varname '-' 'control'];
    printcell(fid,header,p);
    header = [{'Image'} data.foldername'];
    printcell(fid,header,p);
    Bc = [mat2celldirect(n) mat2celldirect(Bc)];
    printcell(fid,Bc,p);
    fprintf(fid,'\n');

    clear header;
    header{1} = ['tottransf'  '-' varname '-' 'noncontrol'];
    printcell(fid,header,p);
    header = [{'Image'} data.foldername'];
    printcell(fid,header,p);
    Bnc = [mat2celldirect(n) mat2celldirect(Bnc)];
    printcell(fid,Bnc,p);
    fprintf(fid,'\n');
   
    clear header;
    header{1} = ['bcktransf' '-' varname '-' 'control'];
    printcell(fid,header,p);
    header = [{'Image'} data.foldername'];
    printcell(fid,header,p);
    Cc = [mat2celldirect(n) mat2celldirect(Cc)];
    printcell(fid,Cc,p);
    fprintf(fid,'\n');

    clear header;
    header{1} = ['bcktransf'  '-' varname '-' 'noncontrol'];
    printcell(fid,header,p);
    header = [{'Image'} data.foldername'];
    printcell(fid,header,p);
    Cnc = [mat2celldirect(n) mat2celldirect(Cnc)];
    printcell(fid,Cnc,p);
    fprintf(fid,'\n');

end;


%
% Print donor cell characteristics
%

msg{1} = upper('Donor-cell-characteristics');
printcell(fid,msg,p);

for i = 1 : data.nvarnamedonor
    
    varname = data.varnamedonor{i};    
           
    Bc = NaN(data.nimmax,nfolder);
    Bnc = NaN(data.nimmax,nfolder);
    for j = 1 : nfolder                
        varhere = data.folder(j).control.donor.(varname);        
        varhere = varhere(:); 
        n = numel(varhere);
        Bc(1:n,j) = varhere;

        varhere = data.folder(j).noncontrol.donor.(varname);        
        varhere = varhere(:);                        
        n = numel(varhere);
        Bnc(1:n,j) = varhere;

    end;                       
    Bc = mat2celldirect(Bc);
    Bnc = mat2celldirect(Bnc);

    
    clear header;
    header{1} = [varname '-' 'control'];
    printcell(fid,header,p);
    clear header;
    header = [{'Image'} data.foldername'];
    printcell(fid,header,p);
    Bc = [mat2celldirect((1:size(Bc,1))') Bc];
    printcell(fid,Bc,p);
    fprintf(fid,'\n');

    clear header;
    header{1} = [varname '-' 'noncontrol'];
    printcell(fid,header,p);
    clear header;
    header = [{'Image'} data.foldername'];
    printcell(fid,header,p);
    Bnc = [mat2celldirect((1:size(Bnc,1))') Bnc];
    printcell(fid,Bnc,p);
    fprintf(fid,'\n');

    
end;

%
% Print colocalization characteristics
%
msg{1} = upper('Colocalization-characteristics');
printcell(fid,msg,p);

usevarname = {'perc1','perc2','n1','n2','meanvol1','meanvol2'};
for i = 1 : numel(usevarname)
    varname = usevarname{i};
    Bc = NaN(data.nimmax,nfolder);
    Bnc = NaN(data.nimmax,nfolder);
    for j = 1 : nfolder
        varhere = data.folder(j).control.coloc.(varname);    
        n = numel(varhere);
        Bc(1:n,j) = varhere;

        varhere = data.folder(j).noncontrol.coloc.(varname);    
        n = numel(varhere);
        Bnc(1:n,j) = varhere;
    end;
    Bc = mat2celldirect(Bc);
    Bnc = mat2celldirect(Bnc);
    
    clear header;
    header{1} = [varname '-' 'control'];
    printcell(fid,header,p);
    clear header;
    header = [{'Image'} data.foldername'];
    printcell(fid,header,p);
    Bc = [mat2celldirect((1:size(Bc,1))') Bc];
    printcell(fid,Bc,p);
    fprintf(fid,'\n');

    clear header;
    header{1} = [varname '-' 'noncontrol'];
    printcell(fid,header,p);
    clear header;
    header = [{'Image'} data.foldername'];
    printcell(fid,header,p);
    Bnc = [mat2celldirect((1:size(Bnc,1))') Bnc];
    printcell(fid,Bnc,p);
    fprintf(fid,'\n');
    
end;

fclose(fid);

%------------------------------------------------

function [data] = rearrangedata(stats)


data.varnamecells = {'intsignal','volsignal','pos','ncells'};
data.nvarnamecells = numel(data.varnamecells);
data.varnamedonor = {'intsignal','volsignal','maxdist','surfsignal'};
data.nvarnamedonor = numel(data.varnamedonor);
data.varnametottr = {'intsignal','volsignal'};
data.nvarnametottr = numel(data.varnametottr);
data.varnamecoloc = {'n1','n2','perc1','perc2','meanvol1','meanvol2'};
data.nvarnamecoloc = numel(data.varnamecoloc);

nfolder = numel(stats.folder);
nrow = 3;
for i = 1 : nfolder
    
    if isempty(stats.folder(i).im)
        data.folder(i).nim = 0;
        continue;
    end;
    data.folder(i).name = stats.folder(i).name;
    nim = numel(stats.folder(i).im);
    ncellsmax = 0;
    for j = 1 : nim
        if isempty(stats.folder(i).im(j).n)
            stats.folder(i).im(j).n = 0;
        end;
        ncellsmax = max(ncellsmax,stats.folder(i).im(j).n);   
    end;

    for j = 1 : nrow
        for k = 1 : data.nvarnamecells
            varname = data.varnamecells{k};
            data.folder(i).control.cells.(varname).row{j} = NaN(ncellsmax,nim);
            data.folder(i).noncontrol.cells.(varname).row{j} = NaN(ncellsmax,nim);
        end;
    end;
    
    for j = 1 : data.nvarnamedonor
        varname = data.varnamedonor{j};
        data.folder(i).control.donor.(varname) = NaN(1,nim);
        data.folder(i).noncontrol.donor.(varname) = NaN(1,nim);
    end;
    for j = 1 : data.nvarnametottr
        varname = data.varnametottr{j};
        data.folder(i).control.tottr.(varname) = NaN(1,nim);    
        data.folder(i).noncontrol.tottr.(varname) = NaN(1,nim);
        data.folder(i).control.bcktr.(varname) = NaN(1,nim);    
        data.folder(i).noncontrol.bcktr.(varname) = NaN(1,nim);        
    end;
    % the 'voloverlap' are cell arrays of colocalization, therefore we start at 3
    for j = 1 : data.nvarnamecoloc
        varname = data.varnamecoloc{j};
        data.folder(i).control.coloc.(varname) = NaN(1,nim);        
        data.folder(i).noncontrol.coloc.(varname) = NaN(1,nim);    
    end;

    varname = 'voloverlap1';
    data.folder(i).control.coloc.(varname) = cell(1,nim);        
    data.folder(i).noncontrol.coloc.(varname) = cell(1,nim);    
    varname = 'voloverlap2';
    data.folder(i).control.coloc.(varname) = cell(1,nim);        
    data.folder(i).noncontrol.coloc.(varname) = cell(1,nim);    
    varname = 'vol1';
    data.folder(i).control.coloc.(varname) = cell(1,nim);        
    data.folder(i).noncontrol.coloc.(varname) = cell(1,nim);    
    varname = 'vol2';
    data.folder(i).control.coloc.(varname) = cell(1,nim);        
    data.folder(i).noncontrol.coloc.(varname) = cell(1,nim);    
    
    
    clear ind
    for j = 1 : nim
        
        % the rows!
        ind{1} = logical(stats.folder(i).im(j).ind1);
        ind{2} = logical(stats.folder(i).im(j).ind2);
        ind{3} = logical(stats.folder(i).im(j).ind);        

        % control?
        control = stats.folder(i).im(j).control;
        c = 'noncontrol';
        if control == 1
            c = 'control';
        end;
        
        % no cell data, why?
        if isempty(stats.folder(i).im(j).cells)
            continue;
        end;
        
        for k = 1 : nrow            
            try
                data.folder(i).(c).cells.intsignal.row{k}(ind{k},j) = stats.folder(i).im(j).cells.intsignal(ind{k});
                data.folder(i).(c).cells.volsignal.row{k}(ind{k},j) = stats.folder(i).im(j).cells.volsignal(ind{k});                
                data.folder(i).(c).cells.pos.row{k}(ind{k},j) = stats.folder(i).im(j).cells.pos(ind{k});                                    
                data.folder(i).(c).cells.ncells.row{k}(1,j) = sum(ind{k});
            catch
                data.folder(i).(c).cells.ncells.row{k}(1,j) = 0;
            end;                       
        end;
                
        data.folder(i).(c).donor.intsignal(1,j) = stats.folder(i).im(j).donor.intsignal;
        data.folder(i).(c).donor.volsignal(1,j) = stats.folder(i).im(j).donor.volsignal;
        data.folder(i).(c).donor.surfsignal(1,j) = stats.folder(i).im(j).donor.surfsignal;
        data.folder(i).(c).donor.maxdist(1,j) = stats.folder(i).im(j).donor.maxdist;

        
        % total transfer
        data.folder(i).(c).tottr.intsignal(1,j) = stats.folder(i).im(j).tottr.intsignal;
        data.folder(i).(c).tottr.volsignal(1,j) = stats.folder(i).im(j).tottr.volsignal;

        % background transfer is total transfer minus row3 transfer (all
        % cells)
        data.folder(i).(c).bcktr.intsignal(1,j) = stats.folder(i).im(j).tottr.intsignal - nansum(data.folder(i).(c).cells.intsignal.row{3}(:,j));
        data.folder(i).(c).bcktr.volsignal(1,j) = stats.folder(i).im(j).tottr.volsignal - nansum(data.folder(i).(c).cells.volsignal.row{3}(:,j));
        
        % colocalization        
        data.folder(i).(c).coloc.n1(1,j) = stats.folder(i).im(j).coloc.n1;
        data.folder(i).(c).coloc.n2(1,j) = stats.folder(i).im(j).coloc.n2;
        data.folder(i).(c).coloc.perc1(1,j) = stats.folder(i).im(j).coloc.perc1;
        data.folder(i).(c).coloc.perc2(1,j) = stats.folder(i).im(j).coloc.perc2;
        data.folder(i).(c).coloc.meanvol1(1,j) = stats.folder(i).im(j).coloc.meanvol1;
        data.folder(i).(c).coloc.meanvol2(1,j) = stats.folder(i).im(j).coloc.meanvol2;
        data.folder(i).(c).coloc.voloverlap1{1,j} = stats.folder(i).im(j).coloc.voloverlap1;
        data.folder(i).(c).coloc.voloverlap2{1,j} = stats.folder(i).im(j).coloc.voloverlap2;
        data.folder(i).(c).coloc.vol1{1,j} = stats.folder(i).im(j).coloc.vol1;
        data.folder(i).(c).coloc.vol2{1,j} = stats.folder(i).im(j).coloc.vol2;
        
    end;
    data.folder(i).nim = nim;
end;
data.nim = 0;
data.nimmax = 0;
for i = 1 : nfolder
    data.foldername{i,1} = stats.folder(i).name;
    data.nim = data.nim + data.folder(i).nim;
    data.nimmax = max(data.nimmax,data.folder(i).nim);
end;

% Remove all data points if there are no (< 4) cells found in either first or
% second row, in neither control or noncontrol
rem = zeros(nfolder,1);
for i = 1 : nfolder

    if data.folder(i).nim == 0
        rem(i) = 1;
        continue;
    end;

    
    % remove from ncells since its only a scalar
    for j = 1 : 3
        data.folder(i).control.cells.ncells.row{j} = data.folder(i).control.cells.ncells.row{j}(1,:);
        data.folder(i).noncontrol.cells.ncells.row{j} = data.folder(i).noncontrol.cells.ncells.row{j}(1,:);
    end;

    % no cells in one of the rows? Remove the data since the data are
    % considered corrupt
    ind = data.folder(i).control.cells.ncells.row{1}(1,:) < 4 | data.folder(i).control.cells.ncells.row{2}(1,:) == 0; 
    for j = 1 : data.nvarnamecells
        for k = 1 : 3
            data.folder(i).control.cells.(data.varnamecells{j}).row{k}(:,ind) = NaN;
        end;
    end;
    for j = 1 : data.nvarnamedonor
        varname = data.varnamedonor{j};
        data.folder(i).control.donor.(varname)(ind) = NaN;        
    end;

    
    ind = data.folder(i).noncontrol.cells.ncells.row{1}(1,:) < 4 | data.folder(i).noncontrol.cells.ncells.row{2}(1,:) == 0; 
    for j = 1 : data.nvarnamecells
        for k = 1 : 3
            data.folder(i).noncontrol.cells.(data.varnamecells{j}).row{k}(:,ind) = NaN;
        end;
    end;
    for j = 1 : data.nvarnamedonor
        varname = data.varnamedonor{j};
        data.folder(i).noncontrol.donor.(varname)(ind) = NaN;        
    end;


end;

% remove empty folders
data.folder(rem == 1) = [];
data.foldername(rem == 1) = [];


%-----------------------------------------------

function [valB,valC] = getvarcells(stats,varhere,inds)

nfolder = numel(stats.folder);
for j = 1 : nfolder

    nim = numel(stats.folder(j).im);
    valB{j} = [];
    valC{j} = [];
    for k = 1 : nim
        imhere = stats.folder(j).im(k);        
        ind = imhere.(inds);        
        % non-control cells
        if imhere.ndonorcell >= 1
            valB{j} = [valB{j}; sum(imhere.cells.(varhere)(ind))];
        % control cells
        else
            valC{j} = [valC{j}; sum(imhere.cells.(varhere)(ind))];
        end;                
            
    end;
            
end;

% convert to matrix
valB = cell2mat(valB);
valC = cell2mat(valC);

%-----------------------------------------

function [valB,valC] = getvardid(stats,varhere)

nfolder = numel(stats.folder);
for j = 1 : nfolder

    nim = numel(stats.folder(j).im);
    valB{j} = [];
    valC{j} = [];    
    for k = 1 : nim
        imhere = stats.folder(j).im(k);                
        if imhere.ndonorcell >= 1
            valB{j} = [valB{j}; sum(imhere.did.(varhere))];
        else
            valC{j} = [valC{j}; sum(imhere.did.(varhere))];
        end;                

            
    end;
            
end;
% convert to matrix
valB = cell2mat(valB);
valC = cell2mat(valC);




%------------------------------------------

function [rownum] = cell2cellcontact(didbw,cellwat,cellval,ncells)

% was at this
% load ball3;se = getball(ball,3,1);
% change 20120810 after discussion with Tanja
% load ball22;se = getball(ball,22,1);
% dildidbw = imdilate(didbw,se);

% find average cell diamater to measure neighbor relationships
properties = {'MajorAxisLength','MinorAxisLength'};
avcelldiam = NaN(ncells,1);
for i = 1 : ncells
    bw = cellwat == cellval(i);
    bw = bw(:,:,4);
    if isempty(find(bw,1))
        continue;
    end;
    props = regionprops(double(bw),properties);
    avcelldiam(i) = (props.MajorAxisLength + props.MinorAxisLength)/2;
end;
avcelldiam = nanmean(avcelldiam);

% find distance from the DiD cell
d = bwdist(didbw);
dildidbw = d < 0.75*avcelldiam;
dildidbw(didbw == 1) = 0;

% find the neighbours
val = cellwat(dildidbw == 1);
val(val == 0) = [];
val = unique(val);
[inters,inda,indb] = intersect(val,cellval);

% 1 is direct neighbour and 2 is distant neighbour
% default is 2, distant neighbour
rownum = 2*ones(ncells,1);
rownum(indb) = 1;

% val
% 
% val
% cellval
% rownum
% showall(didbw,dildidbw,cellwat,cellwat == 2)

%------------------------------------------




