function [cmd] = addtocamlist(lowres,jobname,cm)
%
% ADDTOCAMLIST Adding a highre image job to the CAM list.
%
%   CMD = ADDTOCAMLIST(LOWRES,JOBNAME,CM) adds the list of coordinate in CM
%   to the camlist of highres jobs. The LOWRES struct is a from GETIMAGEOBJ
%   and JOBNAME is the highres jobname.
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
%
cmd = [matrixcompute.addcmd('cli','administrator') ...
    matrixcompute.addcmd('app','matrix') ...
    matrixcompute.addcmd('sys','1') ...
    matrixcompute.addcmd('cmd','add')  ...
    matrixcompute.addcmd('tar','camlist') ...
    matrixcompute.addcmd('exp',jobname) ...
    matrixcompute.addcmd('slide',num2str(lowres.attributen.S)) ...
    matrixcompute.addcmd('wellx',num2str(lowres.attributen.U)) ...
    matrixcompute.addcmd('welly',num2str(lowres.attributen.V)) ...
    matrixcompute.addcmd('fieldx',num2str(lowres.attributen.X+1)) ...
    matrixcompute.addcmd('fieldy',num2str(lowres.attributen.Y+1)) ...
    matrixcompute.addcmd('dxpos',num2str(cm(2))) ...
    matrixcompute.addcmd('dypos',num2str(cm(1))) ...    
    ];

cmd = matrixcompute.addcarriage(cmd);