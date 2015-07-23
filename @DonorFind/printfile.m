function [] = printfile(A,name,prm,i)
%
% PRINTFILE Printing tif files of the super images.
%
%   PRINTFILE(A,NAME,PRM,I) prints a tif image of the super image in A,
%   with the specifica name NAME of line I of the screening template. The
%   PRM struct is defined in the SCREEN script.
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

A = uint8(A);
pathsave = fullfile(prm.foldersave,['line' int2str(prm.lowres.line(i)) '-' name '.tiff']);    
msg = ['Printing ' pathsave];
disp(msg);
imwrite(A,pathsave,'tiff');
