function [cmd] = deletelist
%
% DELETELIST Deleting list for highres imaging in the MATRIXSCREENER
%
%   CMD = DELETELIST sends the command for deleting the list for highres
%   screening. Should always be done prior to start of highres to ensure
%   the list is empty.
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
cmd = [matrixcompute.addcmd('cli','administrator') ...
    matrixcompute.addcmd('app','matrix') ...
    matrixcompute.addcmd('cmd','deletelist') ...    
    ];

cmd = matrixcompute.addcarriage(cmd);

