function [] = sendcmd(obj,cmd)
%
% SENDCMD Sending a command to the CAM server
%
%   SENDCMD(OBJ,CMD) sends a command to the CAM server defined by the 
%   TCP/IP object OBJ as defined by OPENCONNECTION. The command is
%   given as a string in CMD.
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

pause(0.5);
msg = ['Sending ' cmd];
disp(msg);
fprintf(obj,cmd);
