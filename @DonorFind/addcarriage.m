function [cmd] = addcarriage(cmd)
%  ADDCARRIAGE Adding a carriage and new line to the command, otherwise
%  its not correctly read by MatrixScreener
%
%   CMD = ADDCARRIAGE(CMD) adding carriage and new line to the command. 
%   Returning the modified command. 
%
%     =======================================================================================
%     Copyright (C) 2014  Erlend Hodneland
%     Email: erlend.hodneland@gmail.com 
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%     =======================================================================================
%


cmd(end) = [];
cmd = sprintf('%s\r\n',cmd);