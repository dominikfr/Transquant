function [lowres] = getimageinfo(mes)
%
% GETIMAGEINFO Gets image information from the string that is received from
% MAtrixScreener
%
%   LOWRES = GETIMAGEINFO(MES) takes a string coming from MatrixScreener
%   and returns a struct with variables, both as string and numeric, if
%   possible
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


% % remove backslash and replace by slash
% ind.slash = strfind(mes,'\');
% for i = 1 : numel(ind.slash)
%     mes(ind.slash(i)) = '/';
% end;


% find the path
str = 'relpath:';
ind.relpath = strfind(mes,str);
d = numel(str);
lowres.relpath = mes(ind.relpath+d:end);
ind2 = strfind(lowres.relpath,'/jobname');
lowres.relpath = lowres.relpath(1:ind2-1);
[lowres.folder, lowres.filename, lowres.ext] = fileparts(lowres.relpath);

% find the jobname
ind = strfind(mes,'/jobname');
lowres.jobnameread = deblank(mes(ind+9:end));

% find the attributes
str = {'L','S','U','V','J','E','O','X','Y','T','Z','C'};
[att,attn] = getcmdpiece(lowres.filename,str);
lowres.attribute = att;
lowres.attributen = attn;
lowres.attributename = str;

%-------------------------------------------------------------

function [piece,piecen] = getcmdpiece(mes,str)

for i = 1 : numel(str)
    ind = strfind(mes,str{i});
    n = numel(str{i});
    a = mes(ind+n:end);
    
    ind2 = strfind(a,'--');
    if isempty(ind2)
        ind2 = strfind(a,'.');
    end    
    b = a(1:ind2-1);
    piece.(str{i}) = b;
    try
        piecen.(str{i}) = str2double(b);
    catch
    end;
 
end;

%-----------------------
