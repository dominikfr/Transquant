%    CELLSEGM Class definition file
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
classdef DonorFind
       
   methods (Static)
      % main functions
      [cmd] = addcarriage(cmd);
      [v] = addcmd(target,value);
      [cmd] = addtocamlist(lowres,jobname,cm);
      [cmcontrol,controlim] = controlroi(imsegm,labelim,nuclei,rim,prm);
      [cmd] = deletelist;
      [lowres] = getimageinfo(mes);
      [imsegm,imnucleus,labelim,did,cm,nuclei,rim,controlim] = getimageobj(plate,prm);
      [plate] = getlowres(prm,obj);
      [obj] = openconnection;
      printfile(A,name,prm,i);
      screen;
      sendcmd(obj,cmd);
      [cmd] = startcamscan;
      [cmdsave] = starthighres(obj,cm,cmlabelim,info,labelim,jobname);
      [cmd] = stopwaitcamscan;
      
   end 
      
end