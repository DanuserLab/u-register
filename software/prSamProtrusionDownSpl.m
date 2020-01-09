function [output] = prSamProtrusionDownSpl(inputParam, dl_rate)
% prSamProtrusionDownSpl is funciton to perform new Protrusion model with downsampling based
% on constrainted optimization, the input edge pixels are the cell boundary
% after image segmentation, and the mask at t-1 is used to calculating unit
% normal vectors along the edge.  The output consists of out egde pixels 
%
% Input:    inputParam.list_last        input edge pixels at t-1
%           inputParam.list             input edge pixels at t 
%           inputParam.maskTM1          mask at t-1 for calculating unit normal vectors along the edge
%
% Output:   output.pixel_tm1_output     output edge pixels at t-1 (after spline parameterization)
%           output.translate_output     translation (protrusion) from t-1 to t
%           output.pixel_t_output       output edge pixels at t (after spline parameterization)
%           output.xy_normal            unit normal vectors along the edge
%           output.time_elapsed         computational time for protrusion calculation 
%
% Last updated: August 26, 2008 by Shann-Ching Chen, LCCB
% See also: protrusionAnalysis, prSamProtrusion, prSamProtrusionOptSp
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
%
% This file is part of WindowingPackage.
% 
% WindowingPackage is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% WindowingPackage is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with WindowingPackage.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% parameter for spline smoothing
TOLERANCE = inputParam.TOLERANCE;

pixel_list_last = inputParam.pixel_list_last;
pixel_list = inputParam.pixel_list;
start_position = pixel_list_last(1,:);

ISCLOSE = inputParam.ISCLOSE;
l1 = size(pixel_list_last,1);
l2 = size(pixel_list,1);

if ISCLOSE == 1
    tmp = sum((repmat(start_position, [l2 1]) - pixel_list).^2,2);
    [t1, idx2] = min(tmp);
    pixel_list = [pixel_list(idx2:end,:); pixel_list(1:idx2-1,:)];
end

clear global D;
global D

len_last_old = size(pixel_list_last,1);
len_old = size(pixel_list,1);

PIXEL_SAMP_RATE = dl_rate;
OVERLAP_PIXEL = 6;
if ISCLOSE == 1
    pixel_list_last_extended = [pixel_list_last; pixel_list_last; pixel_list_last];
    sp_s = spaps([1:3*len_last_old],pixel_list_last_extended',TOLERANCE);
else
    sp_s = spaps([1:len_last_old],pixel_list_last',TOLERANCE);
end

spacing = len_last_old/round(len_last_old/PIXEL_SAMP_RATE);
last_ds = [1:spacing:len_last_old len_last_old];
normal_ds = 1:1:len_last_old;

if ISCLOSE == 1
    pixel_tm1=fnval(sp_s,[len_last_old + last_ds])';
    xy_normal = prSamCalculateNormal(inputParam.maskTM1, sp_s,normal_ds+len_last_old);    
else
    last_ds(end) = len_last_old;
    pixel_tm1=fnval(sp_s,[last_ds])';
    xy_normal = prSamCalculateNormal(inputParam.maskTM1, sp_s,normal_ds);        
end
len_last = size(pixel_tm1,1);
D.len_last = len_last;

if ISCLOSE == 1
    % make replicate of boundary pixel at time t
    pixel_t = [pixel_list(end-OVERLAP_PIXEL*PIXEL_SAMP_RATE+1:end,:); pixel_list(:,:); pixel_list(1:OVERLAP_PIXEL*PIXEL_SAMP_RATE,:)];
else
    pixel_t = pixel_list;
end

len = size(pixel_t,1);
D.len = len;

% initial solution
if ISCLOSE == 1
    D.sp_s = spaps([1:len]-1,pixel_t',TOLERANCE);
    d_s = ([1:len_last]-1+OVERLAP_PIXEL)*(size(pixel_list,1))/len_last;
    pixel_old = fnval(D.sp_s,d_s')';
else
    D.sp_s = spaps([1:len],pixel_t',TOLERANCE);
    spacing = (len-1)/(len_last-1);
    d_s = 1:spacing:len;
    d_s(1) = 1; d_s(end) = len;
end

D.pixel_t = pixel_t;
D.pixel_tm1 = pixel_tm1;
D.d_s = d_s;
D.iter = 0;

% A and B are the contraints for s1<s2, s2<s3 ...
A = zeros(len_last, len_last);
for i=1+1:len_last
    A(i,i-1:i) = [1 -1];
end
A(1,1) = -1;
B = zeros(len_last,1);
w = d_s';

if isempty(strfind(version,'R2008'))
    options = optimset('largescale','off','Display','off','TolCon', 1e-1,'MaxIter', 50000000);
else
    options = optimset('largescale','off','Display','off','TolCon', 1e-1,'MaxIter', 50000000,'Algorithm','active-set');
end

D.ratio = 1;
D.iter = 0;
D.ISCLOSE = ISCLOSE;
%constrained minimum of a function of several variables

tic
if ISCLOSE == 1
    [w_new, fval, exitflag] = fmincon(@prSamProtrusionOptSp, w, A, B ,[],[],[],[],[],options);
else
    %constrained minimum of a function of several variables
    Aeq = zeros(len_last, len_last);
    Beq = zeros(len_last,1);
    Aeq(1,1) = 1;
    Aeq(end,end) = 1;
    Beq(1) = 1;
    Beq(end) = len;
    [w_new, fval, exitflag] = fmincon(@prSamProtrusionOptSp, w, A, B ,Aeq,Beq,[],[],[],options);
end

pixel=fnval(D.sp_s,w_new)';
translate = pixel - pixel_tm1;

elasped_time = toc;

if ISCLOSE == 1
    pixel_tm1_upsamp=fnval(sp_s,[len_last_old + [1:len_last_old]])';
else
    pixel_tm1_upsamp=fnval(sp_s,[1:len_last_old])';
end
samp = [floor(w_new):floor(w_new)+len_old-1];
j = 1;
clear w_t_upsamp;
if ISCLOSE == 1
    for i=1:len_last_old
        if j<length(last_ds)-1
            if (i>=last_ds(j)) & (i<last_ds(j+1))
            else
                j=j+1;
            end
            w_t_upsamp(i) = ( w_new(j+1)*(i-last_ds(j)) + w_new(j)*(last_ds(j+1)-i) )/ (last_ds(j+1) - last_ds(j));
        else
            if (i>=last_ds(j)) & (i<last_ds(j+1))
                w_t_upsamp(i) = ( w_new(j+1)*(i-last_ds(j)) + w_new(j)*(last_ds(j+1)-i) )/ (last_ds(j+1) - last_ds(j));
            else
                if len_last_old - last_ds(j+1) == 0
                    w_t_upsamp(i) = w_new(end);
                else
                    w_t_upsamp(i) = ( samp(end)*(i-last_ds(j+1)) + w_new(j+1)*(len_last_old-i) )/ (len_last_old - last_ds(j+1));
                end
            end
        end
    end
else    
    for i=1:len_last_old
        if j<length(last_ds)-1
            if (i>=last_ds(j)) & (i<last_ds(j+1))
            else
                j=j+1;
            end
            w_t_upsamp(i) = ( w_new(j+1)*(i-last_ds(j)) + w_new(j)*(last_ds(j+1)-i) )/ (last_ds(j+1) - last_ds(j));
        else
            w_t_upsamp(i) = ( w_new(j+1)*(i-last_ds(j)) + w_new(j)*(last_ds(j+1)-i) )/ (last_ds(j+1) - last_ds(j));
        end
    end
    w_t_upsamp(end) = l2;
%     
%         else
%             if (i>=last_ds(j)) & (i<last_ds(j+1))
%                 w_t_upsamp(i) = ( w_new(j+1)*(i-last_ds(j)) + w_new(j)*(last_ds(j+1)-i) )/ (last_ds(j+1) - last_ds(j));
%             else
%                 if len_last_old - last_ds(j+1) == 0
%                     w_t_upsamp(i) = len_last_old;
%                 else
%                     w_t_upsamp(i) = ( samp(end)*(i-last_ds(j+1)) + w_new(j+1)*(len_last_old-i) )/ (len_last_old - last_ds(j+1));
%                 end
%             end
%         end
%     end
end
pixel_t_upsamp=fnval(D.sp_s,w_t_upsamp')';
translate_upsamp = pixel_t_upsamp - pixel_tm1_upsamp;

output.pixel_tm1_output = pixel_tm1_upsamp;
output.translate_output = translate_upsamp;
output.pixel_t_output = pixel_t_upsamp;

output.xy_normal = xy_normal;
output.time_elapsed = elasped_time;
