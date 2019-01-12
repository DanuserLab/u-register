function [output] = prSamProtrusion(inputParam)
% prSamProtrusion is funciton to perform new Protrusion model based
% on constrainted optimization. The input edge pixels are the cell boundary
% after image segmentation, and the mask at t-1 is used to calculating unit
% normal vectors along the edge.  The output consists of out egde pixels 
%
% Input:    inputParam.list_last        input edge pixels at t-1
%           inputParam.list             input edge pixels at t 
%           inputParam.maskTM1          mask at t-1 for calculating unit normal vectors along the edge
%           inputParam.batch_processing if true, disables progress bars etc.
%
% Output:   output.pixel_tm1_output     output edge pixels at t-1 (after spline parameterization)
%           output.translate_output     translation (protrusion) from t-1 to t
%           output.pixel_t_output       output edge pixels at t (after spline parameterization)
%           output.xy_normal            unit normal vectors along the edge
%           output.time_elapsed         computational time for protrusion calculation 
%
% Last updated: August 26, 2008 by Shann-Ching Chen, LCCB
% See also: protrusionAnalysis, prSamProtrusionDownSpl, prSamProtrusionOptSp
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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

tic
pixel_list_last = inputParam.pixel_list_last;
pixel_list = inputParam.pixel_list;
start_position = pixel_list_last(1,:);
ISCLOSE = inputParam.ISCLOSE;

if ~isfield(inputParam,'batch_processing') %Added batch mode to disable output - HLE
    inputParam.batch_processing = false;
end

l1 = size(pixel_list_last,1);
l2 = size(pixel_list,1);

if ISCLOSE == 1
    % For close contour, rotate the pixels so that the head of the splines is at the intersection of two splines

    [x0_rot,y0_rot,iout_rot,jout_rot] = intersectionsSAM(pixel_list_last(:,1),pixel_list_last(:,2),pixel_list(:,1),pixel_list(:,2),1);
    
    if isempty(x0_rot) == 0
        idx1 = iout_rot(1)+0.0001;    idx2 = jout_rot(1)+0.0001;
    else
        % if no intersection, use downsample
        dl_rate = ceil(size(inputParam.pixel_list_last,1)/100);       
        output = prSamProtrusionDownSpl(inputParam, dl_rate);
        return;
    end
    
    pixel_list_last = [pixel_list_last(ceil(idx1):end,:); pixel_list_last(1:floor(idx1),:)];
    pixel_list = [pixel_list(ceil(idx2):end,:); pixel_list(1:floor(idx2),:)];
    
    pixel_list_last_extended = [pixel_list_last; pixel_list_last; pixel_list_last];
    sp_tm1 = spaps([1:3*l1],pixel_list_last_extended',TOLERANCE);
    spline_pixel_list_last=fnval(sp_tm1,[l1+1:2*l1])';
    
    pixel_list_extended = [pixel_list; pixel_list; pixel_list];
    sp_t = spaps([1:3*l2],pixel_list_extended',TOLERANCE);
    spline_pixel_list=fnval(sp_t,[l2+1:2*l2])';
else
    sp_tm1 = spaps([1:l1],pixel_list_last',TOLERANCE);
    spline_pixel_list_last=fnval(sp_tm1,[1:l1])';    

    sp_t = spaps([1:l2],pixel_list',TOLERANCE);
    spline_pixel_list=fnval(sp_t,[1:l2])';    
end

% What happened if no intersection?
x1 = spline_pixel_list_last(:,1);   y1 = spline_pixel_list_last(:,2);
x2 = spline_pixel_list(:,1);        y2 = spline_pixel_list(:,2);
[x0,y0,iout_i,jout_i] = intersectionsSAM(x1,y1,x2,y2,1);

remove_idx = union(find(diff(iout_i) < 0.01),find(diff(jout_i) < 0.01));
x0(remove_idx) = [];
y0(remove_idx) = [];
iout_i(remove_idx) = [];
jout_i(remove_idx) = [];


% figure(1);  plot(x1, y1, x2, y2);
% axis equal;
% hold on;
% for i=1:length(iout_i)
%     %        plot(x1(round(iout(i))), y1(round(iout(i))),'bo');
%     plot(x0(i), y0(i),'ro');
% %    pause;
% end
% plot(x0, y0,'bo');
% axis equal;
% hold off;
% 
% keyboard;
if isempty(x0)
    % proceed with Downsampling;   
    dl_rate = ceil(size(inputParam.pixel_list_last,1)/100);        
    output = prSamProtrusionDownSpl(inputParam, dl_rate);   
    return;
else
    iout = [];   jout = [];
    if iout_i(1) == 1
        iout = [iout_i];
        jout = [jout_i];
    else
        iout = [1; iout_i];
        jout = [1; jout_i];
    end
    if iout_i(end) < l1
        iout = [iout; l1];
        jout = [jout; l2];
    end

    segLenTM1 = 0; segLenT = 0;
    segNum = 1;
    segIdxTM1 = cell(100,1); segIdxT = cell(100,1);
    i = 1;
    
    while 1
        startIdx = i;
        segLenTM1(segNum) = iout(i+1) - iout(i);
        segLenT(segNum) = jout(i+1) - jout(i);
        segIdxTM1{segNum} = unique([iout(i), ceil(iout(i)):floor(iout(i+1)), iout(i+1)]);
        segIdxT{segNum}  = unique([jout(i),  ceil(jout(i)):floor(jout(i+1)), jout(i+1)]);
        i = i + 1;

        if i==length(iout)
            segIdxTM1(segNum+1:end) = [];
            segIdxT(segNum+1:end) = [];
            break;
        else
            segNum = segNum + 1;
        end
    end

    idx = union(find(segLenTM1 <= 3),find(segLenT <= 3));
    for i=length(idx):-1:1
        if idx(i) == segNum
            % combine the last segment to the previous one
            segLenTM1(segNum-1) = segLenTM1(segNum-1) + segLenTM1(segNum);
            segIdxTM1{segNum-1} = [segIdxTM1{segNum-1}(1:end) segIdxTM1{segNum}(2:end)];
            segIdxT{segNum-1} = [segIdxT{segNum-1}(1:end) segIdxT{segNum}(2:end)];
            segIdxT(segNum) = [];
            segIdxTM1(segNum) = [];
            segNum = segNum - 1;
            segLenTM1(segNum) = [];
        else
            % combine idx(i)th segment to idx(i)-1 th segment
            segLenTM1(idx(i)+1) = segLenTM1(idx(i)+1) + segLenTM1(idx(i));
            segIdxTM1{idx(i)+1} = [segIdxTM1{idx(i)}(1:end) segIdxTM1{idx(i)+1}(2:end)];
            segIdxT{idx(i)+1} = [segIdxT{idx(i)}(1:end) segIdxT{idx(i)+1}(2:end)];
            segIdxT(idx(i)) = [];
            segIdxTM1(idx(i)) = [];
            segNum = segNum - 1;
            segLenTM1(idx(i)) = [];
        end
    end

    segLenTM1 = zeros(1,segNum);
    segLenT = zeros(1,segNum);    
    for j = 1:segNum  
        segLenTM1(j) = segIdxTM1{j}(end) - segIdxTM1{j}(1);
        segLenT(j) = segIdxT{j}(end) - segIdxT{j}(1);
        segIdxTM1{j} = [segIdxTM1{j}(1) unique(round(segIdxTM1{j}(2:end-1))) segIdxTM1{j}(end)];
        segIdxT{j} = [segIdxT{j}(1) unique(round(segIdxT{j}(2:end-1))) segIdxT{j}(end)];        
        if segIdxTM1{j}(2) - segIdxTM1{j}(1) <= 0.5
            segIdxTM1{j}(2) = [];
        end
        if segIdxTM1{j}(end) - segIdxTM1{j}(end-1) < 0.5
            segIdxTM1{j}(end-1) = [];
        end      
        if segIdxT{j}(2) - segIdxT{j}(1) <= 0.5
            segIdxT{j}(2) = [];
        end
        if segIdxT{j}(end) - segIdxT{j}(end-1) < 0.5
            segIdxT{j}(end-1) = [];
        end             
    end        
end
pixel_tm1_output = [];
translate_output = [];
pixel_t_output = [];
d_s_all = [];

if ~inputParam.batch_processing
    fprintf(1,'Processing Frame %d-%d (%d fragments)......',inputParam.time_idx-1,inputParam.time_idx,segNum);
    progressbar;
end

for j = 1:segNum
    segNumTM1 = length(segIdxTM1{j});
    if ~inputParam.batch_processing
        progressbar(j/segNum);
    end
    clear global D;
    global D;
    D.sp_s = sp_t;
    D.iter = 0;
    D.ISCLOSE = 0;
    if ISCLOSE == 1
        d_s = segIdxTM1{j} + l1;
    else
        d_s = segIdxTM1{j};
    end
    
    Samp = inputParam.dl_rate;      Intvl = Samp - 1;    
    
    if segNumTM1 < Samp
        ds_down = d_s;
    else
        spacing = segLenTM1(j)/Intvl;
        ds_down = [d_s(1):spacing:d_s(end)];
        if length(ds_down) < Samp
            ds_down = [ds_down d_s(end)];
        end
    end

    pixel_tm1_seg=fnval(sp_tm1,d_s)';
    pixel_tm1_seg_down=fnval(sp_tm1,ds_down)';
    D.pixel_tm1 = pixel_tm1_seg_down;
    D.ratio = 1;
    if ISCLOSE == 1    
%        w = (segIdxT{j}(1):(segIdxT{j}(end) - segIdxT{j}(1))/(length(d_s)-1):segIdxT{j}(end))' + l2;
        w = (segIdxT{j}(1):(segIdxT{j}(end) - segIdxT{j}(1))/(length(ds_down)-1):segIdxT{j}(end))' + l2;
    else
%        w = (segIdxT{j}(1):(segIdxT{j}(end) - segIdxT{j}(1))/(length(d_s)-1):segIdxT{j}(end))';
        w = (segIdxT{j}(1):(segIdxT{j}(end) - segIdxT{j}(1))/(length(ds_down)-1):segIdxT{j}(end))';
    end
    
%    varNum = length(segIdxTM1{j});
    varNum = length(ds_down);
    A = zeros(varNum, varNum);
    for i=1+1:varNum
        A(i,i-1:i) = [1 -1];
    end
    A(1,1) = -1;
    B = zeros(varNum,1);
    if isempty(strfind(version,'R2008'))
        options = optimset('largescale','off','Display','off','TolCon', 1e-1,'MaxIter', 50000000);
    else
        options = optimset('largescale','off','Display','off','TolCon', 1e-1,'MaxIter', 50000000,'Algorithm','active-set');
    end

    %constrained minimum of a function of several variables
    Aeq = zeros(varNum, varNum);
    Beq = zeros(varNum,1);
    Aeq(1,1) = 1;
    Aeq(end,end) = 1;
    Beq(1) = w(1);
    Beq(end) = w(end);   
    
    [w_new, fval, exitflag] = fmincon(@prSamProtrusionOptSp, w, A, B ,Aeq,Beq,[],[],[],options);
    
    if segNumTM1 < Samp
        w_t_upsamp = w_new;
    else
        clear w_t_upsamp;
        w_t_upsamp = zeros(segNumTM1,1);
        w_t_upsamp(1) = w_new(1);
        w_t_upsamp(segNumTM1) = w_new(end);

        k=1;
        for i=2:segNumTM1-1
            if (d_s(i)>=ds_down(k)) & (d_s(i)<ds_down(k+1))

            else
                k=k+1;
            end
            w_t_upsamp(i) = ( w_new(k)*(ds_down(k+1) - d_s(i)) + w_new(k+1)*(d_s(i) - ds_down(k)) )/ (ds_down(k+1) - ds_down(k));
        end
    end
    
    pixel_t_seg=fnval(sp_t,w_t_upsamp)';   
    translate_seg = pixel_t_seg - pixel_tm1_seg;

    if j == segNum
        pixel_tm1_output = [pixel_tm1_output; pixel_tm1_seg];
        translate_output = [translate_output; translate_seg];
        pixel_t_output = [pixel_t_output; pixel_t_seg];
        d_s_all = [d_s_all d_s];
    else
        pixel_tm1_output = [pixel_tm1_output; pixel_tm1_seg(1:end-1,:)];
        translate_output = [translate_output; translate_seg(1:end-1,:)];
        pixel_t_output = [pixel_t_output; pixel_t_seg(1:end-1,:)];
        d_s_all = [d_s_all d_s(1:end-1)];  
    end
end

%    h = figure; plot(pixel_tm1_output(:,1), pixel_tm1_output(:,2),'ro-', ...
%         pixel_t_output(:,1), pixel_t_output(:,2),'gx-');
%             hold on;
%             quiver(pixel_tm1_output(:,1), pixel_tm1_output(:,2), translate_output(:,1), translate_output(:,2), -1);
%             hold off;
%             axis equal
            
if ISCLOSE == 1
    % Rotate the pixels back to the starting position
    [tmp, tmpI] = min(sum((repmat(start_position, [size(pixel_tm1_output,1) 1]) - pixel_tm1_output).^2,2));
    if tmpI > 1
        pixel_tm1_output = [pixel_tm1_output(tmpI:end,:); pixel_tm1_output(1:tmpI-1,:)];
        translate_output = [translate_output(tmpI:end,:); translate_output(1:tmpI-1,:)];
        pixel_t_output = [pixel_t_output(tmpI:end,:); pixel_t_output(1:tmpI-1,:)];
        d_s_all = [d_s_all(tmpI:end) d_s_all(1:tmpI-1)];
    end
end
xy_normal = prSamCalculateNormal(inputParam.maskTM1, sp_tm1, d_s_all);

output.pixel_tm1_output = pixel_tm1_output;
output.translate_output = translate_output;
output.pixel_t_output = pixel_t_output;
output.xy_normal = xy_normal;
output.time_elapsed = toc;
output.spline_pixel_list = spline_pixel_list;
output.spline_pixel_list_last = spline_pixel_list_last;
if ~inputParam.batch_processing
    fprintf(1, '\tDone in %.3f seconds\n', output.time_elapsed);
end
