
function NoiseRegressingOut(Filename,Fileout,NRun,NumScantoIgnore,RemoveMean,TrendOrder,WMMaskName,CSFMaskName,VesselMaskName,MotionFile,Normalize_To_Std)

% Created by Babak Afshinpour
% Feb-2009

% Filename        :  The name of Input NIFTI file
% Fileout         :  The name of output NIFTI file
% NRun            :  The number of run which Filename contains. The number
%                        of scans in Filename must be dividable by the number of runs.
% NumScantoIgnore :  The Number of scans to be ignored.
% RemoveMean      :    0 or empty matrix [] : It doesn't remove the mean 
%                      1: It removes the mean from
%                      each voxels time-series, and all voxels time-series will have  zero mean.
% TrendOrder      :    0 : No detreding, if it is greater than zero, it determines the
%                      order of polynomial trends that you want to remove.
%                      A proper choice is 0<TrendOrder<5.
% WMMaskeName     : A binary NIFTI or ANALYZE binary image, it must be in
%                        the same dimension as Filename. if you don't want to
%                        remove mean of white-matter times-series, put an
%                        empty matrix [].
% CSFMaskName     : A binary NIFTI or ANALYZE binary image, it must be in
%                        the same dimension as Filename
%                        if you don't to remove mean of CSF times-series, put an
%                        empty matrix [].
% VesselMaskName  : A binary NIFTI or ANALYZE binary image, it must be in
%                        the same dimension as Filename
%                        if you don't to remove mean of Vessel times-series, put an
%                        empty matrix [].
% MotionFile      : A test file which contains six (maybe less or more) motion parameters. if
%                        you don't remove six motion parameters put an empty matrix [].
% Normalize_To_Std: If you want to normalize the times-series to their
%                   standard deviation, set this parameter 1, else set it 0 or empty
%                   matrix [].

% % AFNI_Path = '/software/Afni/';
% % addpath([cd '/NIFTI_20081205']);
% % ffcode = which('NoiseRegressingOut.m');
% % [pathstr] = fileparts(ffcode);
% % addpath([pathstr '/NIFTI_20081205']);
% % 
% % addpath /data/grady_lab/dnichol/Doug_GERM_Preprocess/NIfTI_20140122
% % 

if ~isempty(CSFMaskName)
    
    nii_csf = load_untouch_nii(CSFMaskName);
    CSF     = nii_csf.img;
else
    CSF = [];
end

if ~isempty(VesselMaskName)
    nii_v = load_untouch_nii(VesselMaskName);
    Vessel     = nii_v.img;
else
    Vessel = [];
end

if ~isempty(WMMaskName) && ischar(WMMaskName)
    nii_wm = load_untouch_nii(WMMaskName);
    WM     = nii_wm.img;
else
    if ~isempty(WMMaskName) && ~ischar(WMMaskName)
        WM = WMMaskName;
    else
        WM = [];
    end
end



[hdr, filetype, fileprefix, machine] = load_nii_hdr(Filename);
r = hdr.dime.dim(2);
c = hdr.dime.dim(3);
k = hdr.dime.dim(4);
t = hdr.dime.dim(5);


if fix(t/NRun)~=t/NRun
    error('The number of scans in Filename must be dividable by the number of runs (NRun)');
end
numscans  = t/NRun;

if ~isempty(MotionFile)
    MotionReg = load(MotionFile,'-ASCII');
else
    MotionReg = [];
end


if isempty(NRun)
    NRun = 1;
end
if NRun<=0
    NRun = 1;
end

if isempty(RemoveMean)
    RemoveMean = 0;
end

if isempty(Normalize_To_Std)
    Normalize_To_Std = 0;
end

if ~isempty(WM)
    [r1,c1,k1] = size(WM);
    if ((r1~=r) || (c1~=c) || (k1~=k))
        error('WhiteMatter Mask and Filename must have same dimension');
    end
end

if ~isempty(CSF)
    [r1,c1,k1] = size(CSF);
    if ((r1~=r) || (c1~=c) || (k1~=k))
        error('CSF Mask and Filename must have same dimension');
    end
end

if ~isempty(Vessel)
    [r1,c1,k1] = size(Vessel);
    if ((r1~=r) || (c1~=c) || (k1~=k))
        error('Vessel Mask and Filename must have same dimension.');
    end
end

if ~isempty(MotionReg)
    [t1,n1] = size(MotionReg);
    if (t1~=t)
        error('Error in Motion parameter file size.');
    end
end

[path,name,ext] = fileparts(Fileout);
if isempty(path)
Fileout = name;
else
    Fileout = [path '/' name];
end
for i = 1:NRun

    clear Ts nii D Ts1
    
    display(sprintf('please wait loading run %d',i));
    nii = load_untouch_nii(Filename,numscans*(i-1)+1:numscans*i);
    display(sprintf('run %d is loaded.',i));
    D   = double(nii.img);
    Ts  = reshape(D,r*c*k,t/NRun);
    
    X = [];
    Ts1 = Ts(:,1+NumScantoIgnore:end);
    if TrendOrder~=0
        X = polynomial_regressors(TrendOrder,numscans - NumScantoIgnore);
    end 
    if RemoveMean~=0
        X = [X ones(numscans-NumScantoIgnore,1)];
    end
    if ~isempty(MotionReg)
        temp = MotionReg(numscans*(i-1)+1+NumScantoIgnore:numscans*i,:);
        temp = temp - repmat(mean(temp),numscans - NumScantoIgnore,1);
        temp = temp ./ repmat(sqrt(sum(temp.^2)),numscans - NumScantoIgnore,1);
        X = [X temp];
    end
    if ~isempty(WM);
        temp = noise_reg(Ts1,WM);
        if TrendOrder~=0
            temp = regressout(temp',X(:,1:TrendOrder))';
        end
        X = [X temp];
    end
    if ~isempty(CSF)
        temp = noise_reg(Ts1,CSF);
        if TrendOrder~=0
            temp = regressout(temp',X(:,1:TrendOrder))';
        end
        X = [X temp];
    end
    if ~isempty(Vessel)
        temp = noise_reg(Ts1,Vessel);
        if TrendOrder~=0
            temp = regressout(temp',X(:,1:TrendOrder))';
        end
        X = [X temp];
    end
    [Ts1,EDOF(i)] = regressout(Ts1,X);

    %%% Normalize to standard deviation
    if Normalize_To_Std == 1
        temp_std = std(Ts1,0,2);
        temp_std(temp_std==0) = 1;
        Ts1 = Ts1 ./ repmat(temp_std,1,numscans - NumScantoIgnore);    
    end

    display(sprintf('Efficeint degree of freedom for run %d is %.2f',i,EDOF(i)));
    nii_out = nii;
    Ts(:,1+NumScantoIgnore:end) = Ts1;
    D = reshape(Ts,r,c,k,t/NRun);
    nii_out.img = double(D);
    if NRun>1
        save_untouch_nii(nii_out,[Fileout num2str(i) '.nii']);
    else
        save_untouch_nii(nii_out,[Fileout '.nii']);
    end
end

if NRun>1
    str = [];
    for i = 1:NRun
        str = [str Fileout num2str(i) '.nii '];
    end
    afni_path = fileparts(AFNI_Path);
    dos([afni_path '/3dTcat -prefix ' Fileout ' '  str]);
    for i = 1:NRun
        delete([Fileout num2str(i) '.nii']);
    end
    dos([afni_path '/3dAFNItoNIFTI ' Fileout '+orig']);
    delete([Fileout '+orig.BRIK']);
    delete([Fileout '+orig.HEAD']);
end

function X = polynomial_regressors(TrendOrder,numscan)
lin = [1:numscan]';
X   = [];
for i = 1:TrendOrder
    temp = lin.^i;
    temp = temp - mean(temp);temp = temp/norm(temp);
    X = [X temp];
end


function X = noise_reg(Ts1,Mask)

index = find(Mask~=0);
Mask_Ts = Ts1(index,:);
X = mean(Mask_Ts)';
X = X - mean(X);X = X./norm(X);


function [Tsout,EDOF] = regressout(Ts,X)

[t,n] = size(X);

% [U,S,V] = svd(X,0);
% dd = abs(diag(S));
% rcond = min(dd)/max(dd);
% if rcond<1e-12
%     [dds,sindex] = sort(dd);
%     for i = 1:length(dd)
%         if (dds(i)/dds(end))>1e-10
%             break;
%         end
%     end
% 
%     X = U(:,sindex(i:end))*S(sindex(i:end),sindex(i:end))*V(sindex(i:end),sindex(i:end))';
% end
       
SignalSpace = eye(t) - X*inv(X'*X)*X';
EDOF = trace(SignalSpace);
Tsout = Ts*SignalSpace';

