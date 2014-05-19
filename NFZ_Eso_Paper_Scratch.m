function NFZ_Eso_Paper_Scratch
%% Scratch code for esophageal paper
% created 2014_05_10

tic;close all;

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
cur_fig_ctr = 1;
do_print = true;
do_acute = true;


fig_loc = 'Z:/elw/MATLAB/nfz_analy/slides/figures/latest/';

structures = 'ESOPHAGUS';
str_colors = 'k';

toxicities = 'esotox';

comp_rate = 0.2;

fp = 'Z:\elw\MATLAB\nfz_analy\meta_data\';

a2b = '10';
cur_a2b = 10;

fprintf('\n');
disp(['Tox: ',toxicities]);
disp(['Struct: ',structures]);

fprintf('\n');

fig_basename = [fig_loc,'nfz_',...
    structures,'_',...
    toxicities,'_a2b',...
    a2b];

if do_acute,
    
    fn = ['NFZ_',structures,'_',toxicities,'_a2b',a2b,'_acute_data.mat'];
    disp(['Loading ACUTE data ',fn,'...']);
else
    fn = ['NFZ_',structures,'_',toxicities,'_a2b',a2b,'_data.mat'];
    disp(['Loading: ',fn,'...']);
end


load(strcat(fp,fn),'CGobj_org');
CGobj = CGobj_org;
clear CGobj_org;

%
LymanN = log10(CGobj.mLymanN);
CGobj.mLymanN = LymanN;

pttotal = ones(CGobj.mNumInGrp,1);
ptcomp = ones(CGobj.mNumInGrp,1); ptcomp([CGobj.mGrp.mFlgCensor])=0;

disp(['Rate of esophagitis: ',...
    '(',num2str(sum(ptcomp)),'/',num2str(sum(pttotal)),') = ',...
    num2str(100*(sum(ptcomp)./sum(pttotal)),7),'%']);






d_bins = [CGobj.mGrp.mDoseBins_LQ];
v_bins = [CGobj.mGrp.mVolCum];

d05s = zeros(size(v_bins,2),1);
d35s = zeros(size(v_bins,2),1);
dmax = zeros(size(v_bins,2),1);

end

%         for k=1:size(v_bins,2) % for each patient
%             vol = v_bins(:,k);
%
%             %vol = vol./max(vol);
%             %v05 = vol<0.05;
%
%
%             %v05 = vol<5; %% < 5 cc
%             %v05 = vol<3.5; %% < 5 cc
%
%             v05 = vol<5; %% < 5 cc
%             d05_inds = find(v05);
%             dose = d_bins(:,k);
%             d05 = min(dose(d05_inds));
%             d05s(k) = d05;
%
%             v35 = vol<3.5; %% < 5 cc
%             d35_inds = find(v35);
%             %dose = d_bins(:,j);
%             d35 = min(dose(d35_inds));
%             d35s(k) = d35;
%
%             vol(~vol)=nan;
%             nan_inds = find(isnan(vol));
%             if ~isempty(nan_inds)
%                 min_ind = nan_inds(1)-1;
%             else
%                 min_ind = length(vol);
%             end
%             %[~,min_ind] = min(vol);
%             dmax(k) = dose(min_ind);
%         end
%
%         allFx = [CGobj.mGrp.mFxNum];
%         nFx = unique(allFx);
%
%         for i=1:length(nFx)
%             cur_fx = nFx(i);
%             fx_inds = logical(allFx==cur_fx);
%             cur_d35s = d35s(fx_inds);
%             [~,~,s]=glmfit(cur_d35s,[ptcomp(fx_inds) pttotal(fx_inds)],'binomial','link','logit');
%             disp(['Fraction: ',num2str(cur_fx),' LR p-val: ',num2str(s.p(2))]);
%         end
%     end
% end

