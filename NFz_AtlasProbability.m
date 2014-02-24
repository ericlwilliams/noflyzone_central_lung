function NFz_AtlasProbability

tic; %close all;

screen_size=get(0,'ScreenSize');
ss_two2two = [screen_size(3)/2 0 screen_size(4) screen_size(4)];
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_print = false;
fig_loc = 'Z:/elw/MATLAB/nfz_analy/slides/figures/latest/';

%structures = {'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'PBT' 'LUNGS'};
structures = {'ESOPHAGUS'};
%toxicities = {'rp','pultox','esotox'};
toxicities = {'esotox'};

fp = 'Z:\elw\MATLAB\nfz_analy\meta_data\';

%a2b = {'Inf' '3'};
%a2b = {'10'};
a2b = {'Inf'};


for i=1:length(toxicities)
    
    for j=1:length(structures)
        cur_fig_ctr = (10*i)+j-1;
        
        fprintf('\n');
        disp(['Tox: ',toxicities{i}]);
        disp(['Struct: ',structures{j}]);    
        disp(['Counter: ',num2str(cur_fig_ctr)]);
        fprintf('\n');
      
        fig_basename = [fig_loc,'nfz_',...
                        structures{j},'_',...
                        toxicities{i},'_a2b',...
                        a2b{1}];
    
     %% load data
        fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2b',a2b{1},'_data.mat'];
        load(strcat(fp,fn),'CGobj_org');
        CGobj = CGobj_org;
        LymanN = log10(CGobj.mLymanN);
        CGobj.mLymanN = LymanN;
        f = [CGobj.mGrp.mFlgCensor];
        
        grp = CGobj.mGrp;
        comps = ~f;
        eud = [grp(:).mEUD];

    CGobj = CGobj.fCrudeAtlas_DVH(-1);
    CGobj = CGobj.fBetaCumulativeProbability_DVH();

cur_fig=figure(1); clf reset;
  %set(gcf,'Position',ss_four2three);
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGobj.fProbabilityFig_DVH('');
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'msk_rp20pct'],'-pdf');
   disp(['Saving ',fig_loc,'msk_rp20pct.pdf...']);
  end

    end
end
    