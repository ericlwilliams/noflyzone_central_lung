function NFz_AtlasProbability

tic; %close all;

screen_size=get(0,'ScreenSize');
ss_two2two = [screen_size(3)/2 0 screen_size(4) screen_size(4)];
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_print = true;
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
        fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2b',a2b{1},'_acute_data.mat'];

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
  hold on;
  
  
  set(gca,'FontSize',20)
  xlabel('Physical Dose [Gy]','FontSize',22);
  ylabel('Volume [cc]','FontSize',22);
  
%  ylim([0.1,15]);
    %set(gca,'YTick',[1:2:21]); set(gca,'YTickLabel',0:10:100);
%   title('Probability of true Esophageal Toxicity $\geq$ 20\%','interpreter','latex');
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_basename,'_allfx_prob'],'-pdf');
   disp(['Saving ',fig_basename,'_allfx_prob.pdf...']);
  end

  % 5 fraction
    CGobj = CGobj.fCrudeAtlas_DVH(5);
    CGobj = CGobj.fBetaCumulativeProbability_DVH();

  cur_fig=figure(2); clf reset;
  %set(gcf,'Position',ss_four2three);
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGobj.fProbabilityFig_DVH('');
  hold on;
  % bin locations determined by hand
  % eg [~,med_bin] = min(abs(CGobj.mBinsDose - 11.7))
  h_fit=plot(83,42,'xk','MarkerSize',20,'LineWidth',4);
  h_med=plot(52,42,'+k','MarkerSize',20,'LineWidth',4);
    hold off;
  lgnd=legend([h_med h_fit],...
      ['Median split* of D$_{5\rm{cc}}$ [Gy$_{10}$] for complete dataset',10,'in physical dose'],...
      'Split* at 20\% rate from regression using D$_{5\rm{cc}}$ [Gy$_{10}$]',...
      'Location','North');
  set(lgnd,'FontSize',16);
  set(lgnd,'interpreter','latex');
  %ylim([0.1,15]);
  %text(0,-0.11,['*Back calculated to physical',10', 'dose delivered in 5 fractions'],...
  text(0.45,0.73,['*Back calculated to physical dose',10,'delivered in 5 fractions'],...
      'Units','normalized','interpreter','latex','fontsize',14);
   set(gca,'FontSize',20)
  xlabel('Physical Dose [Gy]','FontSize',22);
  ylabel('Volume [cc]','FontSize',22);
%   title('Probability of true Esophageal Toxicity $\geq$ 20\%','interpreter','latex')
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_basename,'_5fx_prob'],'-pdf');
   disp(['Saving ',fig_basename,'_5fx_prob.pdf...']);
  end
  
  % 4 fraction
    CGobj = CGobj.fCrudeAtlas_DVH(4);
    CGobj = CGobj.fBetaCumulativeProbability_DVH();

  cur_fig=figure(3); clf reset;
  %set(gcf,'Position',ss_four2three);
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGobj.fProbabilityFig_DVH('');
  hold on;
  % bin locations determined by hand
  h_fit=plot(79,42,'xk','MarkerSize',20,'LineWidth',4);
  h_med=plot(50,42,'+k','MarkerSize',20,'LineWidth',4);
  hold off;
  lgnd=legend([h_med h_fit],...
      ['Median split* of D$_{5\rm{cc}}$ [Gy$_{10}$] for complete dataset',10,'in physical dose'],...
      'Split* at 20\% rate from regression using D$_{5\rm{cc}}$ [Gy$_{10}$]',...
      'Location','North');
  set(lgnd,'FontSize',16);
  set(lgnd,'interpreter','latex');
  %ylim([0.1,15]);
%   text(0,-0.11,['*Back calculated to physical',10', 'dose delivered in 4 fractions'],...
%       'Units','normalized','interpreter','latex');
 text(0.45,0.73,['*Back calculated to physical dose',10,'delivered in 4 fractions'],...
      'Units','normalized','interpreter','latex','fontsize',14);
     set(gca,'FontSize',20)
  xlabel('Physical Dose [Gy]','FontSize',22);
  ylabel('Volume [cc]','FontSize',22);
%   title('Probability of true Esophageal Toxicity $\geq$ 20\%','interpreter','latex');
  
%  ylim([0.1,15]);
  
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_basename,'_4fx_prob'],'-pdf');
   disp(['Saving ',fig_basename,'_4fx_prob.pdf...']);
  end
  
    % 3 fraction
    CGobj = CGobj.fCrudeAtlas_DVH(3);
    CGobj = CGobj.fBetaCumulativeProbability_DVH();

  cur_fig=figure(4); clf reset;
  %set(gcf,'Position',ss_four2three);
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGobj.fProbabilityFig_DVH('');
  hold on;
  % bin locations determined by hand
  h_fit=plot(74,42,'xk','MarkerSize',20,'LineWidth',4);
   h_med=plot(47,42,'+k','MarkerSize',20,'LineWidth',4);
  hold off;
    lgnd=legend([h_med h_fit],...
      ['Median split* of D$_{5\rm{cc}}$ [Gy$_{10}$] for complete dataset',10,'in physical dose'],...
      'Split* at 20\% rate from regression using D$_{5\rm{cc}}$ [Gy$_{10}$]',...
      'Location','North');
  set(lgnd,'FontSize',16);
  set(lgnd,'interpreter','latex');

  %ylim([0.1,15]);
%   text(0,-0.11,['*Back calculated to physical',10', 'dose delivered in 3 fractions'],...
%       'Units','normalized','interpreter','latex');
   text(0.45,0.73,['*Back calculated to physical dose',10,'delivered in 3 fractions'],...
      'Units','normalized','interpreter','latex','fontsize',14);
     set(gca,'FontSize',20)
  xlabel('Physical Dose [Gy]','FontSize',22);
  ylabel('Volume [cc]','FontSize',22);
%   title('Probability of true Esophageal Toxicity $\geq$ 20\%','interpreter','latex');
  
%  ylim([0.1,15]);
  
  
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_basename,'_3fx_prob'],'-pdf');
   disp(['Saving ',fig_basename,'_3fx_prob.pdf...']);
  end
  
 
    end
end
    