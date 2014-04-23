function NFz_Analysis
tic;
save_result = true;
mm_do_print = false;
%beta2alpha=[1/3];a2b_corr = 'BED';
beta2alpha=[1/10];a2b_corr = 'BED';
%beta2alpha=[0];a2b_corr = 'PHYS';

tox_grade = 2;% toxicity defined as >= tox_grade

do_tcp_exclude = false;
do_lbed_exclude = false;
do_late_exclude = true;

dosestep=1;
volstep=10; % volume step in cc
timestep=3; % time step in month

if do_tcp_exclude
    filepathname='Z:/elw/MATLAB/original_data/NFZ/Central_tumor_dataset_tcp_4-22-13';
else
    filepathname='Z:/elw/MATLAB/original_data/NFZ/Central_tumor_dataset_4-22-13';
end

if do_lbed_exclude
    filepathname='Z:/elw/MATLAB/original_data/NFZ/Central_tumor_dataset_lbed_4-22-13';
else
    filepathname='Z:/elw/MATLAB/original_data/NFZ/Central_tumor_dataset_4-22-13';
end

 %full   
%structures = {'CLUNG' 'ILUNG' 'ESOPHAGUS' 'HEART'...
%                'NFZ' 'PBT' 'TRACHEA' 'LUNGS'};
                
%structures = {'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'PBT' 'LUNGS'};

%structures = {'PTV' 'GTV'};

structures = {'ESOPHAGUS'};
%structures = {'ILUNG' 'HEART' 'NFZ' 'PBT' 'LUNGS' 'ESOPHAGUS'};
%structures = {'HEART'};

%toxicities = {'lclfail'};
%tox_columns = {'Local Failure'};
%toxdate_columns = {'Local Failure Date'};


%toxicities = {'pultox'};
%tox_columns = {'Highest pulmonary toxicity or RT pneumonitis grade'};
%toxdate_columns = {'Pulmonary tox date'};
 
%toxicities = {'rp'};
%tox_columns = {'RT pneumonitis grade'};
%toxdate_columns = {'RP date'};

toxicities = {'esotox'};
tox_columns = {'Highest esophagitis grade'};
toxdate_columns = {'Highest esophagitis date'};

do_exclude_fu=false;


try
    load(filepathname,'xlsraw');
catch
    error(['Error loading NFz DATASET']);
end

% pick up related data
PtInfo = classDataFromXls();
PtInfo.mXlsRaw = xlsraw;

%% MRN
PtInfo.mLabel = 'MRN';
PtInfo = PtInfo.fExtractColData();
flgptcode = PtInfo.mFlg;
ptcode = PtInfo.mData;

%% gender
PtInfo.mLabel = 'Sex';
PtInfo = PtInfo.fExtractColData();
flggender = PtInfo.mFlg;
ptgender = PtInfo.mData;
ptgender = strcmp(ptgender,'Male');

%% Number of Fractions
PtInfo.mLabel = 'Number of Fractions';
PtInfo = PtInfo.fExtractColData();
flgfx = PtInfo.mFlg;
fx = zeros(size(flgfx));
fx(flgfx) = cell2mat(PtInfo.mData(flgfx));

%% Delivered dose
PtInfo.mLabel = 'Total Dose (cGy)';
PtInfo = PtInfo.fExtractColData();
flgtx = PtInfo.mFlg;
tx = zeros(size(flgtx));
tx(flgtx) = cell2mat(PtInfo.mData(flgtx));

%% NFz (0- not in NFz, 1 - NFz)
PtInfo.mLabel = 'Central-Mediastinum ONLY';
PtInfo = PtInfo.fExtractColData();
flgnfz = PtInfo.mFlg;
nfz = zeros(size(flgnfz));
nfz(flgnfz) = ~cell2mat(PtInfo.mData(flgnfz)); % 1 - in NFz, 0 - not in NFz

%% Date of Birth
PtInfo.mLabel = 'Date of Birth';
PtInfo = PtInfo.fExtractColData();
flgdatebirth = PtInfo.mFlg;
datebirth = zeros(size(flgdatebirth));
f = PtInfo.mData(flgdatebirth);
datebirth(flgdatebirth) = datenum(f);

%% IGRT End Date
PtInfo.mLabel = 'IGRT End Date';
PtInfo = PtInfo.fExtractColData();
flgdateIGRT = PtInfo.mFlg;
dateIGRT = zeros(size(flgdateIGRT));
f = PtInfo.mData(flgdateIGRT);
dateIGRT(flgdateIGRT) = datenum(f);

%% IGRT Start Date
PtInfo.mLabel = 'IGRT Start Date';
PtInfo = PtInfo.fExtractColData();
flgdateStartIGRT = PtInfo.mFlg;
dateStartIGRT = zeros(size(flgdateStartIGRT));
f = PtInfo.mData(flgdateStartIGRT);
dateStartIGRT(flgdateStartIGRT) = datenum(f);

%% age (at start of treatment)
flgAge = (flgdateStartIGRT & flgdatebirth);
age = (dateStartIGRT-datebirth)./365;

%% Date of last follow-up
PtInfo.mLabel = 'Date of last follow-up';
PtInfo = PtInfo.fExtractColData();
flgdatefu = PtInfo.mFlg;
datefu = zeros(size(flgdatefu));
datefu(flgdatefu) = datenum(PtInfo.mData(flgdatefu));

%% Survival/Death info
PtInfo.mLabel = 'Survival (Dead=1, Alive=0)';
PtInfo = PtInfo.fExtractColData();
f1 = PtInfo.mFlg;
flgdeath = false(size(f1));
flgdeath(f1) = cell2mat(PtInfo.mData(f1));

PtInfo.mLabel = 'Date of death';
PtInfo = PtInfo.fExtractColData();
f = PtInfo.mFlg;
datedeath = inf(size(f));
datedeath(f&flgdeath) = datenum(PtInfo.mData(f&flgdeath));
flgdeath = f1&f;


PtInfo.mLabel = 'SurvivalDate';
PtInfo = PtInfo.fExtractColData();
f3 = PtInfo.mFlg;
datedeath(f3&~flgdeath) = datenum(PtInfo.mData(f3&~flgdeath));
%flgdeath = f3&flgdeath&~flgdeath;

%% Pulmonary toxicity

for i=1:length(toxicities)
    
    exclude_fu=false;
    if isequal(toxicities{i},'rp') && do_exclude_fu
        exclude_fu=true;
    end
    
    PtInfo.mLabel = tox_columns{i};
    PtInfo = PtInfo.fExtractColData();
    toxgrade = inf(size(PtInfo.mData));
    flgcensor = true(size(toxgrade));
    
    flgtox = PtInfo.mFlg;
    toxgrade(flgtox) = cell2mat(PtInfo.mData(flgtox));
    if isequal(toxicities{i},'lclfail')
        flgcensor(flgtox) = toxgrade(flgtox)~=1;
    else
        
        
    flgcensor(flgtox) = toxgrade(flgtox)<tox_grade;

    
    PtInfo.mLabel = toxdate_columns{i};
    PtInfo = PtInfo.fExtractColData();
    datepain = inf(size(PtInfo.mData));
    f = PtInfo.mFlg;
    datepain(f&~flgcensor) = datenum(PtInfo.mData(f&~flgcensor));
    
    
    %%ptflg
    
    flg = flgptcode & flggender & flgfx & flgtx & flgnfz & ...
        flgdatebirth & flgdatefu & flgdateIGRT & flgdateStartIGRT & ...
        flgAge & flgtox;
    
    ptcode=ptcode(flg);
    containsNumbers = cellfun(@isnumeric,ptcode);
    ptcode(containsNumbers) = cellfun(@num2str,ptcode(containsNumbers),...
        'UniformOutput',false);
    fx=fx(flg);
    tx=tx(flg);
    nfz=nfz(flg);
    datebirth=datebirth(flg);
    dateIGRT=dateIGRT(flg);
    datepain=datepain(flg);
    dateStartIGRT=dateStartIGRT(flg);
    datefu=datefu(flg);
    datedeath=datedeath(flg);
    toxgrade=toxgrade(flg);
    flgcensor=flgcensor(flg);
    ptgender = ptgender(flg);
    age = age(flg);
  
    
    %% Load DVHs
    for m=1:length(structures) % iterate with each definition
        
        if do_tcp_exclude
            fn=['Z:\elw\MATLAB\nfz_analy\meta_data\NFZ_TCP_',structures{m},'_DVHs.mat'];
        else
            fn=['Z:\elw\MATLAB\nfz_analy\meta_data\NFZ_',structures{m},'_DVHs.mat'];
            %%tmp
        end
        
        if do_lbed_exclude
            fn=['Z:\elw\MATLAB\nfz_analy\meta_data\NFZ_LBED_',structures{m},'_DVHs.mat'];
        else
            fn=['Z:\elw\MATLAB\nfz_analy\meta_data\NFZ_',structures{m},'_DVHs.mat'];
            %%tmp
        end
        
        
        load(fn,'DVH');
        disp(['Loading ',fn]);
        %% match patients DVH and .xls
        
        %%TMP
        %dvhcodes = DVH(:,4);
        %dvhcodes(cellfun(@(x) isempty(x), dvhcodes))='';
        %[f1,g1]=ismember(ptcode,dvhcodes);
        %[f2,g2]=ismember(dvhcodes,ptcode);
        
        [f1,g1]=ismember(ptcode,DVH(:,4));
        [f2,g2]=ismember(DVH(:,4),ptcode);
        if length(unique(g2(g2~=0)))~=length(find(g2))% || length(unique(g2~=0))~=length(find(g2))
            error('The patient ids are not unique');
        end
        if ~( all(f1) && all(f2) )
            disp('The patients in the spread sheet do not match with those in DVH curves, take the common part');
            % use common part only
            DVH=DVH(f2,:);
            ptcode=ptcode(f1);
            ptgender=ptgender(f1);
            fx=fx(f1);
            tx=tx(f1);
            nfz=nfz(f1);
            datebirth=datebirth(f1);
            dateStartIGRT=dateStartIGRT(f1);dateIGRT=dateIGRT(f1);
            datepain=datepain(f1);
            datefu=datefu(f1);
            flgcensor=flgcensor(f1);
            toxgrade=toxgrade(f1);
            age=age(f1);
            
            
            %tmp
            %dvhnames = DVH(:,4);
            %dvhnames(cellfun(@(x) isempty(x), dvhnames))='';
            [~,g1]=ismember(ptcode,dvhcodes);
            
            %[~,g1]=ismember(ptcode,DVH(:,1));
        end
        CIobjs = classOutcomeIndividual();
        
        
        CIobjs = repmat(CIobjs,[size(DVH,1),1]);
        cur_i=1;
        for n = 1:size(DVH,1)
            
            if do_exclude_fu && exclude_fu && flgcensor(n)
                if (datefu(n)-dateIGRT(n))<6*30
                    continue;
                end
            end
               
                
            CIobjs(cur_i,1).mID=ptcode{cur_i};
            CIobjs(cur_i,1).mGender = ptgender(cur_i);
            CIobjs(cur_i,1).mAgeAtTx = age(cur_i);
            CIobjs(cur_i,1).mFxNum=fx(cur_i);
            CIobjs(cur_i,1).mDoseTx=tx(cur_i);
            CIobjs(cur_i,1).mFlgNFz=nfz(cur_i);
            CIobjs(cur_i,1).mDosePerFx=tx(cur_i)/fx(cur_i);
            CIobjs(cur_i,1).mDateBirth = datebirth(cur_i);
            %CIobjs(cur_i,1).mDateDeath = datedeath(cur_i);
            CIobjs(cur_i,1).mDateBaseline = dateIGRT(cur_i);
            CIobjs(cur_i,1).mDateStartTx = dateStartIGRT(cur_i);
            CIobjs(cur_i,1).mDateComp = datepain(cur_i);
            CIobjs(cur_i,1).mDateLastFollowup = datefu(cur_i);
            %CIobjs(cur_i,1).mDateRelapse = datefailure(cur_i);
            
            if do_late_exclude, % mark late patients as censored
                if ~flgcensor(cur_i) && (datepain(cur_i)-dateIGRT(cur_i))>90,
                    CIobjs(cur_i,1).mFlgCensor = ~flgcensor(cur_i);
                    disp(['late: ',num2str(datepain(cur_i)-dateIGRT(cur_i)),10,...
                        'toxgrade: ',num2str(toxgrade(cur_i))])
                else
                    CIobjs(cur_i,1).mFlgCensor = flgcensor(cur_i);
                end
            else
            CIobjs(cur_i,1).mFlgCensor = flgcensor(cur_i);
            end
            
            
            %CIobjs(cur_i,1).mDoseBins_LQ=DVH{g1(cur_i),2}(:,1);
            
            
            vols = DVH{g1(cur_i),2}(:,3);
            diff_vols = DVH{g1(cur_i),2}(:,2);
            doses = DVH{g1(cur_i),2}(:,1);
           
            CIobjs(cur_i,1).mMaxVol=max(vols);
            
            if isequal(structures{m},'PTV') ||...
                    isequal(structures{m},'GTV'),
                
                %vols = vols./vols(1);
                %doses = doses./doses(end);
                
                vols = vols./max(vols);
                diff_vols = diff_vols./max(diff_vols);
                %doses = doses./max(doses);
            end
            
            CIobjs(cur_i,1).mDoseBins_org=doses;
            %CIobjs(cur_i,1).mVolDiff=DVH{g1(cur_i),2}(:,2);
            CIobjs(cur_i,1).mVolDiff=diff_vols;
            
            CIobjs(cur_i,1).mVolCum= vols;
            
            CIobjs(cur_i,1).mBeta2AlphaCorrection = a2b_corr;

            %CIobjs(cur_i,1).mBeta2AlphaCorrection = a2b_corr;    

          cur_i=cur_i+1;
            
        end
        
        CIobjs(cur_i:end)=[];
        
        if any([CIobjs.mFxNum]==0)
            error('Fraction number is zeros');
        end
        CGobj_org = classOutcomeAnalysis();
        CGobj_org.mLymanN = 10.^(-1:0.1:1)';
        
        CGobj_org.mStepDose = dosestep;
        CGobj_org.mStepVol = volstep;
        CGobj_org.mStepTime = timestep;
        CGobj_org = CGobj_org.fAddPatient(CIobjs);
        
        CGobj_org.mBeta2Alpha = beta2alpha(1);
        CGobj_org = CGobj_org.fCalculateEUD();

        %     %% Analyses
        
    
      
        

      
          

        %% TMP
%         %% Cox model
         disp(['CoxModel_DVH...']);
         CGobj_org = CGobj_org.fCoxModel_DVH();
%         
         disp(['Logistic Regression Analysis ...']); 
         CGobj_org = CGobj_org.fLogisticRegressionExact_EUD();  
  
         %% Full gEUD grid scan
%          LogisticBetaRange = {(-10:0.1:10)'; (-1:0.01:1)'};
%          CGobj_org.mLogisticRegressionGridBetaRange = LogisticBetaRange;
%         
         
        
        disp(['LogRankVDx_DVH...']); % V_{x} with empty patients set to zero
        CGobj_org = CGobj_org.fLogRankVDx_DVH();
         

        disp(['OveralCompCurve...']);
        CGobj_org = CGobj_org.fOverallCompCurve();

        %LogisticBetaRange = {(-10:0.1:10)'; (-1:0.01:1)'};
        LogisticBetaRange = {(-10:0.05:10)'; (-1:0.005:1)'};
        CGobj_org.mLogisticRegressionGridBetaRange = LogisticBetaRange;
        disp(['Logistic Grid']);
        CGobj_org = CGobj_org.fLogisticRegressionGridExact_EUD();
        
        
%         disp(['gEUD Mixture Model...']);
%         fig_loc='Z:\elw\MATLAB\nfz_analy\slides\figures\latest\';
%         CGobj_org = CGobj_org.fLogisticRegressionGridExactMixtureModel_2_EUD(mm_do_print,fig_loc);
%         
%             
      
        
        % save result
        CGstrct = ObjToStruct(CGobj_org);
        
        if do_tcp_exclude
            fn=['Z:\elw\MATLAB\nfz_analy\meta_data\NFZ_',...
                structures{m},'_',...
                toxicities{i},'_',...
                'a2b',num2str(1/beta2alpha(1)),'_lowrx_data.mat'];
        elseif do_lbed_exclude,
                 fn=['Z:\elw\MATLAB\nfz_analy\meta_data\NFZ_',...
                structures{m},'_',...
                toxicities{i},'_',...
                'a2b',num2str(1/beta2alpha(1)),'_lbed_data.mat'];
        elseif do_late_exclude,
    	fn=['Z:\elw\MATLAB\nfz_analy\meta_data\NFZ_',...
                structures{m},'_',...
                toxicities{i},'_',...
                'a2b',num2str(1/beta2alpha(1)),'_acute_data.mat'];
            
         else
                
        %
        %tmp
        %
        
        %fn=['C:\Documents and Settings\williae1\nfz_meta_data\NFZ_',...
        fn=['Z:\elw\MATLAB\nfz_analy\meta_data\NFZ_',... 
            structures{m},'_',...
            toxicities{i},'_',...
            'a2b',num2str(1/beta2alpha(1)),'_data.mat'];
        end
        if save_result
            save(fn,'CGobj_org','CGstrct');
            disp(['Saving ',fn]);
        end
    end
end
end
