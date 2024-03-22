%% Analysis code
% Accompanying code for the publication 
% "Exploring the Link between Brain Topological Resilience and ...
% Cognitive Performance in the Context of Aging and Vascular ...
% Risk Factors: A Cross-Ethnicity Population-Based Study"

% Please cite our work if you use the code provided here

% Please NOTE:

% This MATLAB script is to reproduce the statistical analysis of this study. 
% Since the data of this study can be obtained after application, 
% the variables in this script are dummy variables. 
% This script is only for understanding the implementation method of 
% our statistical analysis results and cannot be run directly.
%% Calculate CI for Spearman correlation
mycorr = @(x1,x2) corr(x1,x2,'type','spearman');
nIterations = 10000;
%% Analysis code for Figure 2
% Variables:
% BTR_PRECISE: BTR of PRECISE dataset. A vector of 2220*1.
% Age_PRECISE: Age of PRECISE dataset. A vector of 2220*1.
% MoCA_PRECISE: MoCA score of PRECISE dataset. A vector of 2220*1.
% BTR_MAS: BTR of MAS dataset. A vector of 264*1.
% Age_MAS: Age of MAS dataset. A vector of 264*1.
% Cog_MAS: Cognitive score of MAS dataset. A vector of 264*1.

% Spearman correlation for age and BTR 
[r,p] = corr(Age_PRECISE,BTR_PRECISE,'type','spearman');
[bci] = bootci(nIterations,{mycorr,Age_PRECISE,BTR_PRECISE});

[r,p] = corr(Age_MAS,BTR_MAS,'type','spearman');
[bci] = bootci(nIterations,{mycorr,Age_MAS,BTR_MAS});

% Spearman correlation for cognition and BTR
[r,p] = corr(MoCA_PRECISE,BTR_PRECISE,'type','spearman');
[bci] = bootci(nIterations,{mycorr,MoCA_PRECISE,BTR_PRECISE});

[r,p] = corr(Cog_MAS,BTR_MAS,'type','spearman');
[bci] = bootci(nIterations,{mycorr,Cog_MAS,BTR_MAS});
%% Analysis code for Figure 3
% Variables:
% BTR_PRECISE: BTR of PRECISE dataset. A vector of 2220*1.
% AS_PRECISE: Weighted atherosclerosis score of PRECISE dataset. A vector of 2220*1.
% VRF_PRECISE: Weighted vascular risk factor score of PRECISE dataset. A vector of 2220*1.
% BTR_MAS: BTR of MAS dataset. A vector of 264*1.
% VRF_MAS: Weighted vascular risk factor score of MAS dataset. A vector of 264*1.

% Spearman correlation for VRF and BTR 
[r,p] = corr(VRF_PRECISE,BTR_PRECISE,'type','spearman');
[bci] = bootci(nIterations,{mycorr,VRF_PRECISE,BTR_PRECISE});

[r,p] = corr(VRF_MAS,BTR_MAS,'type','spearman');
[bci] = bootci(nIterations,{mycorr,VRF_MAS,BTR_MAS});

% Spearman correlation for cognition and BTR
[r,p] = corr(AS_PRECISE,BTR_PRECISE,'type','spearman');
[bci] = bootci(nIterations,{mycorr,AS_PRECISE,BTR_PRECISE});

%% Analysis code for Figure S1
% Variables:
% Abnormal_edges: Number of abnormal edges. A vector of 2220*1.
% Edges: Number of edges. A vector of 2220*1.
figure;
subplot(1,2,1);histogram(Abnormal_edges(Abnormal_edges>0));
subplot(1,2,2);histogram(Edges);

%% Analysis code for Figure S3
% Variables:
% BTR_MAS: BTR of MAS dataset. A vector of 264*1.
% MMSE_MAS: MMSE score of MAS dataset. A vector of 264*1.

[r,p] = corr(MMSE_MAS,BTR_MAS,'type','spearman');
[bci] = bootci(nIterations,{mycorr,MMSE_MAS,BTR_MAS});

%% Analysis code for Figure S4
% Variables:
% BTR_PRECISE: BTR of PRECISE dataset. A vector of 2220*1.
% Age_PRECISE: Age of PRECISE dataset. A vector of 2220*1.
% Gender_PRECISE: MoCA score of PRECISE dataset. A vector of 2220*1.

[h,p,ci,stats] = ttest2(BTR_PRECISE(Age_PRECISE>=65),BTR_PRECISE(Age_PRECISE<65));
[h,p,ci,stats] = ttest2(BTR_PRECISE(Gender_PRECISE==0),BTR_PRECISE(Gender_PRECISE==1));

%% Analysis code for Figure S5 & S6
% Variables:
% BTR_PRECISE: BTR of PRECISE dataset. A vector of 2220*1.
% AS_PRECISE: Weighted atherosclerosis score of PRECISE dataset. A vector of 2220*1.
% VRF_PRECISE: Weighted vascular risk factor score of PRECISE dataset. A vector of 2220*1.
% Age_PRECISE: Age of PRECISE dataset. A vector of 2220*1.
% Gender_PRECISE: MoCA score of PRECISE dataset. A vector of 2220*1.

rloc=randperm(2220);
n0=100;k=0;cir=1;N=[];xes=[];xci=[];pv=[];
while cir
    locend=n0+20*k;k=k+1;N(k,:)=locend;
    if locend>2220
        locend=2220;
        cir=0;
    end
    %%%% choose one to analysze %%%%
    mdl=fitlm([AS_PRECISE(rloc(1:locend)),Age_PRECISE(rloc(1:locend)),Gender_PRECISE(rloc(1:locend))],BTR_PRECISE(rloc(1:locend)));
    %mdl=fitlm([VRF_PRECISE(rloc(1:locend)),Age_PRECISE(rloc(1:locend)),Gender_PRECISE(rloc(1:locend))],BTR_PRECISE(rloc(1:locend)));
    %%%%
    ci=coefCI(mdl);
    xes(k,:)=mdl.Coefficients.Estimate(2);
    xci(k,:)=ci(2,:);
    pv(k,:)=mdl.Coefficients.pValue(2);
end
figure;% CI & p-value
xN=[N;N(end:-1:1)];y=[xci(:,1);xci(end:-1:1,2)];
subplot(2,1,1);fill(xN,y,'r','FaceColor','#99DAEE','EdgeColor','#5A90DA','LineStyle','-','LineWidth',0.5);
xlim([N(1) 2220]);
hold on;plot([N(1) N(end)],[0 0],'Color','#5A90DA','MarkerFaceColor','#5A90DA','LineStyle','--','LineWidth',2);
plot(N,xes,'Marker','.','Color','#5A90DA','MarkerSize',20,'MarkerFaceColor','#5A90DA');xlim([N(1) 2220]);
subplot(2,1,2);plot(N,(pv),'Marker','.','MarkerSize',20,'MarkerFaceColor','#5A90DA');
hold on;plot([N(1) N(end)],[0.05 0.05],'Color','#5A90DA','LineStyle','--','LineWidth',2);
hold on;plot([N(1) N(end)],[0.01 0.01],'Color','#5A90DA','LineStyle','--','LineWidth',2);

figure;% CPS
subplot(2,1,1);
plot(N,xes,'Marker','.','MarkerSize',20,'MarkerFaceColor','#5A90DA');xlim([N(1) 2220]);
subplot(2,1,2);plot(N,(pv),'Marker','.','MarkerSize',20,'MarkerFaceColor','#5A90DA');
hold on;plot([N(1) N(end)],[0.05 0.05],'Color','#5A90DA','LineStyle','--','LineWidth',2);
hold on;plot([N(1) N(end)],[0.01 0.01],'Color','#5A90DA','LineStyle','--','LineWidth',2);

%% Analysis code for Figure S7
% Variables:
% BTR_PRECISE: BTR of PRECISE dataset. A vector of 2220*1.
% AS_PRECISE: Weighted atherosclerosis score of PRECISE dataset. A vector of 2220*1.
% VRF_PRECISE: Weighted vascular risk factor score of PRECISE dataset. A vector of 2220*1.
% Age_PRECISE: Age of PRECISE dataset. A vector of 2220*1.
% Gender_PRECISE: MoCA score of PRECISE dataset. A vector of 2220*1.

xes=[];pv=[];
for i=1:1000
	rloc=randperm(2220);
    n0=100;k=0;cir=1;N=[];
    while cir
        locend=n0+100*k;k=k+1;N(k,:)=locend;
        if locend>2220
            locend=2220;
            cir=0;
        end
        %%%% choose one to analysze %%%%
        mdl=fitlm(AS_PRECISE(rloc(1:locend)),BTR_PRECISE(rloc(1:locend)));
        %mdl=fitlm(VRF_PRECISE(rloc(1:locend)),BTR_PRECISE(rloc(1:locend)));
        %%%%
        ci=coefCI(mdl);
        xes(k,i)=mdl.Coefficients.Estimate(2);
        pv(k,i)=mdl.Coefficients.pValue(2);
    end
end
[hh, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pv(1:end-1,:),0.05);
figure;
subplot(2,1,1);boxplot((xes(1:end-1,:)'),'Symbol','');
xticklabels({'100','200','300','400','500','600','700','800','900','1000',...
    '1100','1200','1300','1400','1500','1600','1700','1800','1900','2000','2100','2200'});
subplot(2,1,2);boxplot((adj_p)','Symbol','');
xticklabels({'100','200','300','400','500','600','700','800','900','1000',...
    '1100','1200','1300','1400','1500','1600','1700','1800','1900','2000','2100','2200'});
hold on;plot(xlim,[0.05 0.05],'Color','#5A90DA','LineStyle','--','LineWidth',1);
hold on;plot(xlim,[0.01 0.01],'Color','#5A90DA','LineStyle','--','LineWidth',1);

%% Analysis code for Figure S10
% Variables:
% BTR_PRECISE_AAL3: BTR of PRECISE dataset for AAL-3 atlas. A vector of 2220*1.
% BTR_PRECISE_Schaefer17: BTR of PRECISE dataset for Schaefer17 atlas. A vector of 2220*1.
% Age_PRECISE: Age of PRECISE dataset. A vector of 2220*1.
% MoCA_PRECISE: MoCA score of PRECISE dataset. A vector of 2220*1.
% AS_PRECISE: Weighted atherosclerosis score of PRECISE dataset. A vector of 2220*1.
% VRF_PRECISE: Weighted vascular risk factor score of PRECISE dataset. A vector of 2220*1.

% Spearman correlation for age and BTR 
[r,p] = corr(Age_PRECISE,BTR_PRECISE_AAL3,'type','spearman');
[bci] = bootci(nIterations,{mycorr,Age_PRECISE,BTR_PRECISE_AAL3});
[r,p] = corr(Age_PRECISE,BTR_PRECISE_Schaefer17,'type','spearman');
[bci] = bootci(nIterations,{mycorr,Age_PRECISE,BTR_PRECISE_Schaefer17});

[r,p] = corr(MoCA_PRECISE,BTR_PRECISE_AAL3,'type','spearman');
[bci] = bootci(nIterations,{mycorr,MoCA_PRECISE,BTR_PRECISE_AAL3});
[r,p] = corr(MoCA_PRECISE,BTR_PRECISE_Schaefer17,'type','spearman');
[bci] = bootci(nIterations,{mycorr,MoCA_PRECISE,BTR_PRECISE_Schaefer17});

[r,p] = corr(AS_PRECISE,BTR_PRECISE_AAL3,'type','spearman');
[bci] = bootci(nIterations,{mycorr,AS_PRECISE,BTR_PRECISE_AAL3});
[r,p] = corr(AS_PRECISE,BTR_PRECISE_Schaefer17,'type','spearman');
[bci] = bootci(nIterations,{mycorr,AS_PRECISE,BTR_PRECISE_Schaefer17});

[r,p] = corr(VRF_PRECISE,BTR_PRECISE_AAL3,'type','spearman');
[bci] = bootci(nIterations,{mycorr,VRF_PRECISE,BTR_PRECISE_AAL3});
[r,p] = corr(VRF_PRECISE,BTR_PRECISE_Schaefer17,'type','spearman');
[bci] = bootci(nIterations,{mycorr,VRF_PRECISE,BTR_PRECISE_Schaefer17});
%% Analysis code for Figure S11
% Variables:
% BTR_MAS_AAL3: BTR of MAS dataset for AAL-3 atlas. A vector of 264*1.
% BTR_MAS_Schaefer17: BTR of MAS dataset for Schaefer17 atlas. A vector of 264*1.
% Age_MAS: Age of MAS dataset. A vector of 264*1.
% MoCA_MAS: MoCA score of MAS dataset. A vector of 264*1.
% VRF_MAS: Weighted vascular risk factor score of MAS dataset. A vector of 264*1.

% Spearman correlation for age and BTR 
[r,p] = corr(Age_MAS,BTR_MAS_AAL3,'type','spearman');
[bci] = bootci(nIterations,{mycorr,Age_MAS,BTR_MAS_AAL3});
[r,p] = corr(Age_MAS,BTR_MAS_Schaefer17,'type','spearman');
[bci] = bootci(nIterations,{mycorr,Age_MAS,BTR_MAS_Schaefer17});

[r,p] = corr(MoCA_MAS,BTR_MAS_AAL3,'type','spearman');
[bci] = bootci(nIterations,{mycorr,MoCA_MAS,BTR_MAS_AAL3});
[r,p] = corr(MoCA_MAS,BTR_MAS_Schaefer17,'type','spearman');
[bci] = bootci(nIterations,{mycorr,MoCA_MAS,BTR_MAS_Schaefer17});

[r,p] = corr(VRF_MAS,BTR_MAS_AAL3,'type','spearman');
[bci] = bootci(nIterations,{mycorr,VRF_MAS,BTR_MAS_AAL3});
[r,p] = corr(VRF_MAS,BTR_MAS_Schaefer17,'type','spearman');
[bci] = bootci(nIterations,{mycorr,VRF_MAS,BTR_MAS_Schaefer17});
%% Analysis code for Table S4
% Variables:
% BTR_PRECISE: BTR of PRECISE dataset. A vector of 2220*1.
% MoCA_PRECISE: MoCA score of PRECISE dataset. A vector of 2220*1.
% AS_m_PRECISE: Atherosclerosis scores of PRECISE dataset. A matrix of 2220*7.
% VRF_m_PRECISE: Vascular risk factors of PRECISE dataset. A matrix of 2220*8.

% BTR_MAS: BTR of MAS dataset. A vector of 264*1.
% Cog_MAS: Cognitive score of MAS dataset. A vector of 264*1.
% VRF_m_MAS: Vascular risk factors of MAS dataset. A matrix of 264*8.

weit=[];rep=[];redof=[];reci=[];
for i=1:size(VRF_m_PRECISE,2)
    mdl=fitlm(VRF_m_PRECISE,MoCA_PRECISE);
    weit(i,:)=-mdl.Coefficients.Estimate(2);
    rep(i,:)=mdl.Coefficients.pValue(2);
    ci=coefCI(mdl);
    reci(i,:)=ci(2,:);
    redof(i,:)=mdl.DFE;
end
VRF_PRECISE=VRF_m_PRECISE*weit;

weit=[];rep=[];redof=[];reci=[];
for i=1:size(AS_m_PRECISE,2)
    mdl=fitlm(AS_m_PRECISE,MoCA_PRECISE);
    weit(i,:)=-mdl.Coefficients.Estimate(2);
    rep(i,:)=mdl.Coefficients.pValue(2);
    ci=coefCI(mdl);
    reci(i,:)=ci(2,:);
    redof(i,:)=mdl.DFE;
end
AS_PRECISE=AS_m_PRECISE*weit;

weit=[];rep=[];redof=[];reci=[];
for i=1:size(VRF_m_MAS,2)
    mdl=fitlm(VRF_m_MAS,Cog_MAS);
    weit(i,:)=-mdl.Coefficients.Estimate(2);
    rep(i,:)=mdl.Coefficients.pValue(2);
    ci=coefCI(mdl);
    reci(i,:)=ci(2,:);
    redof(i,:)=mdl.DFE;
end
VRF_MAS=VRF_m_MAS*weit;
%% Analysis code for Table S5
% Variables:
% Age_PRECISE: Age of PRECISE dataset. A vector of 2220*1.
% MoCA_PRECISE: MoCA score of PRECISE dataset. A vector of 2220*1.
% AS_PRECISE: Weighted atherosclerosis score of PRECISE dataset. A vector of 2220*1.
% VRF_PRECISE: Weighted vascular risk factor score of PRECISE dataset. A vector of 2220*1.
% Age_MAS: Age of MAS dataset. A vector of 264*1.
% Cog_MAS: Cognitive score of MAS dataset. A vector of 264*1.
% MMSE_MAS: MMSE score of MAS dataset. A vector of 264*1.
% VRF_MAS: Weighted vascular risk factor score of MAS dataset. A vector of 264*1.

mdl=fitlm(MoCA_PRECISE,Age_PRECISE);
mdl=fitlm(AS_PRECISE,Age_PRECISE);
mdl=fitlm(VRF_PRECISE,Age_PRECISE);
mdl=fitlm(Cog_MAS,Age_MAS);
mdl=fitlm(MMSE_MAS,Age_MAS);
mdl=fitlm(VRF_MAS,Age_MAS);

%% Analysis code for Table S6
% Variables:
% BTR_PRECISE: BTR of PRECISE dataset. A vector of 2220*1.
% Age_PRECISE: Age of PRECISE dataset. A vector of 2220*1.
% Gender_PRECISE: Gender of PRECISE dataset. A vector of 2220*1.
mdl=fitlm(Age_PRECISE,BTR_PRECISE);
mdl=fitlm(Gender_PRECISE,BTR_PRECISE);
mdl=fitlm([Age_PRECISE,Gender_PRECISE],BTR_PRECISE,'y ~ x1 * x2');

%% Analysis code for Table S7
% Variables:
% p_values: P-values of all analysis in this study.
% you should add this toolbox in your matlab path:
% https://ww2.mathworks.cn/matlabcentral/fileexchange/27418-fdr_bh
[hh,crit_p,adj_ci_cvrg,adj_p]=fdr_bh(p_values,0.05);
