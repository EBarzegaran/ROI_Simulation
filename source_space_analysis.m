% This script is for preliminary analysis and visualization of vernier-text
% scramble task in source space
% Elham Barzegaran, 5.22.2018

clear;clc;

PATH = '/Users/kohler/code/git/mrC';
addpath(genpath(PATH));
path2 = '/Users/kohler/code/git/gardner/mrTools/mrUtilities/mrFlatMesh';
addpath(genpath(path2));

addpath(genpath('/Users/kohler/code/git/sweepAnalysis/functions/helper'));

%% Convert spectrum data into source space
Path = '/Volumes/svndl/mrC_Projects/VernierLetters/source';
Inverse = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangROIsCorr.inv';
Inverse = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_wangkgsROIsCorr.inv';
%Inverse = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6.inv';
%Inverses = 'mneInv_bem_snr_10.inv';
[outData,~,FreqFeatures,subIDs] = mrC.SourceBrain(Path,Inverse,'domain','frequency','doSmooth' , true,'template','nl-0014');


%% surface average 
for cond = 1:10
    %outDataM{cond,1} = mean(cat(3,outData{cond,:})./repmat(max(abs(cat(3,outData{cond,:}))),[20484 1 1]),3)*10;% normalized over each condition
    %outDataM{cond,1} = mean(cat(3,outData{cond,:})./repmat(max(abs(cat(3,outData{10,:}))),[20484 1 1]),3)*10;% normalized over all condition
    datacond = cat(3,outData{cond,:});
%     cond
%     for Har = 1:2
%         for s = 1:size(datacond,1)
%             sourcedata = squeeze(datacond(s,(FreqFeatures.i1F1)+1,:));
%             %[ampDiff,phaseDiff,zSNR,errorEllipse] =fitErrorEllipse([real(sourcedata) imag(sourcedata)],[],1);
%             [results] = tSquaredFourierCoefs([real(sourcedata) imag(sourcedata)],[],0.01);
%             p1(s,Har)=results.pVal;
%             tsqrd1(s,Har)= results.tSqrd;
%         end
%     end
%     p{cond,1} = p1<0.0001;
%     tsqrd{cond,1}= tsqrd1;
    outDataM{cond,1} = mean(datacond,3);% not-normalized
end
%% plot amplitude and phase group level

anatDir = '/volumes/svndl/anatomy';
plottype = 'amp';%'amp';%'phase';
direction='ventral';
CondLabels = {'Text-1','Text-2','Text-3','Text-4','Text-5','Vernier-1','Vernier-2','Vernier-3','Vernier-4','Vernier-5'};

PlotSourceData(outDataM,FreqFeatures,subIDs(2),CondLabels,anatDir,plottype,direction);% parameters to add: Frequency of plots, color limits, view angles
PlotSourceData(outDataM,FreqFeatures,subIDs(2),CondLabels,anatDir,plottype,'anterior');% parameters to add: Frequency of plots, color limits, view angles
freq.i1F1 = 0;
%PlotSourceData(tsqrd,freq,subIDs(2),CondLabels,anatDir,plottype,'anterior');% parameters to add: Frequency of plots, color limits, view angles

%ROI plots
[~,RoiList] = mrC.Simulate.VisualizeSourceRoi('nl-0014',[],'wangkgs',[1 7 29 31 45],'anterior');


%% plot amplitude and phase individuals
% PlotSourceData(outData,FreqFeatures,subIDs,CondLabels,anatDir,plottype,'anterior');% parameters to add: Frequency of plots, color limits, view angles
 
%% Roi extraction: calculate amplitude and phase for each Roi
RoiType = 'wangkgs';
F1 = FreqFeatures.i1F1+1;
for sub = 1:numel(subIDs)
   RoiDir = fullfile(anatDir,subIDs{2},'Standard','meshes',[RoiType '_ROIs']); 
   [chunks{sub} roiList{sub}] = mrC.ChunkFromMesh(RoiDir,size(outData{1,sub},1));
   
   for cond = 1:size(outData,1)
       % Average the signal in each ROI
       roidata = chunks{sub}'*outData{cond,sub}./repmat(sum(chunks{sub})',[1 size(outData{cond,sub},2)]);
       %RoiData(sub,cond,:)= (roidata(:,F1))./abs(mean(roidata(:,[F1-2:F1-1 F1+1:F1+2]),2));
       RoiData(sub,cond,:)= roidata(:,F1);
       
       % Estimate the significance of each ROI using T-circ
       % Do we need to use another method rather than mean over ROI (like PCA)
       % Estimate the average phase of each ROI
       
   end
end 

%% plot group level results
% for sub = 1:numel(subIDs)
%     RData = RoiData (sub,:,:);
%     RoiData (sub,:,:) = RoiData (sub,:,:)./max(RData(:));
% end
sub = 1:6;%

CondLabels = {'Text-1','Text-2','Text-3','Text-4','Text-5','Ver-1','Ver-2','Ver-3','Ver-4','Ver-5'};

% plot amplitude and phase results
figure,
subplot(1,2,2),imagesc(squeeze(wrapTo2Pi(angle(nanmean(RoiData(sub,:,:)))))'*57);
h = colorbar; set(get(h,'title'),'string','Phase (degree)');
RoiLabels = (cellfun(@(x) x(5:end-4),roiList{1,1},'uni',false));%unique(cellfun(@(x) x(11:end-4),roiList{1},'uni',false));
set(gca,'YTick',1:numel(roiList{1,1}),'YTickLabel',RoiLabels);
set(gca,'XTick',1:numel(roiList{1,1}),'XTickLabel',CondLabels);
title('ROI Phase');

subplot(1,2,1),imagesc(squeeze(wrapTo2Pi(abs(nanmean(RoiData(sub,:,:)))))');
h = colorbar; set(get(h,'title'),'string','amplitude (µAmp/mm2)')
set(gca,'YTick',1:numel(roiList{1,1}),'YTickLabel',RoiLabels);
set(gca,'XTick',1:numel(roiList{1,1}),'XTickLabel',CondLabels);
title('ROI amplitude')

