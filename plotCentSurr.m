function plotCentSurr(MSrc, roiChunk, SubIDs, List,direction)
for i = 1:12
    figure,
    subplot(2,2,1), mrC.Simulate.VisualizeSourceData(SubIDs{2},roiChunk(:,i),[],jmaColors('coolhot'),direction); 
    caxis([-1 1]);
    title([strrep(List(i),'_',' ') 'Original'],'fontsize',14);
    
    subplot(2,2,2),mrC.Simulate.VisualizeSourceData(SubIDs{2},MSrc(:,i),[],jmaColors('coolhot'),direction); 
    Data = MSrc(:,i);
    caxis([-max(Data) max(Data)]);
    title([strrep(List(i),'_',' ') 'Reconst'],'fontsize',14);
    
    subplot(2,2,3), mrC.Simulate.VisualizeSourceData(SubIDs{2},roiChunk(:,i+12),[],jmaColors('coolhot'),direction); 
    caxis([-1 1]);
    title([strrep(List(i+12),'_',' ') 'Original'],'fontsize',14);
    
    subplot(2,2,4),mrC.Simulate.VisualizeSourceData(SubIDs{2},MSrc(:,i+12),[],jmaColors('coolhot'),direction); 
    Data = MSrc(:,i+12);
    caxis([-max(Data) max(Data)]);
    title([strrep(List(i+12),'_',' ') 'Reconst'],'fontsize',14);
    
    set(gcf, 'Units', 'inches','PaperPositionMode','Manual', 'PaperPosition', [3, 3, 5, 8]);
    print(['/Users/babylab/Documents/Elham/text-scramble/Figs/CenterSurr/' List{i} '_' direction '.tif'],'-dtiff','-r300');close all
end
end