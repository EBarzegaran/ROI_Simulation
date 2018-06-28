function PlotSourceData(outData,FreqFeatures,subIDs,CondLabels,anatDir,plottype,direction)

    for sub = 1:numel(subIDs)
        figure,
        for cond= 1:size(outData,1)
            % phase data
            Data1 = angle(outData{cond,sub}(:,(FreqFeatures.i1F1)+1));
            Data1 = wrapTo2Pi(Data1);
            phasdata(:,cond,sub) = Data1;


            % amplitude data
            Data2 = abs(outData{cond,sub}(:,(FreqFeatures.i1F1)+1));
            ampdata(:,cond,sub) = Data2;

            % plot the results
            h = subplot(2,5,cond);
            switch plottype
                case 'phase'
                    mrC.Simulate.VisualizeSourceData(subIDs{sub},Data1,anatDir,jet(64),direction);
                    %colormap('jet');caxis([0 2*pi]);
                case 'amp'
                    mrC.Simulate.VisualizeSourceData(subIDs{sub},Data2,anatDir,jmaColors('hotcortex'),direction);
                    %mrC.Simulate.VisualizeSourceData(subIDs{sub},Data2,anatDir,pink(64),direction);
                    caxis([0 max(abs(outData{10,sub}(:)))/2]);
            end
            set(h,'position',get(h,'position')+[-.08 -.08 .03 .03]);
            axis tight
            title(CondLabels{cond},'fontweight','bold','fontsize',14);

        end
    end
end