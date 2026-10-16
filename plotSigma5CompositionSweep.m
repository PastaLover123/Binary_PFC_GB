function fig = plotSigma5CompositionSweep(outputDir)
% Replot without rerunning the simulations.
if nargin<1 || isempty(outputDir)
    outputDir=uigetdir(pwd,'Select a sigma5_composition_sweep folder');
    if isequal(outputDir,0), fig=[]; return; end
end
S=load(fullfile(outputDir,'composition_sweep.mat'),'cases');
cases=S.cases;
colors=lines(numel(cases));
fig=figure('Name','Sigma5 composition sweep','Color','w', ...
    'Position',[100 100 1250 850]);
layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
labels={'GB 1 displacement','GB 2 displacement','GB 1 speed','GB 2 speed'};
variables={'displacement1','displacement2','speed1','speed2'};
summary=zeros(0,7); notes=cell(0,1);
for panel=1:4
    ax=nexttile(layout); hold(ax,'on');
    for j=1:numel(cases)
        r=cases{j}; if isempty(r) || isempty(r.tracking), continue; end
        T=r.tracking; keep=T.timeSinceShear>=0;
        plot(ax,T.timeSinceShear(keep),T.(variables{panel})(keep), ...
            'Color',colors(j,:),'LineWidth',1.5, ...
            'DisplayName',sprintf('c = %.2f',r.composition));
    end
    xlabel(ax,'Time since shear start');
    if panel<=2, ylabel(ax,'Signed displacement');
    else, ylabel(ax,'|Normal velocity|');
    end
    title(ax,labels{panel}); grid(ax,'on');
    legend(ax,'show','Location','best');
end
for j=1:numel(cases)
    r=cases{j}; if isempty(r), continue; end
    v=[NaN NaN]; nValid=0; finalGamma=NaN;
    if ~isempty(r.tracking)
        T=r.tracking;
        loadingDuration=r.params.shear.maxStrain/abs(r.params.shear.rate);
        keep=T.gamma>=0.05 & T.gamma<=0.1 & T.timeSinceShear>=0 ...
            & T.timeSinceShear<=loadingDuration;
        nValid=nnz(keep & isfinite(T.displacement1) & isfinite(T.displacement2));
        finalGamma=T.gamma(end);
        for b=1:2
            d=T.(sprintf('displacement%d',b));
            ii=keep & isfinite(d);
            if nnz(ii)>=3
                fit=polyfit(T.timeSinceShear(ii),d(ii),1);
                v(b)=fit(1);
            end
        end
    end
    summary(end+1,:)=[r.composition v abs(v) nValid finalGamma];
    notes{end+1,1}=[r.status ': ' r.message];
end
summaryTable=array2table(summary,'VariableNames', ...
    {'composition','signedVelocity1','signedVelocity2','speed1','speed2', ...
    'validLateLoadingSamples','finalGamma'});
summaryTable.status=notes;
writetable(summaryTable,fullfile(outputDir,'velocity_summary.csv'));
title(layout,'Plain Sigma5 boundaries: same box and imposed shear rate');
savefig(fig,fullfile(outputDir,'composition_comparison.fig'));
exportgraphics(fig,fullfile(outputDir,'composition_comparison.pdf'),'ContentType','vector');
save(fullfile(outputDir,'velocity_summary.mat'),'summaryTable');
end
