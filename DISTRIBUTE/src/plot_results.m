xmin = 0;
xmax=0;
ymin=0;
ymax=0;
zmin=0;
zmax=0;
for ict = 1:mm.Ncmpt
    xx=max(mm.Pts_cmpt_reorder{icmpt}(1,:));
    yx=max(mm.Pts_cmpt_reorder{icmpt}(2,:));
    zx=max(mm.Pts_cmpt_reorder{icmpt}(3,:));
    xn=min(mm.Pts_cmpt_reorder{icmpt}(1,:));
    yn=min(mm.Pts_cmpt_reorder{icmpt}(2,:));
    zn=min(mm.Pts_cmpt_reorder{icmpt}(3,:));
    xmin = min(xmin,xn);
    xmax = max(xmax,xx);
    ymin = min(ymin,yn);
    ymax = max(ymax,yx);
    zmin = min(zmin,zn);
    zmax = max(zmax,zx);
end
figure; hold on;
cmptvec = [1:ncell];
for ict = 1:length(cmptvec)
    icmpt = cmptvec(ict);
    Fac = [];
    for iboundary = 1:mm.Nboundary
        Fac = [Fac,mm.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    h = trisurf(Fac',mm.Pts_cmpt_reorder{icmpt}(1,:),mm.Pts_cmpt_reorder{icmpt}(2,:),...
        mm.Pts_cmpt_reorder{icmpt}(3,:),real(YOUT{end}{end}{icmpt}(:,end)));
    set(h,'facealpha',0.9);
    axis equal;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]); colorbar;
    view(3);
    title('Solution,last experi,last b-value');
end
figure; hold on;
cmptvec = [ncell+1:mm.Ncmpt-1];
for ict = 1:length(cmptvec)
    icmpt = cmptvec(ict);
    Fac = [];
    for iboundary = 1:mm.Nboundary
        Fac = [Fac,mm.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    h = trisurf(Fac',mm.Pts_cmpt_reorder{icmpt}(1,:),mm.Pts_cmpt_reorder{icmpt}(2,:),...
        mm.Pts_cmpt_reorder{icmpt}(3,:),real(YOUT{end}{end}{icmpt}(:,end)));
    set(h,'facealpha',0.9);
    axis equal;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]); colorbar;
    view(3);
    title('Solution,last experi,last b-value');
end

figure;
cmptvec = [mm.Ncmpt];
for ict = 1:length(cmptvec)
    icmpt = cmptvec(ict);
    Fac = [];
    for iboundary = 1:mm.Nboundary
        Fac = [Fac,mm.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    h = trisurf(Fac',mm.Pts_cmpt_reorder{icmpt}(1,:),mm.Pts_cmpt_reorder{icmpt}(2,:),...
        mm.Pts_cmpt_reorder{icmpt}(3,:),real(YOUT{end}{end}{icmpt}(:,end)));
    set(h,'facealpha',0.9);
    axis equal;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]); colorbar;
    view(3);
    title('Solution,last experi,last b-value');
end

figure; hold on
iplot = 0;
for iexperi = 1:nexperi   
    yvec = real(MF_allcmpts(iexperi,:));
    bvec = bvalues(iexperi,:);
    h = plot(bvec, log10(yvec),...
        [colorvec_cell{1},markervec_cell{iexperi}]);
    set(h,'MarkerSize', 10, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = [' Experi ',mynum2str(iexperi)];
end
yvec = Sig_free*sum(IC_cmpts.*VOL);
h = plot(bvalues(:), log10(yvec),...
    [colorvec_cell{2},'','-']);
set(h,'MarkerSize', 10, 'LineWidth',1);
iplot = iplot + 1;
legend_vec{iplot} = ['free diffusion'];
for iexperi = 1:nexperi
    bvec = bvalues(iexperi,:);
    yvec = ADC_allcmpts_S0(iexperi)*exp(-ADC_allcmpts(iexperi)*bvec);
    h = plot(bvec, log10(yvec),...
        [colorvec_cell{1},'-']);
    set(h,'MarkerSize', 10, 'LineWidth',1);
end
legend(legend_vec{1:iplot},'Location','NorthEastOutside');
set(legend, 'FontSize',10)
set(gca, 'FontSize',10)
xlabel('bvalue')
ylabel('log10(Sig)')
title(['UG = [',mynum2str(UG),']']);

figure;
for iexperi = 1:nexperi
    subplot(nexperi,3,(iexperi-1)*3+1); hold on;
    bar(1:mm.Ncmpt,[ADC(:,iexperi)],'b');
    bar(mm.Ncmpt+1,ADC_allcmpts(iexperi,1),'r');
    title('ADC BT');
    set(gca,'ylim',[0,max(DIFF_cmpts)]);
    set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
    grid on;
end

for iexperi = 1:nexperi
    subplot(nexperi,3,(iexperi-1)*3+2); hold on;
    bar(1:mm.Ncmpt,ADC_PDE_formulation(:,iexperi),'b');
    bar(mm.Ncmpt+1,ADC_PDE_allcmpts(iexperi,1),'r');
    title('ADC DE');
    set(gca,'ylim',[0,max(DIFF_cmpts)]);
    set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
    grid on;
end

for iexperi = 1:nexperi
    subplot(nexperi,3,(iexperi-1)*3+3); hold on;
    bar(1:mm.Ncmpt,ADC_STA(:,iexperi),'b');
    bar(mm.Ncmpt+1,ADC_STA_allcmpts(iexperi,1),'r');
    title('ADC STA');
    set(gca,'ylim',[0,max(DIFF_cmpts)]);
    set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
    grid on;
end

figure; 
subplot(2,3,1);
hold all;
boundary_mat_plot=zeros(size(boundary_mat));
boundary_mat_plot(Cell_cmpt,:)=boundary_mat(Cell_cmpt,:);
spy(boundary_mat_plot,'b');

boundary_mat_plot(:,:)=0;
boundary_mat_plot(Box_cmpt,:)=boundary_mat(Box_cmpt,:);
spy(boundary_mat_plot,'r');

boundary_mat_plot(:,:)=0;
boundary_mat_plot(Nucleus_cmpt,:)=boundary_mat(Nucleus_cmpt,:);
spy(boundary_mat_plot,'k'); 

xlabel('iboundary'); ylabel('icmpt');
title('Connections boundary-compartment');
set(gca,'Ytick',[1:mm.Ncmpt]);
set(gca,'Xtick',[1:mm.Nboundary]);
grid on;

subplot(2,3,2); hold on;
bar(1:mm.Ncmpt,DIFF_cmpts,'b');
bar(mm.Ncmpt+1,ADC_free_allcmpts,'r');
title('DIFF free');
set(gca,'ylim',[0,max(DIFF_cmpts)]);
set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
set(gca,'Xtick',[1:mm.Ncmpt+1]);
xlabel('icmpt');
grid on;

subplot(2,3,3); hold on;
bar(1:mm.Nboundary,kappa_vec,'b');
title('Permeability');
set(gca,'ylim',[0,max(kappa_vec)]);
set(gca,'Ytick',linspace(0,max(kappa_vec),6));
set(gca,'Xtick',[1:mm.Nboundary]);
xlabel('iboundary');
grid on;

subplot(2,3,4); hold on;
bar(1:mm.Ncmpt,VOL,'b');
bar(mm.Ncmpt+1,VOL_allcmpts,'r');
title('VOL');
set(gca,'Xtick',[1:mm.Ncmpt+1]);
xlabel('icmpt');
grid on;

subplot(2,3,5); hold on;
bar(1:mm.Ncmpt,SA,'b');
bar(mm.Ncmpt+1,sum(SA),'r');
title('Surface Area');
set(gca,'Xtick',[1:mm.Ncmpt+1]);
xlabel('icmpt');
grid on;

subplot(2,3,6); hold on;
bar(1:mm.Ncmpt,SAu,'b');
bar(mm.Ncmpt+1,sum(SAu),'r');
title('SA in U_g');
set(gca,'Xtick',[1:mm.Ncmpt+1]);
xlabel('icmpt');
grid on;


figure;
for iexperi = 1:nexperi
    subplot(nexperi,2,iexperi); hold on;
    % Column 1: for deff_PDE
    % Column 2: for BT
    bar(1:mm.Ncmpt,[deff_PDE_elapsed_time(:,iexperi)],'b');
    title(['Timing DIFF, Experiment ',num2str(iexperi)]);
    ylabel('Timing (s)');
    xlabel('Cmpt. Index');
    
    subplot(nexperi, 2,iexperi+1); hold on;
    bar(1:length(bvec),[solve_mag_elapsed_time(:,iexperi)],'b');
    title(['Timing BT, Experiment ',num2str(iexperi)]);
    ylabel('Timing (s)');
    xlabel('b-val. Index');
end


figure; hold on
iplot = 0;
for iexperi = 1:nexperi   
    signal = real(MF_allcmpts(iexperi,:)./M0_allcmpts(iexperi,:));
    bvec = bvalues(iexperi,:);
    h = plot(bvec, signal,...
        [colorvec_cell{1},['-',markervec_cell{iexperi}]]);
    set(h,'MarkerSize', 10, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = [' Experi ',mynum2str(iexperi)];
end


