clear; clc;

addpath SRC
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOMETRY SRC/TETGEN SRC/UTILITIES

addpath COMPUTE_mesh_normals/COMPUTE_mesh_normals

mydefinitions;
ncolor = length(colorvec_cell);

SEQ_DEFINITIONS

fname_domain = 'InputFiles_Simulation/simulation_parameters_domain.in';
fname_experiment = 'InputFiles_Simulation/simulation_parameters_experiment.in';

[cell_shape,Rratio_nucleus,dcoeff_nucleus,dcoeff_cytoplasm,dcoeff_exterior,...
    ic_nucleus,ic_cytoplasm,ic_exterior,ic_llimit,ic_ulimit,kappa_nc,kappa_ce,include_box,box_gap,...
    create_geom,fname_geom,ncell,Hcyl,Rmean,Rmin,Rmax,Htetgen] ...
    = read_simulation_parameters_domain(fname_domain);

[gdir,bvalues,qvalues,sdeltavec,bdeltavec,seqvec,npervec,...
    rtol_bt,atol_bt,rtol_deff,atol_deff,const_q,tetgen_cmd] ...
    = read_simulation_parameters_experiment(fname_experiment);

if (include_box == 1)
    box_str = 'box';
else
    box_str = 'nobox';
end
if (cell_shape == 1)
    cell_shape_name = 'ellipses';
elseif (cell_shape == 2)
    cell_shape_name = 'cylinders';
end
if (create_geom == 0)
    fname = fname_geom;
else
    fname = [cell_shape_name,num2str(ncell),'_R',num2str(Rmean)];
end

fname_cells_description = ['InputFiles_Geometry/',fname,'_description.in'];

if (create_geom == 0)
else
    if (cell_shape == 1)
        create_ellipses_inputfile(ncell,Rmean,Rmax,Rmin,fname_cells_description);

    else
        create_cylinders_inputfile(ncell,Rmean,Rmax,Rmin,Hcyl,fname_cells_description);
    end
end

fname_tetgen = ['InputFiles_Tetgen/',fname,'_',box_str];

[fname_tetgen_femesh] = create_cells_femesh(fname_cells_description,fname_tetgen,...
    include_box,box_gap,Rratio_nucleus,cell_shape_name,Htetgen,tetgen_cmd);

[ADC_PDE_formulation,Ncmpt,DIFF_cmpts] ...
    = deff_PDE_formulation(fname_tetgen_femesh,ncell,Rratio_nucleus,...
    dcoeff_nucleus,dcoeff_cytoplasm,dcoeff_exterior,include_box,...
    gdir,sdeltavec,bdeltavec,seqvec,npervec,rtol_deff,atol_deff);

[TOUT,YOUT,MT,Ncmpt,Nboundary,Cell_cmpt,Box_cmpt,Nucleus_cmpt,difftime,...
    Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder,...
    DIFF_cmpts,IC_cmpts,UG] ...
    = solve_magnetization(fname_tetgen_femesh,ncell,Rratio_nucleus,...
    dcoeff_nucleus,dcoeff_cytoplasm,dcoeff_exterior,...
    ic_nucleus,ic_cytoplasm,ic_exterior,ic_llimit,ic_ulimit,kappa_nc,kappa_ce,include_box,...
    gdir,qvalues,sdeltavec,bdeltavec,seqvec,npervec,rtol_bt,atol_bt);

[ADC_allcmpts,ADC_allcmpts_polydeg,ADC_allcmpts_S0,Deff_STA_allcmpts,...
    ADC,ADC_polydeg,ADC_S0,Deff_STA,MF_allcmpts,M0_allcmpts,S0_allcmpts,...
    MF,M0,S0,VOL,SA,SAu,VOL_frac,SoV]...
    = post_processing(MT,bvalues,Ncmpt,Nboundary,sdeltavec,bdeltavec,seqvec,npervec,...
    Pts_cmpt_reorder,Ele_cmpt_reorder,Fac_boundary_reorder,...
    DIFF_cmpts,UG);

for icmpt = 1:Ncmpt
    figure; hold on;
    Fac = [];
    for iboundary = 1:Nboundary
        Fac = [Fac,Fac_boundary_reorder{icmpt}{iboundary}];
    end 
    h = trisurf(Fac',Pts_cmpt_reorder{icmpt}(1,:),Pts_cmpt_reorder{icmpt}(2,:),Pts_cmpt_reorder{icmpt}(3,:));
    set(h,'facealpha',0.1);
    axis equal;
    view(3);
end
title('Finite Element Mesh');

nexperi = length(difftime);
Sig_free = zeros(size(bvalues(:)));

vol = 0;
for icmpt = 1:Ncmpt
    Sig_free = Sig_free+IC_cmpts(1,icmpt)*VOL_frac(icmpt)*exp(-DIFF_cmpts(icmpt)*bvalues(:));
    vol = vol+IC_cmpts(1,icmpt)*VOL_frac(icmpt);
end
Sig_free = Sig_free/vol;


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
    subplot(nexperi,3,(iexperi-1)*3+1);
    bar(ADC(:,iexperi));
    title('ADC BT');
    set(gca,'ylim',[0,max(DIFF_cmpts)]);
    set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
    grid on;
end
for iexperi = 1:nexperi
    subplot(nexperi,3,(iexperi-1)*3+2);
    bar(ADC_PDE_formulation(:,iexperi));
    title('ADC DE');
    set(gca,'ylim',[0,max(DIFF_cmpts)]);
    set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
    grid on;
end
for iexperi = 1:nexperi
    subplot(nexperi,3,(iexperi-1)*3+3);
    bar(Deff_STA(:,iexperi));
    title('ADC STA');
    set(gca,'ylim',[0,max(DIFF_cmpts)]);
    set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
    grid on;
end





    




