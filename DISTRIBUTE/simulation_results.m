clear all;
format short

addpath SRC msh_files
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOM SRC/TETGEN

mat_files = {'spindle_soma_signal_d10_D43_3bvalues_180gdir.mat','spindle_dendrites_signal_d10_D43_3bvalues_180gdir.mat',...
            'pyramidal_soma_signal_d10_D43_3bvalues_180gdir.mat', 'pyramidal_soma_signal_d10_D43_3bvalues_180gdir.mat'};

for ifile = 1:4
    file_name = mat_files(ifile);
    load(file_name{1});
    nexperi = length(experi_btpde.sdeltavec);
    nb = size(experi_btpde.bvalues,2);
    for iexperi = 1:nexperi
        for ib = 2:nb
            S0 = sum((IC_cmpts.*VOL_cmpts));                 
            bv = experi_btpde.bvalues(iexperi,ib);
            title_str = ['BTPDE. All Cmpts. ','Experi ',...
               num2str(iexperi),', b= ',num2str(bv)];
            PLOT_HARDI_PT(points_gdir,squeeze(real(SIG_BTPDE_allcmpts_alldir(:,iexperi,ib)))/S0,title_str);
            PLOT_HARDI(points_gdir,squeeze(real(SIG_BTPDE_allcmpts_alldir(:,iexperi,ib)))/S0,title_str);
            figure;
            plot(all_diff(:, ib))
        end
    end
end