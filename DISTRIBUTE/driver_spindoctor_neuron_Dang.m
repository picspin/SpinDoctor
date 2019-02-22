clear all;
format short

addpath SRC
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOM SRC/TETGEN


global BDELTA SDELTA SEQ OGSEPER
global PGSE OGSEsin OGSEcos

load msh_files/spindle/03b_spindle5aACC_nodes.txt;
load msh_files/spindle/03b_spindle5aACC_elements.txt;

fname_params_cells = 'params_cells2.in';
fname_params_simul_domain  = 'params_simul_domain2.in';
fname_params_simul_experi = 'params_simul_experi2.in';

[params_cells,fname_cells] = create_geom(fname_params_cells);

%PLOT_CELLS(params_cells.cell_shape,fname_cells);

[params_domain_geom,params_domain_pde,params_domain_femesh] ...
    = read_params_simul_domain(fname_params_simul_domain);

%PLOT_SURFACE_TRIANGULATION(params_cells.cell_shape,fname_cells,params_domain_geom);

[DIFF_cmpts,kappa_bdys,IC_cmpts,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index,Ncmpt,Nboundary] ...
    = PREPARE_PDE(params_cells.ncell,params_cells.cell_shape,params_domain_geom,params_domain_pde);

fname_tetgen = [fname_cells];
[fname_tetgen_femesh] = ...
    create_femesh_fromcells(params_cells,fname_cells,params_domain_geom,params_domain_femesh,fname_tetgen);



[mymesh,cmpts_bdys_mat] = read_tetgen(fname_tetgen_femesh,params_cells.para_deform,Ncmpt,Nboundary);

mymodel = createpde();
elements = X03b_spindle5aACC_elements;  
nodes = X03b_spindle5aACC_nodes;

geometryFromMesh(mymodel,nodes',elements');

sigma0 = DIFF_cmpts(1);
specifyCoefficients(mymodel,'m',0,'d',1,'c',sigma0,'a', 0,'f',0);
nface = mymodel.Geometry.NumFaces; %number of face
%     normal dot (c grad u) + q*u = g
for iface = 1:nface
    applyBoundaryCondition(mymodel,'neumann','face',iface,'g',1,'q',1,'Vectorized','on'); % 3-D geometry
end

mymodel_FEM_matrices = assembleFEMatrices(mymodel);

G = mymodel_FEM_matrices.G;
Q = mymodel_FEM_matrices.Q;

mymesh.Fac_boundary_reorder{1}{1} = [];
Ind_Nodes_Bdy = find(G ~= 0);  
mymesh.Pts_boundary_reorder{1}{1} = Ind_Nodes_Bdy;

nele = size(elements,1);
facets = [];
for iele = 1:nele
    ii1 = sum(Ind_Nodes_Bdy == elements(iele,1)); 
    ii2 = sum(Ind_Nodes_Bdy == elements(iele,2));
    ii3 = sum(Ind_Nodes_Bdy == elements(iele,3));
    ii4 = sum(Ind_Nodes_Bdy == elements(iele,4));
    if (ii1+ii2+ii3+ii4 == 3)
        jj = find([ii1,ii2,ii3,ii4]==1);
        facets = [facets; elements(iele,jj)];
    end
end


%%%% Replaced by Dang %%%%
[elems2faces, faces2nodes]=get_faces(elements);
faces2elems=entryInWhichRows(elems2faces); 
bindex=find(faces2elems(:,2)==0);           %index of boundary facets - belong to one element only
boundary=faces2nodes(bindex,:);             %nodes of boundary facets 
facets=sort(boundary,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%

mymesh.Ncmpt = 1;
mymesh.Nboundary = 1;

mymesh.Pts_cmpt_reorder{1} = mymodel.Mesh.Nodes;
mymesh.Ele_cmpt_reorder{1} = mymodel.Mesh.Elements;
mymesh.Pts_ind{1} = [1:size(mymodel.Mesh.Nodes,2)]';

mymesh.Pts_boundary_reorder{1}{1} = [1:size(mymodel.Mesh.Nodes,2)]';
mymesh.Fac_boundary_reorder{1}{1} = facets';



if (~isempty(mymesh))
    %PLOT_FEMESH(mymesh,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index);
    
    [experi_common,experi_hadc,experi_btpde] ...
        = read_params_simul_experi(fname_params_simul_experi);
    
    [VOL_cmpts,SA_cmpts,SAu_cmpts,VOL_allcmpts,VF_cmpts,SoV_cmpts] ...
        = GET_VOL_SA(mymesh,experi_common.gdir);
    
    %PLOT_GEOMETRY_INFO(cmpts_bdys_mat,OUT_cmpts_index,IN_cmpts_index,ECS_cmpts_index,VOL_cmpts,SA_cmpts,SAu_cmpts);
    
    
    
    if (experi_common.ngdir_total == 1)
        % BTPDE
        if (~isempty(experi_btpde))
            
            [TOUT,YOUT,MF_cmpts,MF_allcmpts,difftime,BTPDE_elapsed_time] ...
                = BTPDE(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts);
            
            PLOT_MAGNETIZATION(mymesh,YOUT,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index);
            
            [ADC_cmpts,ADC_allcmpts,ADC_allcmpts_S0] = FIT_SIGNAL(MF_cmpts,MF_allcmpts,experi_btpde.bvalues);
            
            [Sig_free,ADC_free_allcmpts] = ADCFREE(experi_btpde.bvalues,DIFF_cmpts,VOL_cmpts,IC_cmpts);
            
            PLOT_SIGNAL(experi_btpde.bvalues,MF_allcmpts,Sig_free,ADC_allcmpts_S0,ADC_allcmpts,'BTPDE');
            
            PLOT_ADC(ADC_cmpts,ADC_allcmpts,DIFF_cmpts,'BTPDE');
            
            
        end
        %HADC
        if (~isempty(experi_hadc))
            
            
            [ADC_HADC_cmpts,ADC_HADC_allcmpts,HADC_elapsed_time] = HADC(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts);            
            PLOT_ADC(ADC_HADC_cmpts,ADC_HADC_allcmpts,DIFF_cmpts,'HADC');
                        
            nexperi = length(experi_common.sdeltavec);
            for iexperi = 1:nexperi
                if (~isempty(experi_btpde))
                    bvec = experi_btpde.bvalues(iexperi,:);
                else
                    bvec = linspace(0,2000,5);
                end
                nb = length(bvec);
                for ib = 1:nb
                    for icmpt = 1:mymesh.Ncmpt
                        MF_HADC_cmpts(icmpt,iexperi,ib) = VOL_cmpts(icmpt)*IC_cmpts(icmpt)*exp(-ADC_HADC_cmpts(icmpt,iexperi)*bvec(ib));
                    end
                    MF_HADC_allcmpts(iexperi,ib) = 0;
                    for icmpt = 1:mymesh.Ncmpt
                        MF_HADC_allcmpts(iexperi,ib) = MF_HADC_allcmpts(iexperi,ib) + MF_HADC_cmpts(icmpt,iexperi,ib);
                    end
                end
            end  
            [Sig_free,ADC_free_allcmpts] = ADCFREE(bvec,DIFF_cmpts,VOL_cmpts,IC_cmpts);
 
            PLOT_SIGNAL(bvec,MF_HADC_allcmpts,Sig_free,sum(VOL_cmpts.*IC_cmpts)*ones(nexperi,1),ADC_HADC_allcmpts,'HADC');
        end
        [ADC_STA_cmpts,ADC_STA_allcmpts] = STA(experi_common,DIFF_cmpts,VOL_cmpts,SAu_cmpts,IC_cmpts);
        PLOT_ADC(ADC_STA_cmpts,ADC_STA_allcmpts,DIFF_cmpts,'STA');
        
    else
        % BTPDE
        if (~isempty(experi_btpde))
            [SH_points,ADC_BT_cmpts_alldir,ADC_BT_allcmpts_alldir] ...
                = BTPDE_HARDI(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts);
            PLOT_HARDI(SH_points, ADC_BT_allcmpts_alldir);
        end
        %HADC
        if (~isempty(experi_hadc))
            [points,ADC_HADC_cmpts_alldir,ADC_HADC_allcmpts_alldir] ...
                = HADC_HARDI(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts);
            PLOT_HARDI(points,ADC_HADC_allcmpts_alldir);
        end
    end
    
end