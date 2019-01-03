function [ADC_PDE_formulation, elapsed_time] ...
    = deff_PDE_formulation(gdir,sdeltavec,bdeltavec,seqvec,npervec,rtol,atol,...
    Ncmpt,Pts_cmpt_reorder,Ele_cmpt_reorder,DIFF_cmpts,Nboundary,Fac_boundary_reorder)

% diffusion equation (zero IC) to get the time-dependent diffusion coefficient

global FEM_M FEM_K FEM_A FEM_Q FEM_G
global QVAL UG
global BDELTA SDELTA SEQ OGSEPER

disp(['In PDE formulation']);


UG = gdir';
UG = UG/norm(UG);


is_replace_pde = 1;

if is_replace_pde==0
    [model_FEM_matrices] = deff_PDE_formulation_FEMat(Ncmpt,Pts_cmpt_reorder,Ele_cmpt_reorder,DIFF_cmpts);
end;

for icmpt = 1:Ncmpt
    [VOL(icmpt)] ...
        = get_volume_mesh(Pts_cmpt_reorder{icmpt},Ele_cmpt_reorder{icmpt});
end

nexperi = length(sdeltavec);
elapsed_time=zeros(Ncmpt, nexperi);


for icmpt = 1:Ncmpt
    
    if (is_replace_pde==0)
        % The 6 FE matrices can be obtained in the following way.
        FEM_M = model_FEM_matrices{icmpt}.M;
        FEM_K = model_FEM_matrices{icmpt}.K;
        FEM_A = model_FEM_matrices{icmpt}.A;
        FEM_Q = model_FEM_matrices{icmpt}.Q;
        FEM_G = DIFF_cmpts(icmpt)*model_FEM_matrices{icmpt}.G;
    else 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % To replace the pde matrices
        coordinates = Pts_cmpt_reorder{icmpt}; 
        elements = Ele_cmpt_reorder{icmpt};    
        [FEM_K,volumes]=stiffness_matrixP1_3D(elements',coordinates',DIFF_cmpts(icmpt));
        FEM_M=mass_matrixP1_3D(elements',volumes);   
        FEM_G=sparse(size(FEM_K,1),1);
        FEM_A=sparse(size(FEM_K,1),size(FEM_K,2));
        FEM_Q=sparse(size(FEM_K,1),size(FEM_K,2)); 
        coordinates_t=coordinates';
        one = sparse(size(FEM_K,1),1);
        one(:,:) = 1;

        for iboundary = 1:Nboundary
            neumann = Fac_boundary_reorder{icmpt}{iboundary}';
            if sum(size(neumann))>0
                coordVERTICES=zeros(size(neumann,1),3,3);
                coordVERTICES(1:size(neumann,1),1:3,1)=coordinates_t(neumann(:,1),:);
                coordVERTICES(1:size(neumann,1),1:3,2)=coordinates_t(neumann(:,2),:);
                coordVERTICES(1:size(neumann,1),1:3,3)=coordinates_t(neumann(:,3),:);
                [coordNORMALSnew,coordVERTICESnew] = COMPUTE_mesh_normals(coordVERTICES);

                ALL_VERTICES(1:size(elements,2),1:3,1)=coordinates_t(elements(1,:),:);
                ALL_VERTICES(1:size(elements,2),1:3,2)=coordinates_t(elements(2,:),:);
                ALL_VERTICES(1:size(elements,2),1:3,3)=coordinates_t(elements(3,:),:);        
                ALL_VERTICES(1:size(elements,2),1:3,4)=coordinates_t(elements(4,:),:);        

                midpoint_facet_coords = mean(coordVERTICES(:,:,:),3);
                midpoint_ele_coords = mean(ALL_VERTICES(:,:,:),3);
                for ifacet=1:size(midpoint_facet_coords,1)
                  % find the ending point of the normal of each facet
                  endpoint_of_normal = midpoint_facet_coords(ifacet,:) + coordNORMALSnew(ifacet,:);
                  % compare the position of the ending point of the normal to the
                  % facet
                  endpoint_eval =   coordNORMALSnew(ifacet,1)*(endpoint_of_normal(1)-midpoint_facet_coords(ifacet,1))+...
                                    coordNORMALSnew(ifacet,2)*(endpoint_of_normal(2)-midpoint_facet_coords(ifacet,2))+...
                                    coordNORMALSnew(ifacet,3)*(endpoint_of_normal(3)-midpoint_facet_coords(ifacet,3));

                  % find the closest element to the ending point of the normal  
                  d=sum((midpoint_ele_coords-ones(size(midpoint_ele_coords,1),1)*midpoint_facet_coords(ifacet,:)).^2,2);
                  [mval,mid]=min(d);
                  closestpoint_to_normal = midpoint_ele_coords(mid,:);

                  % compare the position of the closest point to the facet
                  closestpoint_eval = coordNORMALSnew(ifacet,1)*(closestpoint_to_normal(1)-midpoint_facet_coords(ifacet,1))+...
                                      coordNORMALSnew(ifacet,2)*(closestpoint_to_normal(2)-midpoint_facet_coords(ifacet,2))+...
                                      coordNORMALSnew(ifacet,3)*(closestpoint_to_normal(3)-midpoint_facet_coords(ifacet,3));
                  normal_orientation=(endpoint_eval*closestpoint_eval);
                  coordNORMALSnew(ifacet,:) =sign(-normal_orientation)*coordNORMALSnew(ifacet,:);
                end;
                if (1==0)
                    figure;
                    PLOT_3D_stl_patch(coordVERTICESnew,coordNORMALSnew)      
                end;
                mycoeff = (coordNORMALSnew(:,1)*UG(1) + coordNORMALSnew(:,2)*UG(2)+coordNORMALSnew(:,3)*UG(3));
                GG = flux_matrixP1_3D(neumann,coordinates', DIFF_cmpts(icmpt)*mycoeff);     
                FEM_G = FEM_G + GG*one;
            end;
        end;
% %         errorG=sum((MyG - FEM_G).^2);
% %         errorK=sum(sum((MyK - FEM_K).^2));
% %         errorM=sum(sum((MyM - FEM_M).^2));
% %         if errorG+errorK+errorM>1e-7
% %             disp('Matrices are different!'); 
% %             disp(['Error M:', num2str(errorM),', Error K:', num2str(errorK), ', Error G:', num2str(errorG)]);
% %     %         [full(MyG),  full(FEM_G)]
% %             stop;
% %         end;
% %         disp(['Error M:', num2str(errorM),', Error K:', num2str(errorK), ', Error G:', num2str(errorG)]);
        % To replace the pde matrices
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end;
    
    ODEsolve_atol = atol;
    ODEsolve_rtol = rtol;
    
    options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',...
        ODEsolve_rtol,'Vectorized','on','Stats','off',...
        'Jacobian',@odejac_bt_includeb);
    disp('DEFF PDE MODEL ***Uncoupled: start ode solver ode23t');
    nexperi = length(sdeltavec);  % Is it correct to be here???
    for iexperi = 1:nexperi
        iex_start_time = clock;
        SDELTA = sdeltavec(iexperi);
        BDELTA = bdeltavec(iexperi);
        SEQ = seqvec(iexperi);% for choosing case PGSE, OGSEcos or OGSEsin
        omega = 2*pi*npervec(iexperi)/SDELTA;
        OGSEPER = 1./omega*2*pi;%% set up number for OGSE
        QVAL = 0;
        TLIST = [0,SDELTA+BDELTA];
        ICC = zeros(size(FEM_M,1),1);
        sol = ode23t(@odefun_bt_includeb,TLIST,ICC,options);
       
        %deff_PDE_formulation_src{iexperi}{icmpt} = FEM_G.'*sol.y/VOL(icmpt)/VOL(icmpt)/DIFF_cmpts(icmpt);
        deff_PDE_formulation_src{iexperi}{icmpt} = FEM_G.'*sol.y/VOL(icmpt)/VOL(icmpt);

        deff_PDE_formulation_src_time{iexperi}{icmpt} = sol.x;

        hvec = deff_PDE_formulation_src{iexperi}{icmpt};
        tvec11 = deff_PDE_formulation_src_time{iexperi}{icmpt};
        Ftvec11 = seqintprofile(tvec11);
        a = trapz(tvec11,Ftvec11.*hvec*VOL(icmpt))/trapz(tvec11,Ftvec11.^2);
        
        ADC_PDE_formulation(icmpt,iexperi) = DIFF_cmpts(icmpt)-a;

        %ADC_PDE_formulation(icmpt,iexperi) = DIFF_cmpts(icmpt)*(1-a);
        %deff_PDE_formulation_src{iexperi}{icmpt} = FEM_G.'*sol.y/VOL(icmpt)/VOL(icmpt)/DIFF_cmpts(icmpt);
        elapsed_time(icmpt, iexperi) = etime(clock, iex_start_time);
    end
end
