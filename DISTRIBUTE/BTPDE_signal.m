function [grad_dir, signal_allcmpts] ...
    = BTPDE_signal(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts)

% copy and modify from BTPDE_HARDI

Ncmpt = length(DIFF_cmpts);
[points,C,v] = spheresurface_regularpoints(1,experi_btpde.ngdir_total);
ngdir_total = size(points,1);
ii = find(points(:,3) >= 0);
% negii
for j = 1:size(ii,1)
    for k = 1:ngdir_total
        if (norm(points(j,1:2)-points(k,1:2)) < 1e-10  && points(j,3)+points(k,3) < 1e-10)
            negii(ii(j)) = k;
        end
    end
end
jc = 0;
for ic = 1:size(C,1)
    jj = find(C(ic,1) == ii);
    kk = find(C(ic,2) == ii);
    ll = find(C(ic,3) == ii);
    if (~isempty(jj) & ~isempty(kk) & ~isempty(ll))
        Ckeep(jc+1,1:3) = C(ic,1:3);
        jc = jc+1;
    end
end
graddir_index = ii;
ndir = length(graddir_index);

signal_allcmpts = nan*ones(ngdir_total,length(experi_btpde.bvalues));
grad_dir = nan*ones(ngdir_total,3);
for idir = 1:ndir
    experi_btpde.gdir = points(graddir_index(idir),:)';
    [~,~,~,MF_allcmpts,~,~] ...
        = BTPDE(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts);
    signal_allcmpts(idir,:) = MF_allcmpts;
    grad_dir(idir,:) = points(graddir_index(idir),:);
end
signal_allcmpts = rmmissing(signal_allcmpts);
grad_dir = rmmissing(grad_dir);