dim =3;
sigma_0 = .0001;
cluster_num = 2;
cluster_size = 1E3;
X_clusters = [];
col_clusters = repmat(1:cluster_num,cluster_size,1);
col_clusters = col_clusters(:);
for i =1:cluster_num,
    X_clusters = [X_clusters; mvnrnd(rand(dim,1)', diag(sigma_0*ones(dim,1)), cluster_size)];
end
N_clusters = size(X_clusters,1)
figure(19988)
scatter3(X_clusters(:,1), X_clusters(:,2), X_clusters(:,3),10, col_clusters, 'filled');
title('Original Data')


 %%
clear opts
opts.stop_lying_iter = 2E2;
opts.max_iter = 400;
opts.start_late_exag_iter = 300;
opts.late_exag_coeff = 6;
opts.rand_seed = -1;
cluster_firstphase = fast_tsne(X_clusters,opts);
figure(433)
scatter(cluster_firstphase(:,1), cluster_firstphase(:,2),  10, col_clusters, 'filled');
title('After t-SNE')
colormap(lines(cluster_num))





%% Example using user-provided P
dim =3;
sigma_0 = .0001;
cluster_num = 2;
cluster_size = 500;

% Since we will provide the P, the X is entirely unnecessary. Making it random here.
X_clusters = rand(cluster_size*cluster_num,cluster_size*cluster_num);
col_clusters = repmat(1:cluster_num,cluster_size,1);
col_clusters = col_clusters(:);

%Making a block diagonal P and then add noise to it
P =blkdiag(rand(cluster_size,cluster_size), rand(cluster_size,cluster_size));
P(P<0.9) = 0;
P =P+ 0.1*rand(cluster_size*cluster_num,cluster_size*cluster_num);
P(P<0.08) = 0;
write_P(P)

%Run FIt-SNE
clear opts
opts.stop_lying_iter = 2E2;
opts.max_iter = 400;
opts.start_late_exag_iter = 300;
opts.late_exag_coeff = 6;
opts.rand_seed = -1;
opts.use_existing_P = 1;
cluster_firstphase = fast_tsne(X_clusters,opts);

%Plot result
figure(433)
scatter(cluster_firstphase(:,1), cluster_firstphase(:,2),  10, col_clusters, 'filled');
title('After t-SNE')
colormap(lines(cluster_num))
