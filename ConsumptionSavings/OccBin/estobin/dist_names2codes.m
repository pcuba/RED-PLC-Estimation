function codes = dist_names2codes(dist_names)

ndists = size(dist_names,1);
codes = zeros(ndists,1);


supported_dists = {'BETA_PDF','GAMMA_PDF','NORMAL_PDF','INV_GAMMA_PDF','UNIFORM_PDF'};


for i = 1:ndists
    codes(i) = find(strcmp(dist_names(i),supported_dists));
end

