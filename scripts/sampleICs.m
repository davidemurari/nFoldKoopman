function ic=sampleICs(n_samples,n_pendula,plb,pub)
    %plb, pub : upper and lower bounds for the sampling of the momenta p

    theta = 2*pi*rand(n_samples,n_pendula);
    p = plb + (pub-plb)*rand(n_samples,n_pendula);

    ic = [theta,p];