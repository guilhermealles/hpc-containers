require(DoE.base);

ondes3d_DoE <- fac.design (
    nfactors = 2,
    replications=10,
    repeat.only=FALSE,
    blocks=1,
    randomize=TRUE,
    seed=72539,
    nlevels=c(2,3),
    factor.names=list(
        environment=c('singularity','native'),
        parallelism=c(64,128,256)));

export.design(ondes3d_DoE, 
    path=".", 
    filename=NULL, 
    type="csv", 
    replace=TRUE,
    response.names=c('time'));
