require(DoE.base);

ep_DoE <- fac.design (
    nfactors = 2,
    replications=10,
    repeat.only=FALSE,
    blocks=1,
    randomize=TRUE,
    seed=72539,
    nlevels=c(3,4),
    factor.names=list(
        environment=c("docker","singularity","native"),
        parallelism=c(1,4,8,16)));

export.design(ep_DoE, 
    path=".", 
    filename=NULL, 
    type="csv", 
    replace=TRUE,
    response.names=c("time"));