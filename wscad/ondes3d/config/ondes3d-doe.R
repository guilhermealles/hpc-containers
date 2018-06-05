require(DoE.base);

ondes3d_DoE <- fac.design (
    nfactors = 3,
    replications=10,
    repeat.only=FALSE,
    blocks=1,
    randomize=TRUE,
    seed=72539,
    nlevels=c(3,3,3),
    factor.names=list(
        environment=c("docker","singularity","native"),
        context=c("mpi"),
        parallelism=c(1,4,8,16)));

export.design(ondes3d_DoE, 
    path=".", 
    filename=NULL, 
    type="csv", 
    replace=TRUE,
    response.names=c("time"));