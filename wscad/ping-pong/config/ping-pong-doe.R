require(DoE.base);

ping_pong_doe <- fac.design (
    nfactors = 2,
    replications=10,
    repeat.only=FALSE,
    blocks=1,
    randomize=TRUE,
    seed=72539,
    nlevels=c(3,21),
    factor.names=list(
        environment=c("docker","singularity","native"),
        sizeexp=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)));

export.design(ping_pong_doe, 
    path=".", 
    filename=NULL, 
    type="csv", 
    replace=TRUE,
    response.names=c("time"));
