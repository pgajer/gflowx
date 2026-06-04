\name{rdgraph-bridges}
\alias{rdgraph-bridges}
\alias{fit.rdgraph.regression}
\alias{refit.rdgraph.regression}
\alias{permutation.test.rdgraph}
\alias{perm.test.audit}
\alias{bayes.bootstrap.rdgraph}
\alias{lcor.with.posterior}
\alias{extremality.summary}
\alias{label.extremality.3d}
\alias{compute.pextrema.nbhds}
\alias{compute.cluster.summary}
\alias{extract.cluster.representatives}
\alias{compute.dbscan.cluster.summary}
\alias{extract.dbscan.cluster.representatives}
\alias{analyze.dbscan.noise}
\alias{subject.neighborhood.stats}
\alias{rdgraph.neighbor.weights}
\title{Bridge exports for rdgraph regression workflows}
\usage{
fit.rdgraph.regression(...)
refit.rdgraph.regression(...)
permutation.test.rdgraph(...)
perm.test.audit(...)
bayes.bootstrap.rdgraph(...)
lcor.with.posterior(...)
extremality.summary(...)
label.extremality.3d(...)
compute.pextrema.nbhds(...)
compute.cluster.summary(...)
extract.cluster.representatives(...)
compute.dbscan.cluster.summary(...)
extract.dbscan.cluster.representatives(...)
analyze.dbscan.noise(...)
subject.neighborhood.stats(...)
rdgraph.neighbor.weights(...)
}
\arguments{
\item{...}{Arguments forwarded unchanged to the corresponding \pkg{gflow}
function.}
}
\value{
The object returned by the delegated \pkg{gflow} implementation.
}
\description{
These functions are the initial \pkg{gflowx} landing surface for legacy rdgraph
regression workflows. They currently delegate to the implementations in
\pkg{gflow}; this keeps existing behavior fixed while the native rdgraph engine
and its graph utilities are extracted in smaller audited slices.
}
