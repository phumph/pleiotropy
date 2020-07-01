# TODO.md

## Before next update

- [x] Finish and run `call_adapteds.R`
- [x] Compile BCs via `combine_BCs_and_WGS.R`
- [x] Generate `t-SNE` cluster process
- [ ] Generate summary tables and plots for lineage count by cluster
- [ ] Map genotype to phenotype with mutual information

### Notes

Finishing `combine_BCs_and_WGS.R` script and I need to implement some logging in order to keep track of the counts at each of the steps:

* loading WGS
* loading barcodes from BFA files
* which WGS barcodes joined

Will implement the logging functions and write this a command-line executed script. One of the outputs will be a rendered markdown table with the results. In order for this to happen, I need to do the following:

- [x] Implement command-line arguments via `docopt`
- [ ] Write `src` file for the table template
- [ ] Initiate `.md` step-by-step for the analysis as `phumph.github.io/plt` page.

Thinking about the structure of this workflow.

I need to write `call_adapteds.R` to take the fitness files as input, as well as a couple of parameters, and output a file containing each of the barcodes and whether we consider it adapted or not. The script will handle each assay separately.

I'll need to setup some conditional routines to handle the auto-diploid flagging.. or perhaps I'll do that in a separate script.

Here's the file: [file:///Users/phumph/Dropbox/PLT/analysis_PLT/PLT/_obs/FIG_PLANNING_V1.nb.html](file:///Users/phumph/Dropbox/PLT/analysis_PLT/PLT/_obs/FIG_PLANNING_V1.nb.html)

Perhaps I need to compare output from average versus inverse-variance weighted averaging of the fitness between replicates? For now, I'll take the average type as an input variable and grab the correct columns accordingly. The next steps would be:

* calculate distribution without each barcode; calculate z-score
* sum up z-scores across environments for all neutral barcodes
* discard outliers
* re-calculate distribution with retained set
* calculate z-scores for all non-neutral set barcodes
* generate threshold for calling adapted

## Wednesday 13 November 2019

Still haven't made much progress. Need to perform outlier detection on neutral set and generate empirical distribution for generating contrasts.

- [x] Write `detect_outliers()` function to evaluate each fitness value of neutrals in a leave-one-out framework.
- [x] Plot z-scores of each neutral BC across environments
- [x] Decide on exclusion criteria for, e.g., average z-score.
  - [x] decide whether this should be signed or abs

OK so the simple theory here goes like this: the squared $z_i$ should behave as $\sim \chi^2(df = 1)$. So, if we calculate the joint likelihood of the data under the null as

$$
\mathcal{L}(x) = \prod^{n}_{j=1} \texttt{p.chisq}(x_{j},df=1)
$$

for each barcode, we should have a distribution of likelihoods we can play with and generate some kind of cutoff. Then we can flag neutrals as outliers and remove them, first plotting the retained versus excluded barcodes across environments.

Could find the maximum likelihood value of the estimator under the chi-sq function. OK I basically need to just write this all down and see what it looks like.

### New idea

Use random matrixes; generate parametric bootstrapped multivariate distribution and examine empirical distribution of $\delta$ from centroids of all individual barcodes. Identify outliers, then remove, then repeat to generate testing distribution for all barcodes.

Generate density estimates such that I can determine the probability under the distribution of each barcode in multivariate space.

The best way to think about this is to calculate the Mahalanobis distance:
> The Mahalanobis distance is the distance of the test point from the center of mass divided by the width of the ellipsoid in the direction of the test point.

A very nice description is on the wikipedia page [here](https://en.wikipedia.org/wiki/Mahalanobis_distance). Essentially, I need to write down a generative model of $k$-variate random vectors for the neutral set, standardize by environment, then use the transformed $k$-variant vectors to calculate the covariance matrix of the set. Using this, I can calculate the Mahalanobis distance of each component $k$-variate vector and the ellipsoid mean and generate probability estimates for each vector under this distribution.

Ok so here is the plan:

1. generate MVN $k$-variate random vectors for each barcode ($R=100000$) except for the leave-one-out barcode $q$.
2. center and scale each component vector $x_{ij} \rightarrow z_{ij}$.
3. calculate Mahalanobis distance of each component vector and the mean (which, to numerical error, should be all zeros).
4. Compare distance of hold-out vector versus distribution of others; get $p$-value for the barcode.

Check this site out in more detail: [https://rmflight.github.io/post/vignette-analysis/](https://rmflight.github.io/post/vignette-analysis/)
