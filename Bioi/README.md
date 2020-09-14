# Bioi
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/Bioi)](https://cran.r-project.org/package=Bioi)
[![Travis-CI Build Status](https://travis-ci.org/zcolburn/Bioi.svg?branch=master)](https://travis-ci.org/zcolburn/Bioi)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/zcolburn/Bioi?branch=master&svg=true)](https://ci.appveyor.com/project/zcolburn/Bioi)
[![codecov](https://codecov.io/gh/zcolburn/Bioi/branch/master/graph/badge.svg)](https://codecov.io/gh/zcolburn/Bioi)
[![](https://cranlogs.r-pkg.org/badges/Bioi)](https://cran.r-project.org/package=Bioi)
[![DOI](https://zenodo.org/badge/110783607.svg)](https://zenodo.org/badge/latestdoi/110783607)


## Overview
Bioi is an R package containing implementations of solutions to common cell biology image processing problems. In particular, Bioi provides functions to perform connected component labeling on 1, 2, or 3 dimensional arrays, single linkage clustering on n-dimensional arrays, and identification of the points in one data set that are closest to each point in a second data set.


## Installation
The Bioi project can be found on its [GitHub repository](https://github.com/zcolburn/Bioi). It can be installed from that source using functionality available in the devtools package.


```r
install.packages("devtools")
library(devtools)
install_github("zcolburn/Bioi")
library(Bioi)
```


## Key functionality
### N-dimensional single linkage clustering
The objective of single linkage clustering is to place all points into groups such that all points within a group can be reached from any other point in the group by crossing bridges between points that are less than a critical separation distance. A small critical separation distance may result in a larger number of groups being identified. In contrast, a large critical separation distance may result in fewer groups being identified.


Using the function `euclidean_linker`, single linkage clustering can be performed in 1 or more dimensions. The function works in three modes: unpartitioned, partitioned, and parallelized. For small sample sizes (number of points less than `partition_req`) the unpartitioned method is used. For larger sample sizes the partitioned method is used.


The partitioned method works by iteratively dividing the data into smaller and smaller subsections until each partition contains fewer points than `partition_req`. The unpartitioned method is then used on each subsection before combining the data from each subsection.


The parallelized method works similarly to the partitioned method but operates in parallel. The number of cores to use can be specified by `num_cores`. To prevent too many threads from being generated `parallel_call_depth` can be specified. A higher depth will result in more threads being generated.


### Connected component labeling
A common image processing task is to group all connected "object-positive" pixels in an image into single groups. The connected component labeling function implemented here can be used on 1, 2, or 3-dimensional arrays representing 1, 2, or 3-dimensional "images". This functionality can be accessed using the `find_blobs` function.


### Determine the nearest neighbor in an alternate set of points
Photoactivated localization microscopy (PALM) data results in large numbers of protein localizations being identified. A common task when working with dual channel PALM data is to identify the distance separating points in one data set from points in a second data set. The function `find_min_dists` identifies the nearest neighbor to a point in a second data set and its distance from the point of interest.
