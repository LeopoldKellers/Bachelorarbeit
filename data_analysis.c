/*
 * Ising Model Monte Carlo Simulation via Cluster Algorithm
 * Copyright (C) 2013  Leopold Kellers
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */


/*! @file data_analysis.c
 * @brief Contains the implementation of the stuff described in data_analysis.h.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "common.h"
#include "data_analysis.h"

binned_statistics* init_statistics(const int maxbinexp){
        assert(maxbinexp >= 0);
        binned_statistics* stat = xalloc(sizeof(binned_statistics));
        // +1 because maxbinexp = 0 is valid
        stat->S = xalloc((maxbinexp + 1) * sizeof(stat->S[0]));
        stat->n = 0; //No data in here so far
        stat->maxbinsize = 1 << maxbinexp;
        return stat;
}

void free_statistics(binned_statistics* stat) {
        free(stat->S);
        free(stat);
}

void update_statistics(binned_statistics* const stat, double* const restrict val){
        assert(NULL != stat);
        assert (NULL != val);
        double M = stat->M;
        double* const restrict S = stat->S;
        const int64_t n = stat->n;
        const int64_t maxbinsize = stat->maxbinsize;
        //The general idea is to walk throgh the input data (which is assumed
        //to have length equal to the size of the largest bin), and to update
        //the sum of squares S for each bin at each full bin length with
        //the average of the values in the bin. The update for M is 
        //straightforward, however, for binsizes > 1, it must be considered,
        //that the actual T is derived from T / binsize).
        int64_t k = (n * maxbinsize + 1);
        //first the easy case, which is maxbinsize 1
        if(0 == n) {
                M = val[0];
                S[0] = 0.0;
        }
        else{ 
                M += (val[0] - M) / k;
                const double delta = M - val[0];
                S[0] = S[0] + delta * delta * k / (k - 1);
        }
        //then the bin of size one for maxbinsize > 1
        for (int i = 1; i < maxbinsize; ++i) {
                ++k;
                M += (val[i] - M) / k;
                double delta = M - val[i];
                S[0] = S[0] + delta * delta * k / (k - 1);
                //and finally all the remaining bins
                int64_t binsize = 1;
                for (int l = 1; binsize < maxbinsize; ++l) {
                        binsize <<= 1;
                        // if i + 1 is dividable by the current bin size
                        // (works only for bin sizes which are poewers of 2!)
                        if (0 == ((i + 1) & (binsize - 1))) { 
                                // accumulate values
                                val[i] += val[i - (binsize >> 1)];
                                val[i] *= 0.5;
                                //compute the correct number of samples that
                                //have already been processed in this bin
                                //this does basically the same job as the k
                                //above, just a wee bit more fine-grained
                                int64_t j = n * (maxbinsize >> l);
                                j += (i >> l);
                                //check for initialization could be optimized by
                                //copy of paste this code into init_statistics()
                                //and omitting the check here, but this breaks
                                //the api and the current bottleneck is
                                //in the sampling, not the statistics...
                                if (j > 0) {
                                        //and update the sum of squares with a
                                        //method adapted from Youngs and Cramer
                                        delta = M - val[i];
                                        S[l] = S[l] + delta * delta * (j + 1) / (j);
                                        ++j;
                                }
                                else {
                                        //otherwise initialize the current bin
                                        S[l] = 0.0;
                                }
                        }
                }
        }
        stat->M = M;
        stat->n += 1;
}

//for the implementations below, keep in mind that stat->n only counts the
//number of times the update function has been called, and thus the correct
//number of elements processed for the bins with largest size. So that value
//has to be corrected for all other binsizes
double get_mean(const binned_statistics* const stat) {
        assert(NULL != stat);
        return stat->M;
}

double get_variance(const binned_statistics* const stat, const int binindex) {
        assert(NULL != stat);
        assert((1 << binindex) <= stat->maxbinsize);
        if (stat->n * (stat->maxbinsize >> binindex) > 1) {     
                return stat->S[binindex] / (stat->n * (stat->maxbinsize >> binindex) - 1);
        }
        else {
                return 0.0;
        }
}

double get_standard_error(const binned_statistics* const stat, const int binindex) {
        assert(NULL != stat);
        assert((1 << binindex) <= stat->maxbinsize);
        if (stat->n * (stat->maxbinsize >> binindex) > 1) {
                return sqrt(get_variance(stat, binindex) / (stat->n * (stat->maxbinsize >> binindex)));
        }
        else {
                return 0.0;
        }
}

