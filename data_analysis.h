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


/*! @file data_analysis.h
 *  @brief Provides on the fly statistics with binning functionality to the
 *  other modules.
 *
 *  This module is designed to perform the data analysis part in the monte carlo
 *  simulation task. After initialization of the data structure provided below,
 *  there are a couple of functions to update the statistics with new values,
 *  retrieve the quantities mean, variance and standard error and free the
 *  initially allocated memory.
 *  The binning means, that a number (the binsize) in steps of 1, 2, 4, ... 2^N
 *  of input elements are averaged before they are processed in their respective
 *  bins. So assume N = 3. Then we have 4 bins of size 1, 2, 4 and 8.
 *  That means, when a set of 8 data values {1, 2, 3, 4, 5, 6, 7, 8} is put into
 *  the statistics, each bin will be updated like:
 *  bin 1: {1, 2, 3, 4, 5, 6, 7, 8}
 *  bin 2: {1.5,  3.5,  5.5,  7.5 }
 *  bin 3: {2.5,        6.5       }
 *  bin 4: {4.5                   }
 *  The method used is adapted from the method of Welford as recommended by
 *  Knuth in The Art of Computer Programming: Seminumerical Algorithms p. 232
 *  and also mentioned and analyzed with regards to numerical stability in
 *  Algorithms for Computing the Sample Variance: Analysis and Recommendations
 *  Tony F. Chan, Gene H. Golub and Randall J. LeVeque
 *  The American Statistician , Vol. 37, No. 3 (Aug., 1983), pp. 242-247 
 *  (http://www.jstor.org/stable/2683386)
 */
#ifndef DATA_ANALYSIS_H
#define DATA_ANALYSIS_H


/** The data structure that contains all the updateable statistical information.
 * Hand an instance of this to every function below, aside from
 * init statistics, which creates an instance of this.
 */
typedef struct Binned_statistics{
        int64_t n;          /**< The number of updates processed so far. For
                              the largest bin, this is equal to the number of
                              elements processed, the for the other bins, the
                              number of processed elements is equal to
                              n * maxbinsize / binsize */
        int64_t maxbinsize; /**< The size of the largest bin. At the same time,
                              this is the number of elements expected to be
                              handed to the update_statistics() function */
        double M;           /**< The sum of data points processed so far. */
        double* S;          /**< Array, that stores the sum of squares of all
                              elements that have been processed so far for
                              each bin.*/
} binned_statistics;


/** Create a new statistics instance with a specified number of bins.
 * The neccessary memory is allocated, initialized and returned as a pointer.
 * Exits the program in case no memory can be allocated!
 * @param maxbinexp The number of bins - 1,
 *      i.e. bins with sizes 1, 2, 4, ... 2^maxbinexp. Must be non-negative!
 *      A value of 0 obviously means no binning at all.
 * @return a pointer to the structure, that can be used in subsequent calls
 *      to the other functions below.
 */
binned_statistics* init_statistics(const int maxbinexp);

/** Free the memory that has been allocated by a previous call of
 * init_statistics().
 * @param stat the data structure to be freed. Must not be NULL!
 */
void free_statistics(binned_statistics* stat);

/** Update the statistics for all bins with the supplied values.
 *  The binning is carried out as in the module description above by averaging
 *  binsize consecutive values and using only these averages for the sum of
 *  squares.
 *  The number of elements is expected to be the 2^mabinexp, that had been
 *  supplied to the call of init_statistics() before upon creation of the stat
 *  structure.
 *  @param stat A valid pointer to a structure that has originally been returned
 *      by a call to init_statistics().
 *  @param val an array of 2^maxbinexp elements, to process now, where maxbinexp
 *      is the number of bins - 1, as had been supplied in the initial call to
 *      init_statistics(maxbinexp). The values in val may be mangled afterwards!
 */
void update_statistics(binned_statistics* stat, double* restrict val);

/** Get the mean of all the data processed so far.
 * The mean does obviously not depend on any binning.
 * @param stat A valid pointer to a structure that has originally been returned
 *      by a call to init_statistics().
 * @return the mean of the elements that have been processed so far.
 */
double get_mean(const binned_statistics* stat);

/** Get the variance of all the data processed so far for the specified bin
 * of size 2^binindex (That means numbering starts at 0 for no binning at all).
 * @param stat A valid pointer to a structure that has originally been returned
 *      by a call to init_statistics().
 * @param binindex The number of the bin of size 2^binindex to get the variance
 *      from. Must be positive and equal or less than maxbinexp, that has been
 *      supplied to the initial call of init_statistics(maxbinexp).
 * @return The variance of all the elements processed so far for the specified
 *      bin of size 2^binindex.
 */
double get_variance(const binned_statistics* stat, int binindex);
/** Get the standard error of all the data processed so far for the specified
 * bin of size 2^binindex (That means numbering starts at 0 for no binning at all).
 * @param stat A valid pointer to a structure that has originally been returned
 *      by a call to init_statistics().
 * @param binindex The number of the bin of size 2^binindex to get the standard
 *      error for. Must be positive and equal or less than maxbinexp, that has
 *      been supplied to the initial call of init_statistics(maxbinexp).
 * @return The standard error of all the elements processed so far for the 
 *      specified bin of size 2^binindex.
 */
double get_standard_error(const binned_statistics* stat, int binindex);
#endif // DATA_ANALYSIS_H
