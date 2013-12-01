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


/*! @file monte_carlo.h
 *  @brief Contains all neccessary declarations for two implementations of
 *  Markov chain monte carlo cluster algorithms for the Ising model, the
 *  Swendsen-Wang and the Wolff algorithm.
 *
 *  This module contains all the stuff neccessary to sample the two-point
 *  correlation function in the high temperature regime of the Ising Model in
 *  arbitrary dimension n supplied at compile time via -DD=n switch of the
 *  preprocessor. Also at compile time, there is a choice to select between the
 *  Mersenne Twister random number generator (-DRNG=1) or the (slower) WELL RNG
 *  (-DRNG=2) for making sure that the choice of random number generator does
 *  not make a difference.
 *
 *  The two cluster algorithms are implemented in a fundamentally different way,
 *  so the results can be checked against each other for consistency to make
 *  sure, that the code produces correct results.
 *  The implementation of the Swendsen Wang algorithm is according to 
 *  K. Muthy: Monte Carlo Methods in Statistical Physics p.25
 *  so using the Hoshen-Kopelman approach, which is basically a union-find
 *  algorithm for cluster labeling, while the Wolff algorithm is based on the
 *  pseudocode with the pocket-cluster method as described in
 *  W. Krauth: Statistical Mechanics - Algorithms and Computations p. 255.
 *  The Swendsen Wang method, which always computes all clusters in one sweep
 *  and flips some of them, averages over more values and thus has less variance
 *  in the samples of the correlation function, but is slower than the Wolff
 *  single cluster method, which only builds a single cluster and flips it,
 *  which is much faster, especially the higher the temperature is above the
 *  critical point, but also has a higher variance, so that the usefulness
 *  depends on the very application.
 */
#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#define UP ((int8_t) 1) /**< value for spin up on the lattice*/
#define DOWN ((int8_t) -1) /**< Opposite of UP*/

#define DEFAULT_L 6 /**<the default side-length of the lattice to be evaluated*/
#define DEFAULT_J 1.0 /**<the default coupling strength of neighbouring spins*/
#define DEFAULT_BETA 1.0 /**<the default inverse temperature to simulate at*/
#define DEFAULT_NSWEEPS 10000 /**<default number of sweeps to sample in a run*/
#define DEFAULT_EQUILIBRATIONTIME 1000 /**< The number of sweeps to drop before
                                         starting to sample*/

//choose the right random number generator
#ifndef RNG
#define RNG 1  /**<default random number generator - see Makefile for details*/
#warning "choosing dsfmt as default random number generator"
#endif
#if RNG == 1 
#include "dsfmt/dSFMT.h"
#elif RNG == 2
#include "well/WELL19937a.h"
#else 
#error "Unknown random number generator chosen - see Makefile for details"
#endif
#define SEEDLENGTH 624 /**<number of random values to initialize the rng with*/
#define DEFAULT_SEEDFILE "/dev/urandom" /**< of course only avaliable on linux*/
#define FLAGGED 1 /**< is site already in cluster or pocket? */
#define NOT_FLAGGED 0 /**< In this case, it is not */



#ifdef USE_SWENDSEN_WANG
/** Representation of data neccessary for the Swendsen Wang cluster algorithm.
 * This data structure is returned by init_swendsen_wang(). The functions
 * starting with swendsen_wang_* below operate on it.
 * Modelling the data in plain arrays is much more convenient, and probably
 * faster in this case, as the whole lattice has to be traversed on every sweep.
 */
typedef struct Swendsen_Wang {
        index_t N;               /**< The overall number of lattice sites to
                                   work on*/
        sidelength_t L0;         /**< The sidelength in correlation direction*/
        int8_t* spins;           /**< Store the lattice spin configuration*/
        index_t* labels;         /**< Store the labels associating spins here*/
        index_t* helper1;        /**< Same size as labels,
                                   for intermediate values during computation*/
        index_t* helper2;        /**< Same size as helper1*/
        index_t* lookuplabels;   /**< Store cluster sized at label indexes and
                                   use as a lookup table for the root labels*/
        index_t(* neighbours)[D]; /**< Lookup table for neighbouring sites*/
} swendsen_wang;
#endif //USE_SWENDSEN_WANG


#ifdef USE_WOLFF
/** Representation of one lattice site as it is used in the Wolff Single Cluster
 * Algorithm.
 * As unlike in the Swendsen Wang case only parts of the lattice are examined in
 * every update, it makes more sense to group all the information together like
 * this instead of using arrays for every single value, as we can not expect
 * caching benefits from that here.
 * This actually cries for object orientation, maybe I really have to sit down
 * and learn C++ one day...
 */
typedef struct Wolff_Site {
        int8_t spin;             /**< The spin on this lattice site,
                                   can eihter be UP or DOWN */
        int8_t flag;             /**< FLAGGED if already in cluster or pocket,
                                   else NOT_FLAGGED */
        sidelength_t L0_index;   /**< Index in the first direction, that is the
                                   direction in which the correlation length
                                   will be sampled later.*/
        index_t neighbours[2*D]; /**< Location of every neighour site */
} wolff_site;

/** A set of all lattice sites and some arrays neccessary for performing cluster
 * updates in the Wolff algorithm.
 * Contains all the information about all sites and two arrays, used in a list-
 * like manner to store the elements in the pocket and the cluster.
 * Also has a lookup table for the sum of cluster elements in the same time
 * slice, for a convenient access after the update.
 */
typedef struct Wolff {
        index_t N;            /**< The number of all lattice sites */
        sidelength_t L0;           /**< The length in first dimension, wich is
                                the direction in which the correlation function
                                will be sampled*/
        index_t cluster_size; /**< The number of elements already in the
                                cluster*/
        wolff_site* sites;    /**< All lattice sites in the model */
        wolff_site** cluster; /**< The actual sites in the cluster*/
        index_t* pocket;      /**< The sites, that grow the cluster*/
        index_t* time_slice;  /**< Sum of cluster sites in each time slice */
} wolff;
#endif //USE_WOLFF


/** Initializes the random number generator.
 * The random number generator, which has been selected at compile time, is
 * initialized with random seed from the seedfile string (first 4*624 bytes).
 * Seedfile should point to os-specific source (e.g. /dev/random/ to provide
 * for batch jobs) or via fixed seed from file for debugging. After that, it is
 * then warmed so that even worst case mersenne twister scenario (zero-excess in
 * the  seed) is overcome. Then the memory for the data structures is allocated
 * and initialized.
 * Note that the implementation makes use of global variables in every case,
 * so no thread safety is provided! If parallelization is desired, then this
 * has to be reconsidered...
 * @param seedfile path to initialize the random number generator from. Must be
 * valid!
 */
void initialize_rand(const char* seedfile);


#ifdef USE_WOLFF
/** Create a representation of a D-dimensional Ising Model with sidelengths L[i]
 * for the Wolff single cluster algorithm.
 * The necessary memory is allocated, initialized with the proper values (cold
 * start - all spins aligned initially) and returned as a pointer.
 * Will panic and exit the program in case no memory can be allocated!
 * @param L[] an array with D (corresponding to the dimension supplied at
 *      compile time) elements that has the side lengths of the lattice.
 *      Must not be NULL and each element must be greater than zero!
 * @return A pointer to the generated structure that can be used in subsequent
 *      calls to the other functions starting with wolff_* below.
 */
wolff* init_wolff(const sidelength_t L[D]);

/** Free the memory that has been allocated by a previous call to init_wolff().
 * @param w the data structure to be freed. Must not be NULL!
 */
void free_wolff(wolff* w);

/** Perform one cluster update.
 * That means build a cluster of equal spins with bond probability p from a
 * randomly chosen site and flip it. Also update the information in the wolff
 * data structure w for an efficient evaluation of the time slice correlation
 * function in a subsequent call to wolff_measure_time_slice_correlation()
 * before returning the size of the cluster that has just been created.
 * @param w the data structure, which has been returned by a previous call to
 *      init_wolff(), so that it contains the lattice and all other relevant
 *      information. Must not be NULL!
 * @param p The probability for a bond between adjacent spins of same sign.
 *      Must be positive and can in general be computed as 1-exp(-2*beta).
 * @return The site of the cluster that has just been created and flipped in the
 *      call.
 */
index_t wolff_perform_cluster_update( wolff* w, double p);

/** Measure the energy per lattice site of a given spin configuration produced
 * by the Wolff single cluster algorithm in a previous call of
 * wolff_perform_cluster_update().
 * The energy from each spin and its neighbour spin is computed in natural
 * units of J = -1, so an addintional coupling constant has to be multiplied
 * with the return value if desired.
 * @param w The data structure, which has been returned by a previous call to
 *      init_wolff(), so that it contains the lattice and all other relevant
 *      information. Must not be NULL!
 * @return The energy per lattice site in natural units. 
 */
double wolff_measure_energy(const wolff* w);

/** Measure the time slice correlation function in the direction of the first
 * dimension, i.e. the L[0] wich has been passed to init_wolff(L[D]) initially.
 * For the exact definition of what the time slice correlation function is,
 * refer to the accompanying bachelor thesis on Dimensional Crossover in the 
 * Ising Model (in German), for which this software has been initially written.
 * An improved estimator is used in the cluster picture as derived in
 * U. Wolff: ASYMPTOTIC FREEDOM AND MASS GENERATION IN THE 0(3) NONLINEAR 
 * sigma-MODEL
 * which greatly reduces the variance compared to the normal evaluation in the
 * spin picture. To further reduce the variance, the return value is averaged
 * over all time slices at distance t in the cluster 
 * @param w The data structure, which has been returned by a previous call to
 *      init_wolff(), so that it contains the lattice and all other relevant
 *      information. Must not be NULL!
 * @param t The distance between time slices. 0 <= t <= L[0]/2 must be true!
 *
 */
double wolff_measure_time_slice_correlation(const wolff* w, int t);
#endif //USE_WOLFF


#ifdef USE_SWENDSEN_WANG
/** Create a representation of a D-dimensional Ising Model with sidelengths L[i]
 * for the Swendsen Wang multi cluster algorithm.
 * The necessary memory is allocated, initialized with the proper values (cold
 * start - all spins aligned initially) and returned as a pointer.
 * Will panic and exit the program in case no memory can be allocated!
 * @param L[] an array with D (corresponding to the dimension supplied at
 *      compile time) elements that has the side lengths of the lattice.
 *      Must not be NULL and each element must be greater than zero!
 * @return A pointer to the generated structure that can be used in subsequent
 *      calls to the other functions starting with swendsen_wang_* below.
 */
swendsen_wang* init_swendsen_wang(const sidelength_t L[D]);

/** Free the memory that has been allocated by a previous call to 
 * init_swendsen_wang().
 * @param sw the data structure to be freed. Must not be NULL!
 */
void free_swendsen_wang(swendsen_wang* sw);

/** Perform an update of the whole lattice.
 * That means divide the whole lattice into multiple clusters of equal spins
 * with bond probability p and flip each cluster with probability 0.5.
 * @param sw the data structure, which has been returned by a previous call to
 *      init_swendsen_wang(), so that it contains the lattice and all other
 *      relevant information. Must not be NULL!
 * @param p The probability for a bond between adjacent spins of same sign.
 *      Must be positive and can in general be computed as 1-exp(-2*beta).
 */
void swendsen_wang_perform_sweep(const swendsen_wang* sw, double p);

/** Measure the energy per lattice site of a given spin configuration produced
 * by the Swendsen-Wang multi-cluster algorithm.
 * The energy from each spin and its neighbour spin is computed in natural
 * units of J = -1, so an addintional coupling constant has to be multiplied
 * with the return value if desired.
 * @param sw The data structure, which has been returned by a previous call to
 *      init_swendsen_wang(), so that it contains the lattice and all other
 *      relevant information. Must not be NULL!
 * @return The energy per lattice site in natural units. 
 */
double swendsen_wang_measure_energy(const swendsen_wang* sw);

/** Measure the time slice correlation function in the direction of the first
 * dimension, i.e. the L[0] wich has been passed to init_swendsen_wang(L[D])
 * initially.
 * For the exact definition of what the time slice correlation function is,
 * refer to the accompanying bachelor thesis on Dimensional Crossover in the 
 * Ising Model (in German), for which this software has been initially written.
 * An improved estimator is used in the cluster picture as derived in
 * U. Wolff: ASYMPTOTIC FREEDOM AND MASS GENERATION IN THE 0(3) NONLINEAR 
 * sigma-MODEL
 * which greatly reduces the variance compared to the normal evaluation in the
 * spin picture. To further reduce the variance, the return value is averaged
 * over all time slices at distance t in the cluster 
 * @param sw The data structure, which has been returned by a previous call to
 *      init_swendsen_wang(), so that it contains the lattice and all other
 *      relevant information. Must not be NULL!
 * @param t The distance between time slices. 0 <= t <= L[0]/2 must be true!
 *
 */
double swendsen_wang_measure_time_slice_correlation(const swendsen_wang* sw, int t);
#endif //USE_SWENDSEN_WANG

#endif //MONTE_CARLO_H
