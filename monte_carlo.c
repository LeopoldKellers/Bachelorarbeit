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


/*! @file monte_carlo.c
 * @brief Contains the implementation of the stuff described in monte_carlo.h.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include "common.h"
#include "monte_carlo.h"


/** Generate a random number in the interval [0;1).
 * Basically only a wrapper to interface the appropriate functions in the rngs.
 * Note that the implementation makes use of global variables in every case (not
 * my fault, though), so no thread safety is provided! If parallelization is
 * desired, then this has to be reconsidered...
 * @return the random number in the interval [0;1) 
 */
static inline double get_rand(){
#if RNG == 1 //dsfmt
        return dsfmt_gv_genrand_close_open();
        /* The code below did prove a wee bit more efficient than fetching
         * single values like above in this specific case. But the choice of
         * randblocksize appears to be quite delicate in order to maintain that
         * advantage, so the hassle is just not worth it, as the random number
         * generation is just a very tiny fraction of the actual computation
         * time in most cases.*/
        /*
        static double rands[RANDBLOCKSIZE];
        static int i = 0;
        if (RANDBLOCKSIZE == i) { // create a new bunch of random numbers
                i = 0;
                dsfmt_gv_fill_array_close_open(rands, RANDBLOCKSIZE);
        }
        return rands[i++];
        */
#elif RNG == 2 //well
        //accumulation of rands as above appears to be slightly more inefficient
        return (*WELLRNG19937a)();
#else 
#error "This is not supposed to be possible - definition in header file screwed\n"
#endif
}

void initialize_rand(const char* seedfile){
        //seed the random number generator from file
        uint32_t seed_array[SEEDLENGTH];
        FILE* f;
        assert(NULL != seedfile);
        f = fopen(seedfile, "r");
        if (NULL == f){
                fprintf(stderr, "Error: could not open random number generator"
                                " seed file %s - exiting!\n", seedfile);
                exit(-1); 
        }
        if (fread(&seed_array, sizeof(int), SEEDLENGTH, f) != SEEDLENGTH){
                fprintf(stderr, "Error: could not read %lu bytes from %s - "
                                "exiting!\n", sizeof(seed_array), seedfile);
                exit(-1);
        }
        if (0 != fclose(f)){
                fprintf(stderr, "Could not close %s, but we have the data we "
                                "wanted anyways...\n", seedfile);
        }
#if RNG == 1 //dsfmt
        dsfmt_gv_init_by_array(seed_array, SEEDLENGTH);
        debug_msgf("Using random number generator %s\n", dsfmt_get_idstring());
#elif RNG == 2 //well
        InitWELLRNG19937a(seed_array);
        debug_msgf("Using random number generator %s\n", "WELLRNG19937c");
        //19937c because of tempering option in Makefile
#else
#error "This is not supposed to be possible - definition in header file screwed\n"
#endif
        //warm up the random number generator - Mersenne Twister algorithms can
        //suffer from bad (zero-excess) initialization which in this case is
        //unprobabl but theoretically possible. Discarding a million numbers
        //should be sufficient and only takes a fraction of a second.
        int32_t Nrands = 1E6;
        for (int32_t i = 0; i < Nrands; ++i) {
                get_rand();
        }
}

#ifdef USE_WOLFF
wolff* init_wolff(const sidelength_t L[D]) {
	assert(NULL != L);

        wolff* w = xalloc(sizeof(wolff));

        w->L0 = L[0];
        //how many sites in total?
        w->N = 1;
        for (int i = 0; i < D; ++i) {
                assert(0 != L[i]);
                w->N *= L[i];
        }


        //allocate the memory for the arrays
        int64_t size = 0;
        
        size = w->N * sizeof(*w->sites);
        w->sites = xalloc(size);

        size = w->N * sizeof(*w->cluster);
        w->cluster = xalloc(size);

        size = w->N * sizeof(*w->pocket);
        w->pocket = xalloc(size);

        size = w->L0 * sizeof(*w->time_slice);
        w->time_slice = xalloc(size);


        index_t accumulated_L[D + 1]; //helper construct to save calculations
        for (int i = 0; i <= D; ++i) {
                accumulated_L[i] = 1;
                for(int d = 0; d < i; ++d) {
                        accumulated_L[i] *= L[d];
                }
        }
        //initialize the lattice sites
        index_t index1;
        index_t index2;
        for (index_t i = 0; i < w->N; ++i) {
                //cold start with all spins up initially
                wolff_site* s = w->sites + i;
                s->spin = UP;
                s->flag = NOT_FLAGGED;
                s->L0_index = (sidelength_t) (i % w->L0);
                //fill the lookup table with 2 * D neighbours for each spin
                //in all directions and periodic boundary conditions
                for (int d = 0; d < D; ++d) {
                        //fill neighbour tables in a way similar to that in the
                        //Swendsen Wang implementation below
                        index1 = i - (i % accumulated_L[d + 1]);
                        index1 += ((int64_t) i + accumulated_L[d])
                                % accumulated_L[d + 1];
                        index2 = i - (i % accumulated_L[d + 1]);
                        index2 += ((int64_t) i - accumulated_L[d])
                                % accumulated_L[d + 1];
                        //Note: sign dependence of modulo operation only
                        //standardized in c99 so make sure to compile with that!
                        if (index2 < 0) {
                                index2 += accumulated_L[d + 1];
                        }
                        s->neighbours[2 * d] = index1;
                        s->neighbours[2 * d + 1] = index2;
                }

        }
        //uncomment below to print the neighbour table and check the result
        //after making any changes if test fails!
        /*
        for(int64_t i = 0; i < w->N; ++i) {
                printf("%ld:\t", i);
                for (int d = 0; d < 2*D; ++d) {
                        printf("%d\t", (w->sites + i)->neighbours[d]);
                }
                        printf("%ld\n", (w->sites) + i)->L0_index);
        }
        */
        return w;
}

void free_wolff(wolff* w) {
        assert(NULL != w);
        free(w->sites);
        free(w->cluster);
        free(w->pocket);
        free(w->time_slice);
        free(w);
}


index_t wolff_perform_cluster_update(wolff* const w, const double p) {
        assert(NULL != w);
        assert(NULL != w->sites);
        assert(NULL != w->cluster);
        assert(NULL != w->pocket);
        assert(NULL != w->time_slice);
        wolff_site* const sites = w->sites;
        wolff_site** const cluster = (*w).cluster;
        index_t* const restrict pocket = w->pocket;
        index_t* const restrict time_slice = w->time_slice;
        const index_t N = w->N;
        const sidelength_t L0 = w->L0;
        memset(time_slice, 0, L0 * sizeof(*w->time_slice));
        index_t cluster_first_free = 0; //points to the first non-occupied
                                        //element of the array-list cluster
        index_t pocket_first_free = 1;  //points to the first non-occupied
                                        //element of the array-list pocket
        index_t pocket_first_elem = 0;  //the head of the array-list pocket

        //start by putting a randomly chosen site into the pocket to grow
        //the whole cluster from lateron
        index_t site_index = (index_t) (get_rand() * N);
        assert(NOT_FLAGGED == sites[site_index].flag);
        sites[site_index].flag = FLAGGED;
        pocket[0] = site_index;
        //now take the first element from the pocket and move it into the
        //cluster after adding all neighbours with same spin that have not yet
        //been flagged to the cluster with probability p until all pocket
        //elements have been processed.
        while (pocket_first_elem < pocket_first_free) {
                //some assertions left over from fixing memory corruption 
                assert(pocket_first_free <= N);
                assert(cluster_first_free <= N);
                site_index = pocket[pocket_first_elem++];
                wolff_site* const site = sites + site_index;
                cluster[cluster_first_free++] = site;
                assert(NULL != site->neighbours);
                assert(site->spin == UP || site->spin == DOWN);
                assert(site->L0_index < L0);
                time_slice[site->L0_index] += 1;
                for (int i = 0; i < 2 * D; ++i) {
                        const index_t neigh_index = site->neighbours[i];
                        wolff_site* const neigh  = sites + neigh_index;
                        if (site->spin == neigh->spin
                                && NOT_FLAGGED == neigh->flag
                                && get_rand() < p) {
                                neigh->flag = FLAGGED;
                                pocket[pocket_first_free++] = neigh_index;
                        }
                }
        }
        //the index of the first free element in the cluster coincides with
        //the cluster size
        w->cluster_size = cluster_first_free;
        //now flip the spins in the cluster
        const int8_t spin = cluster[0]->spin == UP ? DOWN : UP;
        while (0 < cluster_first_free) {
                //just make sure all elements of the cluster actually have 
                //the same spin...
                assert(cluster[cluster_first_free - 1]->spin != spin);
                //and on top of flipping also unflag it
                cluster[--cluster_first_free]->spin = spin;
                assert(FLAGGED == cluster[cluster_first_free]->flag);
                cluster[cluster_first_free]->flag = NOT_FLAGGED;
        }
        return w->cluster_size;
}

double wolff_measure_time_slice_correlation(const wolff* const w, int t) {
        assert(NULL != w);
        assert(NULL != w->time_slice);
        assert(0 < w->cluster_size && 0 < w->L0);
        assert(t < w->L0);
        const index_t* const restrict time_slice = w->time_slice;
        //improved estimator, averaging the product of the number of sites in
        //the time slices at distance t from each other over the whole lattice
        int64_t sum = 0;
        for (int i = 0; i < w->L0; ++i) { 
                sum += time_slice[i] * time_slice[t];
                // increment of t modulo L0, so that we always start with the
                // correct first site
                if (++t == w->L0) {
                        t = 0;
                }
        }
        return sum / (double) w->L0 / w->cluster_size;// * w->N / w->N = 1
}

double wolff_measure_energy(const wolff* const w) {
        assert(NULL != w);
        assert(NULL != w->sites);
        const wolff_site* const sites = w->sites;
        //the primitive method of summing over all nearest neighbours...
        index_t sum = 0;
        for (index_t i = 0; i < w->N; ++i) {
                //+= 2, because we want to count each neighbour only once!
                for (int d = 0; d < 2 * D; d += 2) {
                        const wolff_site* const site = sites + i;
                        const wolff_site* const neighbour = sites + site->neighbours[d];
                        sum += (site->spin == neighbour->spin) ? 1 : -1;
                }
        }
        return -sum / (double) w->N;
}
#endif //USE_WOLFF

#ifdef USE_SWENDSEN_WANG
swendsen_wang* init_swendsen_wang(const sidelength_t L[D]) {
	assert(NULL != L);

        swendsen_wang* sw = xalloc(sizeof(swendsen_wang));

        sw->L0 = L[0];
        //how many sites do we have in total?
        sw->N = 1;
        for (int i = 0; i < D; ++i) {
                assert(0 != L[i]);
                sw->N *= L[i];
        }


        //allocate the memory for the arrays
        int64_t size = 0;
        
        size = sw->N * sizeof(*sw->spins);
        sw->spins = xalloc(size);

        size = sw->N * sizeof(*sw->helper1);
        sw->helper1 = xalloc(size);

        size = sw->N * sizeof(*sw->helper2);
        sw->helper2 = xalloc(size);

        size = sw->N * sizeof(*sw->labels);
        sw->labels = xalloc(size);

        size = (sw->N + 1) * sizeof(*sw->lookuplabels); //+1 because label 0 invalid
        sw->lookuplabels = xalloc(size);

        size = D * sw->N * sizeof(*sw->neighbours);
        sw->neighbours = xalloc(size);
        
        //initialize the spins all up (cold start)
        memset(sw->spins, UP, sw->N);

        //helper construct to save calculations with neighbours
        index_t accumulated_L[D + 1];
        for (int i = 0; i <= D; ++i) {
                accumulated_L[i] = 1;
                for(int d = 0; d < i; ++d) {
                        accumulated_L[i] *= L[d];
                }
        }
        //fill the lookup table with D neighbours for each spin
        //and periodic boundary conditions
        index_t index;
        for(index_t i = 0; i < sw->N; ++i) {
                for (int d = 0; d < D; ++d) {
                        index = i - (i % accumulated_L[d + 1]);
                        index += ((int64_t) i + accumulated_L[d])
                                % accumulated_L[d + 1];
                        //to make yourself clear about why this formula works,
                        //try it on paper with D=3 and L=3 for example ;-)
                        sw->neighbours[i][d] = index;
                }
        }
        //uncomment below to print the neighbour table for checking after making
        //any changes if test fails!
        /*
        for(int64_t i = 0; i < sw->N; ++i) {
                printf("%ld:\t", i);
                for (int d = 0; d < D; ++d) {
                        printf("%d\t", sw->neighbours[i][d]);
                }
                        puts("");
        }
        */
        return sw;
}

void free_swendsen_wang(swendsen_wang* sw) {
        assert(NULL != sw);
        free(sw->spins);
        free(sw->labels);
        free(sw->lookuplabels);
        free(sw->neighbours);
        free(sw->helper1);
        free(sw->helper2);
        free(sw);
}

/** Looks up the root label of the element at index t in the lookuplabels array.
 * @param lookuplabels An array, that contains the number of lattice sites with
 * the same label (> 0) for the root label or a (multi stage) reference to the
 * actual root label.
 * @param t The label for which the actual root label is desired.
 * @return The original label that t refers to actually
 */
inline static index_t lookup_label(const index_t* const lookuplabels, index_t l) {
        assert(0 < l);
        while (lookuplabels[l] <= 0) {
                l = lookuplabels[l];
                l = 0 < l ? l : -l;
        }
        return l;
}

void swendsen_wang_perform_sweep(const swendsen_wang* const sw, const double p) {
        assert(NULL != sw);
        assert(NULL != sw->spins);
        assert(NULL != sw->labels);
        assert(NULL != sw->lookuplabels);
        assert(NULL != sw->neighbours);
        assert(NULL != sw->helper1);
        int8_t* const restrict spins = sw->spins;
        index_t* const restrict labels = sw->labels;
        index_t* const restrict helper1 = sw->helper1;
        index_t* const restrict lookuplabels = sw->lookuplabels;
        index_t(* const restrict neighbours)[D] = (*sw).neighbours;
        const index_t N = sw->N;
        index_t maxlabel = 0;
        bool bond;
        memset(lookuplabels, 0, (N + 1) * sizeof(*lookuplabels));
        memset(labels, 0, N * sizeof(*labels));
        memset(helper1, 0, N * sizeof(*helper1));
        //Implementation of the Swendsen Wang cluster algorithm using the
        //Hoshen Kopelman cluster counting approach
        for (index_t i = N - 1; i >= 0; --i) {
//                nbonds = 0;
                for (int d = 0; d < D; ++d) {
                        const index_t n = neighbours[i][d]; //current neighbour

                        //Ok, what is going on below?
                        //We have the following possibilities (no label means
                        //that label[i] = 0):
                        //no bond and no label: assign new label to site
                        //no bond and already labeled: do nothing
                        //bond from no label to no label: assign both new label
                        //bond from no label to label:
                        //      assign neighbour label to site
                        //bond from label to no label: 
                        //      assign site label to neighbour
                        //bond from label to same label: do nothing
                        //bond from label to different label:
                        //      unite the labels via the lookuplabel structure
                        if (spins[i] == spins[n]) {
                                //assign bond with probability p
                                if (get_rand() < p) {
                                        bond = true;
                                }
                                else {
                                        bond = false;
                                }
                        }
                        else { //no bonds between opposing spins
                                bond = false;
                        }
                        //ok, indention will get really bad really soon, but it
                        //is also hard to move this part of the code into a 
                        //separate function as we uneccessarily have to pass
                        //and return the value for maxlabel all the time then.
                        //So you will need a widescreen for reading this part
                        //of the algorithm :-)
                        if (true == bond) { //make new bond
                                if (0 == labels[i]) { //no site label
                                        if (0 == labels[n]) {//no neighbour label
                                                // assign both new label
                                                ++maxlabel;
                                                labels[i] = maxlabel;
                                                labels[n] = maxlabel;
                                                //new cluster created
                                                assert(0 == lookuplabels[maxlabel]);
                                                lookuplabels[maxlabel] += 2;
                                        }
                                        else { //neighbour has label
                                                //assign neighbour label to site
                                                labels[i] = lookup_label(lookuplabels, labels[n]);
                                                // increase cluster size by one
                                                assert(0 < lookuplabels[labels[i]]);
                                                lookuplabels[labels[i]] += 1;
                                        }
                                }
                                else { // site has label
                                        if (0 == labels[n]) { //neighbour has not
                                                //assign site label
                                                //to neighbour
                                                labels[n] = lookup_label(lookuplabels, labels[i]);
                                                // increase cluster size by one
                                                assert(0 < lookuplabels[labels[n]]);
                                                lookuplabels[labels[n]] += 1;
                                        }
                                        else{ //both sites are labelled
                                              //the most common case  
                                                index_t rootlabeli = lookup_label(lookuplabels, labels[i]);
                                                index_t rootlabeln = lookup_label(lookuplabels, labels[n]);
                                                assert(0 < lookuplabels[rootlabeli]);
                                                assert(0 < lookuplabels[rootlabeln]);
                                                if (rootlabeli != rootlabeln) {
                                                        //bond between different
                                                        //labels: unite the
                                                        //labels via the
                                                        //rootlabel structure
                                                        index_t newlabel = (rootlabeli > rootlabeln) ? rootlabeln : rootlabeli;
                                                        index_t oldlabel = (rootlabeli > rootlabeln) ? rootlabeli : rootlabeln;
                                                        lookuplabels[newlabel] += lookuplabels[oldlabel];
                                                        lookuplabels[oldlabel] = -newlabel;
                                                }
                                                else {//same label on both sites
                                                        //do nothing
                                                        ;
                                                }
                                        }
                                }
                        }
                        else { // no new bond
                                if (0 == labels[i]) { // site has no own label
                                        //assign new label to site
                                        ++maxlabel;
                                        labels[i] = maxlabel;
                                        //new cluster created
                                        assert(0 == lookuplabels[maxlabel]);
                                        lookuplabels[maxlabel] += 1;
                                }
                                else{ //site already has own label
                                        //do nothing
                                        assert(0 < lookuplabels[lookup_label(lookuplabels, labels[i])]);
                                }
                        }

                }
        }

        //reduce all labels to their root label and flip each cluster with 
        //probability 0.5
        index_t total_clustersize = 0;
        for (index_t i = 0; i < N; ++i) {
                labels[i] = lookup_label(lookuplabels, labels[i]);
                assert(0 == spins[i] || UP == spins[i]
                                || DOWN == spins[i]);
                //associate each cluster with a spin if not yet happened
                if (0 == helper1[labels[i] - 1]) {
                        index_t spin = (get_rand() < 0.5) ? UP : DOWN; 
                        helper1[labels[i] - 1] = spin;
                }
                //and apply that spin
                spins[i] = (int8_t) helper1[labels[i] - 1];
                //some consistency check
                if (lookuplabels[i+1] > 0) {
                        total_clustersize += lookuplabels[i+1];
                }
        }
        assert(N == total_clustersize); //make sure, every spin labelled
}

double swendsen_wang_measure_time_slice_correlation(const swendsen_wang* const sw, int t) {
        assert(NULL != sw);
        assert(NULL != sw->helper1);
        assert(NULL != sw->helper2);
        assert(NULL != sw->labels);
        assert(NULL != sw->lookuplabels);
        assert(t < sw->L0);
        const index_t* const restrict labels = sw->labels;
        index_t* const restrict lookuplabels = sw->lookuplabels;
        index_t* const restrict helper1 = sw->helper1;
        index_t* const restrict helper2 = sw->helper2;
        const index_t N = sw->N;
        const sidelength_t L0 = sw->L0;
        double sum = 0;
        //take average over all time slices at distance t for better statistics
        for(index_t i = 0; i < L0; ++i) {
                memset(lookuplabels, 0, (N + 1) * sizeof(*lookuplabels));
                memset(helper1, 0, N * sizeof(*helper1));
                memset(helper2, 0, N * sizeof(*helper2));
                index_t k = t;
                //store the amount of lattice sites with a label in two arrays
                index_t maxlabel = 0;
                for (index_t j = i; j < N; j += L0) {
                        // - 1 because labels start counting at 1
                        index_t m = labels[j] - 1;
                        index_t n = labels[k] - 1;
                        if (0 == helper2[m] + (helper1[m])++){
                                lookuplabels[maxlabel++] = m;
                        }
                        if (0 == helper1[n] + (helper2[n])++){
                                lookuplabels[maxlabel++] = n;
                        }
                        k += L0;
                }
                //and sum up their products, for which, to make it less horribly
                //inefficient as we would add up zero most of the time track is
                //kept of the touched indexes and only these are evaluated,
                for (k = 0; k < maxlabel; ++k) {
                        index_t j = lookuplabels[k];
                        sum += helper1[j] * helper2[j];
                }
                //increment t and take modulo L0 to make sure we start
                //at the right position all the time
                if (++t == L0) {
                        t = 0;
                }
        }
        return sum / N / L0;
}

double swendsen_wang_measure_energy(const swendsen_wang* const sw) {
        assert(NULL != sw);
        assert(NULL != sw->spins);
        assert(NULL != sw->neighbours);
        int8_t* const restrict spins = sw->spins;
        index_t(* const restrict neighbours)[D] = (*sw).neighbours;
        //use the easiest method by simply summing over all neighbours
        index_t sum = 0;
        for (index_t i = 0; i < sw->N; ++i) {
                for (int d = 0; d < D; ++d) {
                        sum += (spins[i] == spins[neighbours[i][d]]) ? 1 : -1;
                }
        }
        return -sum / (double) sw->N;
}
#endif //USE_SWENDSEN_WANG

