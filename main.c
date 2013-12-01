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


/*! @file main.c
 * @brief Contains the main program, which parses the command line arguments, 
 * initializes the simulation, runs the simulation and data analysis to finally
 * spit out the results.
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include "common.h"
#include "monte_carlo.h"
#include "data_analysis.h"
#include "main.h"


/** Parse the command line arguments.
 * Returns the configuration of the program state with explicitly set or
 * default values.
 * First, the long and short forms of the allowed command line options are
 * evaluated, everything that remains is initialized by appropriate default
 * value.
 *
 * @param argc argc from main
 * @param argv argv from main
 * @return struct configuration of data relevant for program run
 */
static configuration* parse_options(int argc, char** argv) {
        configuration *c;
        c = xalloc(sizeof(configuration));
        memset(c, 0, sizeof(configuration)); //Initialize all values with zero
        //parse command line options - following the code in `man 3 getopt_long`
        int ch = -1;
        while (true) {
                int option_index = 0;
		double dval;
		long ival;
		char* endptr = NULL; //needed for strto*() checking numbers
                static struct option long_options[] = {
                        {"seedfile",    required_argument,      NULL, 's'},
                        {"length",      required_argument,      NULL, 'L'},
                        {"coupling",    required_argument,      NULL, 'J'},
                        {"beta",        required_argument,      NULL, 'b'},
                        {"nsamples",    required_argument,      NULL, 'n'},
                        {"equilibrate", required_argument,      NULL, 'e'},
                        {"binexp",      required_argument,      NULL, 'B'},
                        {"help",        no_argument,            NULL, 'h'},
                        {"calibrate",   no_argument,            NULL, 'c'},
                        {NULL,          0,                      NULL, 0}
                };
                ch = getopt_long (argc, argv, "hcs:L:J:b:B:n:e:",
                                long_options, &option_index);
                if (-1 == ch) {
                        break; // No more options to parse
                }
		errno = 0; // for strtol() and strtod()
                switch (ch) {
                        case 0:
                                fprintf(stderr, "Not supposed to happen! "
                                                "Options badly implemented!\n");
                                exit(EXIT_FAILURE);
                        case 'h':
                                OUTPUT_USAGE(argv[0]);
                                exit(EXIT_SUCCESS);
                        case 's':
                                c->seedfile = optarg;
                                break;
                        case 'c':
                                c->calibrate = true;
                                break;
			case 'L':
                        { 
                                char* l = NULL;
                                char* delim = ",";
                                char* rejectstr = "Rejected %s as side length!"
                                        "Invalid! number!\n";
                                int i = 0;
                                //parse the l1,l2, ... , ln comma separated
                                //argument list into its tokens
                                l = strtok(optarg, delim);
				ival = strtol(l, &endptr, 10);
				if (0 != errno || '\0' != *endptr
                                               || 1 << 15 <= ival
                                               || 1 > ival) {
					fprintf(stderr, rejectstr,
						l);
				}
				else {
                                        c->L[i] = ival;
				}
                                while (NULL != (l = strtok(NULL, delim))) {
                                        ++i;
                                        if (D <= i) {
                                                //Ok, fuck, I really get screwed
                                                //by indenting here. But it is
                                                //too logically coherent for
                                                //splitting it up into another
                                                //function. Then I will make
                                                //mistakes when changing code...
                                                fprintf(stderr, "Can't process"
                                                                " side length "
                                                                "%s - too many"
                                                                " in list (max."
                                                                " %d!)!!!\n",
                                                                l, D);
                                        }
                                        else {
                                                ival = strtol(l, &endptr, 10);
                                                if (0 != errno
                                                        || '\0' != *endptr
                                                        || 1 << 15 <= ival
                                                        || 1 > ival) {
                                                                fprintf(stderr,
                                                                 rejectstr, l); }
                                                else {
                                                        c->L[i] = ival;
                                                }
                                        }

                                }
                                while (i < D) {
                                        c->L[i] = ival;
                                        ++i;
                                }
                        }
				break;
			case 'J':
				dval = strtod(optarg, &endptr);
				if (0 != errno || '\0' != *endptr) {
					fprintf(stderr,"Rejected %s as coupling"
						" constant! Invalid number!\n",
						optarg);
				}
				else {
					c->J = dval;
				}
				break;
			case 'b':
				dval = strtod(optarg, &endptr);
				if (0 != errno || '\0' != *endptr) {
					fprintf(stderr, "Rejected %s as inverse"
						" temp.! Invalid number!\n",
						optarg);
				}
				else {
					c->beta = dval;
				}
				break;
			case 'n':
				ival = strtol(optarg, &endptr, 10);
				if (0 != errno || '\0' != *endptr) {
					fprintf(stderr, "Rejected %s as #"
						"sweeps! Invalid number!\n",
						optarg);
				}
				else {
					c->nsamples = ival;
				}
				break;
			case 'B':
				ival = strtol(optarg, &endptr, 10);
				if (0 != errno || '\0' != *endptr) {
					fprintf(stderr, "Rejected %s as 2^n"
						"bins! Invalid number!\n",
						optarg);
				}
				else {
					c->binexp = ival;
				}
				break;
			case 'e':
				ival = strtol(optarg, &endptr, 10);
				if (0 != errno || '\0' != *endptr) {
					fprintf(stderr, "Rejected %s as equilib"
						"rationtime! Invalid number!\n",
						optarg);
				}
				else {
					c->equilibrationtime = ival;
				}
				break;
                        case '?':
                                /* getopt_long already printed  error message.*/
                                break;
                        default:
                                debug_msgf("%s", "Not supposed to happen! "
                                                "Options badly implemented!\n");
                                exit(EXIT_FAILURE);
                }
        }
        //now check if values set or initialize with default values
        if (true == c->calibrate) {
                debug_msgf("%s\n", "Making a calibration run instead of "
                                "corrlation length samling!");
        }
        if (NULL == c->seedfile) {
                debug_msgf("No seed filename passed - falling back to default:"
                                " %s\n", DEFAULT_SEEDFILE);
                c->seedfile = DEFAULT_SEEDFILE;
        }
        else {
                debug_msgf("Using seedfile: %s\n", c->seedfile);
        }

        if (0 >= c->nsamples) {
                debug_msgf("No number of monte carlo sweeps passed - falling "
                                "back to default: %d\n", DEFAULT_NSWEEPS);
                c->nsamples = DEFAULT_NSWEEPS;
        }
        else{
                debug_msgf("Using %ld monte carlo sweeps\n", c->nsamples);
        }

        if (0 >= c->equilibrationtime) {
                debug_msgf("No number for equilibration sweeps passed - "
                                "falling back to default: %d\n",
                                DEFAULT_EQUILIBRATIONTIME);
                c->equilibrationtime= DEFAULT_EQUILIBRATIONTIME;
        }
        else{
                debug_msgf("Using %ld sweeps to equilibrate the system\n", 
                                c->equilibrationtime);
        }

        c->N = 1;
        for (int i = 0; i < D; ++i) {
                if (0 >= c->L[i]) {
                        debug_msgf("No side length in x_%d-direction passed - "
                                        "falling back to default: %d\n",
                                        i + 1, DEFAULT_L);
                        c->L[i] = DEFAULT_L;
                }
                else{
                        debug_msgf("Using %d as sidelengh in x_%d-direction\n",
                                        c->L[i], i+1);
                }
                c->N *= c->L[i];
        }
        //check that N is positive  and less than(! because rootlabels in
        //swendsen wang algo needs one more element!) the biggest possible
        //number of its datatype
        if(c->N + 1 <= 0) {
                fprintf(stderr, "Severe error: Too many lattice sites: %d",
                                c->N);
                exit(EXIT_FAILURE);
        }

        c->tmax = (sidelength_t) c->L[0] - 1;

        if (0.0 == c->J) {//usually floating point equality comparison is
                //unreliable, but here it is ok, as we just check intialization
                debug_msgf("No coupling strength const. passed - falling back "
                                "to default: %f\n", DEFAULT_J);
                c->J = DEFAULT_J;
        }
        else{
                debug_msgf("Using %f as coupling constant\n", c->J);
        }

        if (0 >= c->beta) {
                debug_msgf("No inverse temperature passed - falling back "
                                "to default: %f\n", DEFAULT_BETA);
                c->beta = DEFAULT_BETA;
        }
        else{
                debug_msgf("Using %f as inverse temperature\n", c->beta);
        }

        //check if nsweep is a power of two, otherwise cap to the next lower one
        int16_t i = 0;
        int64_t nsamples = c->nsamples;
        do {
                nsamples = nsamples >> 1;
                ++i;
        } while (nsamples != 1);
        if ((c->nsamples & (c->nsamples - 1)) != 0) {
                c->nsamples = nsamples << i;
                debug_msgf("number of sweeps is not a power of 2! capping to "
                                "%ld!\n", c->nsamples);
        }

        if (0 == c->binexp) {
                debug_msgf("No bin exponent passed or binning switched off - "
                                "- falling back to default: %d\n", 0);
        }
        else if (c->binexp > (signed) sizeof(index_t) * 8 - 2 || //fit data type
                        (1 << ((index_t) c->binexp + 1)) > c->nsamples) {
                c->binexp = (int16_t) (i - 1);
                debug_msgf("Number of samples must be at least 2^max. "
                                "bin size - falling back to binexp %d\n",
                                c->binexp);
        }
        else{
                debug_msgf("Using %d bins of maximum size %ld \n",
                                c->binexp, ((int64_t) 1) << c->binexp);
        }

        // + 1 because we also want to be able to sample t = 0
        c->lsamples = ((int64_t) 1 << c->binexp) * (c->tmax + 1);

        //and finally give back the program configuration to work with
        return c;
}


/** Runs the main program (who would have thought it?).
 *
 * An overview of the available parameters that can be passed to the main
 * program can be printed by running it with the -h or --help option.
 * Has the command line arguments parsed and then equlibrates the system with
 * the user specified or default value. After that, the correlation length is
 * sampled with the algorithm(s) specified at compile time.
 *
 * There is a calibration mode, in which first the evolution of the energy per
 * site for each sweep during the  equilibration phase is printed, enabling a
 * visual inspection if the system has settled in thermodynamic equilibrium.
 * 
 * @param argc The number of arguments as usual
 * @param argv The options passed from the shell as usual
 * @return EXIT_SUCESS if everything went fine, EXIT_FAILURE if not enough 
 *      memory could be allocated. Also, in debug mode, an assertion may fail.
 *      
 */
int main(int argc, char** argv) {
        debug_msgf(VERSION, D);
        //some stuff to print time information first, before we actually start
        //to get busy
        clock_t tic = 0;
        clock_t toc = 0;
        double elapsed_time = 0;
        struct timeval tv;
        struct tm* ptm;
        char t_str[50];
        long ms;

        gettimeofday (&tv, NULL);
        ptm = localtime (&tv.tv_sec);
        strftime (t_str, sizeof (t_str), "%d.%m.%Y %H:%M:%S", ptm);
        ms = tv.tv_usec / 1000;
        debug_msgf("Run at: %s.%ld\n", t_str, ms);

        const configuration* const restrict c = parse_options(argc, argv);
        /* Ok, so admittedly all theses ifdefs make the code a bit unreadable,
         * although the actual code structure is not that difficult. It has the
         * advantage over if(...) {...} though, that the compiler will tell you,
         * if you use one of the structures when you are not allowed to.
         * If you prefer, you can just throw them all out, then you will just
         * end up with both algorithms in parallel, which is not a problem at
         * all . In that case make sure to define both USE_WOLFF and
         * USE_SWENDSEN_WANG either in monte_carlo.h or via compiler switch! */
        initialize_rand(c->seedfile);
#ifdef USE_WOLFF
        wolff* const restrict w = init_wolff(c->L);
#endif
#ifdef USE_SWENDSEN_WANG
        swendsen_wang* const restrict sw = init_swendsen_wang(c->L);
#endif

        const double p = 1 - exp(-2.0 * c->beta * c->J);

        size_t size = c->lsamples * sizeof(double);

#ifdef USE_WOLFF
        double* const restrict w_samples = xalloc(size);
#endif
#ifdef USE_SWENDSEN_WANG
        double* const restrict sw_samples = xalloc(size);
#endif

        // initialize t + 1 statistics units because we also want to sample t=0
        size = (c->tmax + 1) * sizeof(binned_statistics*);
        //same as with *_samples above
#ifdef USE_WOLFF
        binned_statistics** const restrict w_stat = xalloc(size);
        binned_statistics* w_stat_energy = NULL;
#endif
#ifdef USE_SWENDSEN_WANG
        binned_statistics** const restrict sw_stat = xalloc(size);
        binned_statistics* sw_stat_energy = NULL;
#endif
        for (int i = 0; i <= c->tmax; ++i) {
#ifdef USE_WOLFF
                w_stat[i] = init_statistics(c->binexp);
#endif
#ifdef USE_SWENDSEN_WANG
                sw_stat[i] = init_statistics(c->binexp);
#endif
        }

        if (c->calibrate) {
                toc = clock();
                elapsed_time =  (double)(toc - tic) / CLOCKS_PER_SEC;
                fprintf(stderr, "#Dropped %ld samples in %f seconds!\n",
                                c->equilibrationtime, elapsed_time);
                fprintf(stderr, "#sweep\t");
#ifdef USE_WOLFF
                w_stat_energy = init_statistics(0);
                fprintf(stderr, "e_W\t\t<e_W>\t\t");
#endif
#ifdef USE_SWENDSEN_WANG
                sw_stat_energy = init_statistics(0);
                fprintf(stderr, "e_SW\t\t<e_SW>\t\t");
#endif
                fprintf(stderr, "\n");
                tic = clock(); //do this for measuring execution time
        }
        int64_t j = 0;
        for (int64_t i = 1; i < c->equilibrationtime; ++i) {
#ifdef USE_WOLFF
                //remember that one cluster update is not a complete sweep here!
                while (j < c->N) {
                        j += wolff_perform_cluster_update(w, p);
                }
                j -= c->N;
#endif
#ifdef USE_SWENDSEN_WANG
                swendsen_wang_perform_sweep(sw, p);
#endif
                if (c->calibrate) {// print the energy if we calibrate
                        double e;
#ifdef USE_WOLFF
                        e = wolff_measure_energy(w);
                        update_statistics(w_stat_energy, &e);
                        fprintf(stderr, "%ld\t%e\t%e",
                                        i, e, get_mean(w_stat_energy));
#endif
#ifdef USE_SWENDSEN_WANG
                        e = swendsen_wang_measure_energy(sw);
                        update_statistics(sw_stat_energy, &e);
                        fprintf(stderr, "%ld\t%e\t%e",
                                        i, e, get_mean(sw_stat_energy));
#endif
                        fprintf(stderr, "\n");
                }

        }
        if (c->calibrate) {
#ifdef USE_WOLFF
                free(w_stat_energy);
#endif
#ifdef USE_SWENDSEN_WANG
                free(sw_stat_energy); 
#endif

        }
	debug_msgf("%s", "Finished initial skirmish - start sampling and analyzing (will take some time)\n");

        assert(0 == (c->nsamples & (c->nsamples - 1))); //implementation assumes
        //that the overall number of samples to take is a power of two here
        int64_t samples_per_bin = (1 << c->binexp);
        for (int64_t i = 0; i < (c->nsamples >> c->binexp); ++i) {
                for (int64_t j = 0; j < samples_per_bin; ++j) {
                        for (int t = 0; t <= c->tmax; ++t) {
                                //s should rather be called t'
                                //s does not go thorough t consecutively,
                                //but rather in the order of
                                //0,2,4, ..., tmax, tmax-1 , tmax-3, ... 1
                                int s = 2 * t;
                                if (s > c->tmax) {
                                        s = 2 * c->tmax - s + 1;
                                }
#ifdef USE_WOLFF
                                wolff_perform_cluster_update(w, p); 
                                w_samples[s * samples_per_bin + j] =                                                                                                                          
                                        wolff_measure_time_slice_correlation(w, s); 
#endif
#ifdef USE_SWENDSEN_WANG
                                swendsen_wang_perform_sweep(sw, p); 
                                sw_samples[s * samples_per_bin + j] =
                                        swendsen_wang_measure_time_slice_correlation(sw, s); 
#endif
                        }
                }
                for (int t = 0; t <= c->tmax; ++t) {
#ifdef USE_WOLFF
                        update_statistics(w_stat[t],
                                        w_samples + t * samples_per_bin);
#endif
#ifdef USE_SWENDSEN_WANG
                        update_statistics(sw_stat[t],
                                        sw_samples + t * samples_per_bin);
#endif
                }

        }
        debug_msgf("%s", "Finished sampling - printing data analysis!\n");
        tic = clock();
        elapsed_time =  (double)(tic - toc) / CLOCKS_PER_SEC;
        printf("#Evaluated %ld samples in %f seconds!\n",
                                c->nsamples, elapsed_time);
        printf("#t\tb\t");
#ifdef USE_WOLFF
        printf("G_W(t)\t\t2*std_dev_W\tstd_err_W\t");
#endif
#ifdef USE_SWENDSEN_WANG
        printf("G_SW(t)\t\t2*std_dev_SW\tstd_err_SW\t");
#endif
        printf("\n");
        for (int t = 0; t <= c->tmax; ++t) {
                for (int b = 0; b <= c->binexp; ++b) {
                        printf("%d\t%d\t", t, b);
#ifdef USE_WOLFF
                        printf("%e\t%e\t%e\t",
                                        get_mean(w_stat[t]),
                                        2.0 * sqrt(get_variance(w_stat[t], b)),
                                        get_standard_error(w_stat[t], b));
#endif
#ifdef USE_SWENDSEN_WANG
                        printf("%e\t%e\t%e",
                                        get_mean(sw_stat[t]),
                                        2.0 * sqrt(get_variance(sw_stat[t], b)),
                                        get_standard_error(sw_stat[t], b));
#endif
                        printf("\n");
                }
        }

        //not really neccessary as we exit anyways...
        //free(NULL) does not have any effect
        for (int t = 0; t < c->tmax - 1; ++t) {
#ifdef USE_WOLFF
                free_statistics(w_stat[t]);
#endif
#ifdef USE_SWENDSEN_WANG
                free_statistics(sw_stat[t]);
#endif
        }
#ifdef USE_WOLFF
        free(w);
        free(w_samples);
#endif
#ifdef USE_SWENDSEN_WANG
        free_swendsen_wang(sw);
        free(sw_samples);
#endif
        return EXIT_SUCCESS;
}
