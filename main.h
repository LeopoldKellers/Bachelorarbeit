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


/*! @file main.h
 * @brief Contains only the usage to be displayed when the -h option is passed
 * to the program and the data structure that contains all the relevant
 * information about the parameters for the simulation.
 */
#ifndef MAIN_H
#define MAIN_H


#define VERSION \
"%dD Ising Model Monte Carlo Simulation via Cluster Algorithm v1.0\n"\

#define USAGE VERSION \
"Usage: %s [options]\n  options:\n"\
"\n"\
"    -h, --help:          display this message\n"\
"\n"\
"    -s, --seedfile file: use file to customly seed the random number\n"\
"                         generator - default is %s\n"\
"                         point to file with random constants for debugging\n"\
"\n"\
"    -L,\n"\
"    --length l1, .. ,ln: Comma separated list of ls. Gives th sytem volume.\n"\
"                         Sets length of lattice side 1 to l1, side 2 to l2, "\
"etc.\n"\
"                         l1 is always the direction in which the time slice\n"\
"                         correlation function is sampled. The number of\n"\
"                         dimensions that the program has been compiled with\n"\
"                         specifies the maximum n. If n is less than that,\n"\
"                         then ln will be default for the ones unspecified!\n"\
"                         Each l must be greater than1, default is %d\n"\
"                         The product of all l, or, if only one l given l^%d\n"\
"                         must be less than(!) 2^31!\n"\
"\n"\
"    -J, --coupling j:    set coupling strength of neighbouring spins to j\n"\
"                         should not be 0.0, default is %f\n"\
"\n"\
"    -b, --beta bet:      bet -> inverse temperature the simulation runs at\n"\
"                         must be greater than zero, default is %f\n"\
"\n"\
"    -n, --nsamples #     collect # samples\n"\
"                         note: should be a power of 2 for best efficiency\n"\
"                         with binning, might be capped down otherwise\n"\
"                         must be at least 1, default is %d\n"\
"\n"\
"    -B, --binexp b:      put data in b bins of size up to 2^b to estimate \n"\
"                         the autocorrelation between adjacent samples\n"\
"                         must be nonnegative, default is 0, so no binning\n"\
"\n"\
"    -e, --equilibrate  # drop # sweeps before sampling to eqilibrate system\n"\
"                         must be at least 1, default is %d\n"\
"\n"\
"    -c, --calibrate      make a calibration run instead of sampling the \n"\
"                         time slice correlation function\n"\
"Please report bugs to l.k@wwu.de\n"
#define OUTPUT_USAGE(progname) printf(USAGE, D, progname, DEFAULT_SEEDFILE, DEFAULT_L, \
                              D, DEFAULT_J, DEFAULT_BETA, \
			      DEFAULT_NSWEEPS, DEFAULT_EQUILIBRATIONTIME)


/** A representation of the program state.
 * Contains all the relevant information about paramters for one run.
 */
typedef struct Configuration {
        char* seedfile;         /**< The filename to read the seed from*/
        int64_t nsamples;       /**< Number of samples to analyze in full run*/
        int16_t binexp;         /**< How many bins for error calculation max?*/
        int64_t equilibrationtime; /**< THe number of sweeps to drop initially*/
        double J;               /**< The coupling strength of neighbour spins*/
        double beta;            /**< The inverse temperature of the system*/
        sidelength_t L[D];      /**< The length of all lattice directions*/
        sidelength_t tmax;      /**< maximum correlation distance, aka L[0]/2 */
        index_t N;              /**< The number of lattice sites/spins*/
        int64_t lsamples;       /**< The number of elements to store temporarily
                                  before passing it to the data analysis module
                                  when sampling the correlation length. */
        bool calibrate;         /**< indicate if we are running a calibration */
} configuration;

#endif //MAIN_H
