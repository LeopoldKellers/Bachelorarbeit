#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <common.h>
#include <monte_carlo.h>
#include <data_analysis.h>

double eq(double a, double b, double eps) {
        double diff = a - b;
        diff = diff > 0 ? diff : -diff;
        return diff < eps;
}

void test_statistics() {
        double vals[16];
        binned_statistics* stat = init_statistics(4);

        for (int i = 0; i < 16; ++i) {
                vals[i] = i + 1;
        }
        update_statistics(stat, vals);
        //just check against the corresponding values computed manually
        assert(eq(8.5, get_mean(stat), 0.0001));
        assert(eq(22.6667, get_variance(stat, 0), 0.0001));
        assert(eq(24.0, get_variance(stat, 1), 0.0001));
        assert(eq(26.6667, get_variance(stat, 2), 0.0001));
        assert(eq(32.0, get_variance(stat, 3), 0.0001));
        assert(eq(0.0, get_variance(stat, 4), 0.0001));
        assert(eq(1.19024, get_standard_error(stat, 0), 0.0001));
        assert(eq(1.73205, get_standard_error(stat, 1), 0.0001));
        assert(eq(2.58199, get_standard_error(stat, 2), 0.0001));
        assert(eq(4.0, get_standard_error(stat, 3), 0.0001));
        assert(eq(0.0, get_standard_error(stat, 4), 0.0001));
        for (int i = 0; i < 16; ++i) {
                vals[i] = i + 17;
        }
        update_statistics(stat, vals);
        assert(eq(16.5, get_mean(stat), 0.0001));
        assert(eq(88.0, get_variance(stat, 0), 0.0001));
        assert(eq(90.6667, get_variance(stat, 1), 0.0001));
        assert(eq(96.0, get_variance(stat, 2), 0.0001));
        assert(eq(106.6667, get_variance(stat, 3), 0.0001));
        assert(eq(128.0, get_variance(stat, 4), 0.0001));
        assert(eq(1.65831, get_standard_error(stat, 0), 0.0001));
        assert(eq(2.38048, get_standard_error(stat, 1), 0.0001));
        assert(eq(3.4641, get_standard_error(stat, 2), 0.0001));
        assert(eq(5.16398, get_standard_error(stat, 3), 0.0001));
        assert(eq(8.0, get_standard_error(stat, 4), 0.0001));
        for (int i = 0; i < 16; ++i) {
                vals[i] = i + 33;
        }
        update_statistics(stat, vals);
        for (int i = 0; i < 16; ++i) {
                vals[i] = i + 49;
        }
        update_statistics(stat, vals);
        assert(eq(32.5, get_mean(stat), 0.0001));
        assert(eq(346.6667, get_variance(stat, 0), 0.0001));
        assert(eq(352.0, get_variance(stat, 1), 0.0001));
        assert(eq(362.6667, get_variance(stat, 2), 0.0001));
        assert(eq(384., get_variance(stat, 3), 0.0001));
        assert(eq(426.6667, get_variance(stat, 4), 0.0001));
        assert(eq(2.32737, get_standard_error(stat, 0), 0.0001));
        assert(eq(3.31662, get_standard_error(stat, 1), 0.0001));
        assert(eq(4.76095, get_standard_error(stat, 2), 0.0001));
        assert(eq(6.9282, get_standard_error(stat, 3), 0.0001));
        assert(eq(10.328, get_standard_error(stat, 4), 0.0001));
        stat = init_statistics(3);
        //and also test with some less systematic values
        double vals1[] = {234, 934, 423, 234, 45,6,454,45};
        update_statistics(stat, vals1);
        double vals2[] = {2,7,8,9,342,5,45,870};
        update_statistics(stat, vals2);
        double vals3[] = {3123,31,645,765,86,42,43,65};
        update_statistics(stat, vals3);
        assert(eq(352.625, get_mean(stat), 0.0001));
        assert(eq(436236, get_variance(stat, 0), 1));
        assert(eq(204564, get_variance(stat, 1), 1));
        assert(eq(177190, get_variance(stat, 2), 1));
        assert(eq(50511.3, get_variance(stat, 3), 1));
        assert(eq(134.82, get_standard_error(stat, 0), 0.001));
        assert(eq(130.564, get_standard_error(stat, 1), 0.001));
        assert(eq(171.848, get_standard_error(stat, 2), 0.001));
        assert(eq(129.758, get_standard_error(stat, 3), 0.001));
}

void test_neighbours() {
        sidelength_t L[] = {3,4,2};
        index_t correct_values[24][6] =
        {
                {1, 2, 3, 9, 12, 12}, 
                {2, 0, 4, 10, 13, 13}, 
                {0, 1, 5, 11, 14, 14}, 
                {4, 5, 6, 0, 15, 15}, 
                {5, 3, 7, 1, 16, 16}, 
                {3, 4, 8, 2, 17, 17}, 
                {7, 8, 9, 3, 18, 18}, 
                {8, 6, 10, 4, 19, 19}, 
                {6, 7, 11, 5, 20, 20}, 
                {10, 11, 0,6 , 21, 21}, 
                {11, 9, 1, 7, 22, 22}, 
                {9, 10, 2, 8, 23, 23}, 
                {13, 14, 15, 21, 0, 0}, 
                {14, 12, 16, 22, 1, 1}, 
                {12, 13, 17, 23, 2, 2}, 
                {16, 17, 18, 12, 3, 3}, 
                {17, 15, 19, 13, 4, 4}, 
                {15, 16, 20, 14, 5, 5}, 
                {19, 20, 21, 15, 6, 6}, 
                {20, 18, 22, 16, 7, 7}, 
                {18, 19, 23, 17, 8, 8}, 
                {22, 23, 12, 18, 9, 9}, 
                {23, 21, 13, 19, 10, 10}, 
                {21, 22, 14, 20, 11, 11}, 
        };
        //test the neighbours for the wolff and swendsen wang algorithm:
        wolff* w = init_wolff(L);
        swendsen_wang* sw = init_swendsen_wang(L);
        for (int i = 0; i < w->N; ++i) {
                for (int d = 0; d < D; ++d) {
                        assert((w->sites + i)->neighbours[2 * d] == correct_values[i][2 * d]);
                        assert((w->sites + i)->neighbours[2 * d + 1] == correct_values[i][2 * d + 1]);
                        assert(sw->neighbours[i][d] == correct_values[i][2 * d]);
                }
        }
        free_wolff(w);
        free_swendsen_wang(sw);
}
void test_energies() {
        //check against the values quoted in
        //W. Krauth: Statistical Mechanics Computations and Algorithms, p. 236
        sidelength_t L[D] = {6, 6, 1};
        int N = 6*6*1;
        double beta[8] = {2.0, 1.0, 0.666666666666666666666, 0.5, 0.4,
                0.333333333333333333333333, 0.2857142857, 0.25};
        double e[8] = {-1.999, -1.997, -1.951, -1.747, -1.28, -0.887, -0.683,
                -0.566};
        double c_V[8] = {0.00003, 0.02338, 0.19578, 0.68592, 1.00623,
                0.55665, 0.29617, 0.18704};
        double eps_e[8] = {0.001, 0.001, 0.001, 0.001, 0.001 ,0.002, 0.002 ,
                0.001};
        double eps_c_V[8] = {0.0001, 0.001, 0.005, 0.005, 0.005, 0.005 , 0.003,
                0.002};
        int equilibrationtime = 5000;
        int nsamples = 1 << 20;

        swendsen_wang* sw = init_swendsen_wang(L);
        wolff* w = init_wolff(L);

        binned_statistics* w_stat_energy = NULL;
        binned_statistics* sw_stat_energy = NULL;

        for(int bi = 0; bi < 8; ++bi) {
                double p = 1 - exp(-2.0 * beta[bi]);

                free(w_stat_energy);
                free(sw_stat_energy);
                w_stat_energy = init_statistics(0);
                sw_stat_energy = init_statistics(0);

                for (int64_t i = 1; i < equilibrationtime; ++i) {
                        swendsen_wang_perform_sweep(sw, p);
                }
                int64_t i = 0;
                while (i < N * equilibrationtime) {
                        i += wolff_perform_cluster_update(w, p);
                }

                for (int64_t i = 0; i < nsamples; ++i) {
                        wolff_perform_cluster_update(w, p);
                        swendsen_wang_perform_sweep(sw, p);
                        double energy;
                        // -1 because 3D gives a -1 energy per spin
                        // from exchange with itself on a n x n x 1 lattice
                        energy = wolff_measure_energy(w) + 1;
                        update_statistics(w_stat_energy, &energy);
                        energy = swendsen_wang_measure_energy(sw) + 1;
                        update_statistics(sw_stat_energy, &energy);
                }
                double w_mean = get_mean(w_stat_energy);
                double sw_mean = get_mean(sw_stat_energy);
                double w_variance = beta[bi]*beta[bi]*N*get_variance(w_stat_energy,0);
                double sw_variance = beta[bi]*beta[bi]*N*get_variance(sw_stat_energy,0);
                debug_msgf("Swendsen Wang:\t beta: %f\t <e>: %f\t c_V:%f\terror:%f\n", beta[bi], sw_mean,  sw_variance, sqrt(sw_variance/nsamples));
                debug_msgf("Wolff:\t\t beta: %f\t <e>: %f\t c_V:%f\terror:%f\n", beta[bi], w_mean,  w_variance, sqrt(w_variance/nsamples));
                debug_msgf("Difference:\t \t\t <e>: %f\t c_V:%f\n", fabs(w_mean - sw_mean),  fabs(w_variance - sw_variance));
                assert(eq(w_mean, e[bi], eps_e[bi]));
                assert(eq(sw_mean, e[bi], eps_e[bi]));
                assert(eq(w_variance, c_V[bi], eps_c_V[bi]));
                assert(eq(sw_variance, c_V[bi], eps_c_V[bi]));
        }
        free_swendsen_wang(sw);
        free(w);
}

void test_wolff_vs_swendsen_wang(){
        //test if the correlation lengths match each other
        sidelength_t L[D] = {32, 32, 32};
        int N = 32*32*32;
        int tmax = 16;
        double beta = 0.217;
        int binexp = 10;
        int lsamples = (1 << binexp) * (tmax + 1);
        int equilibrationtime = 5000;
        int nsamples = 1 << 16;

        swendsen_wang* sw = init_swendsen_wang(L);
        wolff* w = init_wolff(L);

        double* w_samples = xalloc(lsamples * sizeof(*w_samples));
        double* sw_samples = xalloc(lsamples * sizeof(*sw_samples));

        double p = 1 - exp(-2.0 * beta);

        // + 1 because we also want to sample t = 0
        size_t size = (tmax + 1) * sizeof(binned_statistics*);
        binned_statistics* w_stat_energy = init_statistics(0);
        binned_statistics* sw_stat_energy = init_statistics(0);
        binned_statistics** w_stat_corr = xalloc(size);
        binned_statistics** sw_stat_corr = xalloc(size);
        for (int i = 0; i <= tmax; ++i) {
                w_stat_corr[i] = init_statistics(binexp);
                sw_stat_corr[i] = init_statistics(binexp);
        }
clock_t tic = clock();
	for (int64_t i = 1; i < equilibrationtime; ++i) {
                swendsen_wang_perform_sweep(sw, p);
	}
        int64_t i = 0;
        while (i < N * equilibrationtime) {
                i += wolff_perform_cluster_update(w, p);
        }
clock_t toc = clock();
debug_msgf("Dropped %d sweeps to equilibrate system in %f seconds\n",
                equilibrationtime,
                (double)(toc - tic) / CLOCKS_PER_SEC);
	debug_msgf("%s", "Finished initial skirmish - start sampling and analyzing (will take some time)\n");


        assert(0 == (nsamples & (nsamples - 1))); //implementation assumes
        //that the overall number of samples to take is a power of two here
        int64_t samples_per_bin = (1 << binexp);
	for (int64_t i = 0; i < (nsamples >> binexp); ++i) {
                for (int64_t j = 0; j < samples_per_bin; ++j) {
                        wolff_perform_cluster_update(w, p);
                        swendsen_wang_perform_sweep(sw, p);
                        for (int t = 0; t <= tmax; ++t) {
                                w_samples[t * samples_per_bin + j] = wolff_measure_time_slice_correlation(w, t);
                                sw_samples[t * samples_per_bin + j] = swendsen_wang_measure_time_slice_correlation(sw, t);
                        }
                        double energy;
                        energy = wolff_measure_energy(w);
                        update_statistics(w_stat_energy, &energy);
                        energy = swendsen_wang_measure_energy(sw);
                        update_statistics(sw_stat_energy, &energy);
                        
                }
                for (int t = 0; t <= tmax; ++t) {
                update_statistics(w_stat_corr[t], w_samples + t * samples_per_bin);
                update_statistics(sw_stat_corr[t], sw_samples + t * samples_per_bin);
                }
                //update_statistics(stat[0], config->samples);

	}
        debug_msgf("%s", "Finished sampling - starting data analysis!\n");
        printf("%s\n", "#t\tb\tG_Wolff\t2sigma\tst_err\tG_S-W\t2sigma\tst_err");
        debug_msgf("Swendsen Wang:\t beta: %f\t <e>: %f\t c_V:%f\terror:%f\n", beta, get_mean(sw_stat_energy),  beta*beta*N*get_variance(sw_stat_energy,0), sqrt(get_variance(sw_stat_energy, 0)/nsamples));
        debug_msgf("Wolff:\t\t beta: %f\t <e>: %f\t c_V:%f\terror:%f\n", beta, get_mean(w_stat_energy),  beta*beta*N*get_variance(w_stat_energy,0), sqrt(get_variance(w_stat_energy, 0)/nsamples));
        debug_msgf("Difference:\t \t\t <e>: %f\t c_V:%f\n", fabs(get_mean(w_stat_energy) - get_mean(sw_stat_energy)),  beta*beta*N*fabs(get_variance(w_stat_energy,0) - get_variance(sw_stat_energy,0)));
        assert(eq(get_mean(w_stat_energy), get_mean(sw_stat_energy), 0.005));
        assert(eq(beta*beta*N*get_variance(w_stat_energy,0), beta*beta*N*get_variance(sw_stat_energy,0), 0.2));
        for (int t = 0; t <= tmax; ++t) {
                for (int b = 0; b <= binexp; ++b) {
                        printf("%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n", t, b, get_mean(w_stat_corr[t]), 2*sqrt(get_variance(w_stat_corr[t], b)), get_standard_error(w_stat_corr[t], b), get_mean(sw_stat_corr[t]), 2*sqrt(get_variance(sw_stat_corr[t], b)), get_standard_error(sw_stat_corr[t], b));
                }
        }

        //not really neccessary as we exit anyways...
        for (int t = 0; t < tmax - 1; ++t) {
                free_statistics(w_stat_corr[t]);
                free_statistics(sw_stat_corr[t]);
        }
        free(w_samples);
        free(sw_samples);
        free_swendsen_wang(sw);
        free(w);
}


int main(int argc, char** argv) {
        debug_msgf("%s\n", "Testing statistics...");
        test_statistics();

        debug_msgf("%s\n", "Testing neighbour table initialization...");
        test_neighbours();

        initialize_rand(argc == 2 ? argv[1] : DEFAULT_SEEDFILE);

        debug_msgf("%s\n", "Testing energies of 6x6x1 lattice...");
        test_energies();

        debug_msgf("%s\n", "Testing Wolff vs. Swendsen Wang implementation");
        test_wolff_vs_swendsen_wang();

        fprintf(stderr, "All static tests passed - check results of algorithm "
                        "comparison manually, please!\n");
}
