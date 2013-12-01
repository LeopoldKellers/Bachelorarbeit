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


/*! @file common.c
 * @brief Contains the implementation of the stuff described in common.h.
 */
#include <stdlib.h>
#include <stdio.h>
#include "common.h"

void* xalloc(const size_t size) {
        void *p;
        static uint64_t total = 0;
        p = malloc(size);
        if (NULL == p) {
                //convert the total size into a bit more readable output
                char ch = ' ';
                if ((total >> 10) > 100) {
                        total = total >> 10;
                        ch = 'k';
                        if ((total >> 10) > 100) {
                                total = total >> 10;
                                ch = 'M';
                                if ((total >> 10) > 100) {
                                        total = total >> 10;
                                        ch = 'G';
                                }
                        }
                }
                fprintf(stderr, "malloc() failed trying to allocate %lu bytes "
                                "after allocating %lu %cB already!!\n",
                                size, total, ch);
                exit(EXIT_FAILURE);
        }
        total += size;
        return p;
}
