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


/*! @file common.h
 * @brief contains some declarations useful for all other modules.
 *
 * In particular, these are a macro for convenient formatted messages that are
 * printed in debug mode with an appended #, but not in production mode and 
 * a function to allocate memory or exit the program on failure.
 */
#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>

#ifndef D
#define D 3/**<The default number of dimensions to simulate the Ising model in*/
#warning "choosing default dimensions: 3 - to disable this warning, add -DD=3"
#endif

#ifndef NDEBUG
#define DEBUGMODE 1 /**< print debug messages or not*/
#else
#define DEBUGMODE 0
#endif
/** use this in order to print formatted debug messages condotionally */
#define debug_msgf(fmt, ...) \
         /* make sure, that all output is in the right order by flushing the buffer */\
        do { if (DEBUGMODE) { fflush(stdout); fputc('#', stderr); fprintf(stderr, fmt, __VA_ARGS__); }} while (0);

typedef int32_t index_t; /**<This allows a maximum number of spins of 2^31 which
                           should be enough on today's computers but may be
                           desired higher in the future*/
typedef int16_t sidelength_t; /**< This allows a side length up to 32000, which
                                is very far away from current computer's
                                capabilities!*/



/** A convenience function to allocate memory and panic on failure.
 * Usage like malloc(), just that an error message is printed and the program
 * exited in case, malloc returns zero. Thus, return value can be assumed to
 * point to valid memory region of requested size.
 * @param size number of bytes to allocate
 * @return address of valid memory region of requested size
 */
void* xalloc(const size_t size);

#endif //COMMON_H
