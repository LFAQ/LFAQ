/*
#	This part is done based on the source code in http://www.rob-mcculloch.org/
#
#  This program is a free software; you can redistribute it and / or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/Licenses/

*/
#ifndef GUARD_bd_h
#define GUARD_bd_h
#include "funs.h"

#ifdef MPIBART
bool bd(tree& x, xinfo& xi, pinfo& pi, RNG& gen, size_t numslaves);
#else
bool bd(std::ofstream& log, tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
#endif

#endif
