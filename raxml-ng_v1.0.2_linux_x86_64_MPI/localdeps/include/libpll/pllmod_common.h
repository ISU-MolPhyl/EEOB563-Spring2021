/*
 Copyright (C) 2016-18 Diego Darriba, Alexey Kozlov

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

#ifndef PLLMOD_COMMON_H_
#define PLLMOD_COMMON_H_

#define PLLMOD_COMMON_BRLEN_LINKED    0
#define PLLMOD_COMMON_BRLEN_SCALED    1
#define PLLMOD_COMMON_BRLEN_UNLINKED  2

#define PLLMOD_COMMON_REDUCE_SUM     0
#define PLLMOD_COMMON_REDUCE_MAX     1
#define PLLMOD_COMMON_REDUCE_MIN     2

#define PLLMOD_UNUSED(expr) do { (void)(expr); } while (0)

#define PLLMOD_ERRMSG_LEN 200

/* common error codes for all libpll modules (1001-2000)*/
#define PLLMOD_ERROR_INVALID_RANGE                1001
#define PLLMOD_ERROR_INVALID_NODE_TYPE            1002
#define PLLMOD_ERROR_INVALID_INDEX                1003
#define PLLMOD_ERROR_NOT_IMPLEMENTED              1004

void pllmod_set_error(int errno, const char* errmsg_fmt, ...);
void pllmod_reset_error();

#endif
