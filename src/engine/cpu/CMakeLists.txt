
# This file is part of MPI3SNP.
# Copyright (C) 2018 by Christian Ponte
# 
# MPI3SNP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# MPI3SNP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MPI3SNP. If not, see <http://www.gnu.org/licenses/>.

add_library(engine_cpu ThreadParams.h ThreadError.h CPUEngine.cpp CPUEngine.h ContTable.h EntropySearch.cpp EntropySearch.h)
target_include_directories(engine_cpu INTERFACE .)
target_link_libraries(engine_cpu PUBLIC engine_interface PRIVATE dataset engine_distributor utils)