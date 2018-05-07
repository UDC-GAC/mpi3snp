
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

cuda_add_library(engine_gpu CUDAError.h GPUEngine.cpp GPUEngine.h GPUContTable.h EntropySearch.cu EntropySearch.h)
target_include_directories(engine_gpu INTERFACE .)
target_link_libraries(engine_gpu dataset engine_distributor engine_interface)