#!/usr/bin/env python

# Copyright 2020 Deniz Dogan, Ezgi Karaca
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Utility functions to parallelize execution of code.
"""

import logging

import pathos
from pathos.multiprocessing import ProcessingPool

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


def create_mp_pool(nproc=None):
    """Creates a multiprocessing pool of processes.

    Arguments
    ---------
        nproc : int, optional
            number of processors to use. Defaults to number of available CPUs
            minus 2.
    """

    n_cpu = pathos.multiprocessing.cpu_count()
    if nproc is None:
        nproc = n_cpu - 2

    assert nproc <= n_cpu, \
        f'Cannot allocate more processes than existing CPUs: {nproc} > {n_cpu}'

    return ProcessingPool(nproc)


def parallelize(mp_pool, func, finput, **kwargs):
    """Distributes execution of a function over a pool of processes.

    Keyword arguments (kwargs) will be passed to the target function.

    Arguments
    ---------
        mp_pool : ProcessingPool
            previously created pool of processes.
        func : function
            target function to execute in parallel
        finput : list
            list of input files for the target function.
    """

    ninput = len(finput)

    # Build argument lists
    arglist = []
    for argname, argvalue in kwargs.items():
        arglist.append(
            [argvalue] * ninput
        )

    logging.info(
        f'Executing {func.__name__} for {ninput} input files'
    )
    mp_pool.map(
        func,
        finput,
        *arglist
    )  # execute
