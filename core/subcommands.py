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
Command-line subcommand executors for minnie top-level script.
"""

import logging

from .parallel import (
    create_mp_pool,
    parallelize
)


# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


def splitpdbs(args):
    """Executes splitpdbs subcommand.

    Arguments
    ---------
        args: argparse.Namespace
            Result of parsing user arguments with argparse.
    """

    from .analysis import split_pdbs

    # Validate files exist
    args.pdbs = [
        f.resolve(strict=True) for f in args.pdbs
    ]

    # Validate project ids
    if args.project_ids:
        assert len(args.pdbs) == len(args.project_ids), \
            'Number of trajectories does not match number of ids: ' + \
            f'{len(args.pdbs)} != {len(args.project_ids)}'
    else:
        logging.info(f'Auto-assigning project ids from input file names')
        args.project_ids = [
            f.stem for f in args.pdbs
        ]

    for pdbfile, project_id in zip(args.pdbs, args.project_ids):
        split_pdbs(pdbfile, project_id)


def findbonds(args):
    """Executes findbonds command.

    Arguments
    ---------
        args: argparse.Namespace
            Result of parsing user arguments with argparse.
    """

    from .analysis import (
        comb_int,
        combine_interfacea_results
    )

    _itypes = [
        'hbonds',
        'ionic',
        'hydrophobic',
        'ring_stacking'
    ]

    if args.pdbfile:  # single structure
        pdblist = [args.pdbfile.resolve(strict=True)]
    elif args.folder:
        folder = args.folder.resolve(strict=True)
        pdblist = [
            f for f in folder.rglob('*.pdb')
        ]
    else:
        raise NotImplementedError(
            f'input namespace requires "pdbfile" or "folder" parameter'
        )

    if not args.project_id:
        logging.info(f'Auto-assigning project id from input file names')
        args.project_id = (args.pdbfile or args.folder).stem

    if 'all' in args.itypes:
        args.types = _itypes

    mp_pool = create_mp_pool(args.nproc)
    for itype in args.itypes:
        parallelize(
            mp_pool,
            comb_int,
            pdblist,
            project_id=args.project_id,
            itypes=itype,
            include_intra=args.intra
        )

    combine_interfacea_results(
        args.project_id,
        args.clean
    )

    mp_pool.close()
