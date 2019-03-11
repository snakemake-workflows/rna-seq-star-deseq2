#!/usr/bin/env python
"""
NAME: gtf2bed
=========

DESCRIPTION
===========

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

0.0.1    20190310    Initial version.

LICENCE
=======
2019, copyright Sebastian Schmeier
schmeier@tuta.io // https://www.sschmeier.com

template version: 2.0 (2018/12/19)
"""
__version__ = '0.0.1'
__date__ = '20190310'
__email__ = 'schmeier@tuta.io'
__author__ = 'Sebastian Schmeier'

import sys
import os
import argparse
import csv
import gzip
import bz2
import zipfile
import time

# non-standard lib
try:
    import gffutils
except ImportError:
    sys.stderr.write('Error: "gffutils" needed. ' + \
                     'Install with "conda install -c bioconda gffutils".\n')
    sys.exit(1)

# non-standard lib: For color handling on the shell
try:
    from colorama import init, Fore
    # INIT color
    # Initialise colours for multi-platform support.
    init()
    reset = Fore.RESET
    colors = {'success': Fore.GREEN,
              'error': Fore.RED,
              'warning': Fore.YELLOW,
              'info': ''}
except ImportError:
    sys.stderr.write('colorama lib desirable. ' + \
                     'Install with "conda install colorama".\n\n')
    reset = ''
    colors = {'success': '', 'error': '', 'warning': '', 'info': ''}


def alert(atype, text, log, repeat=False):
    if repeat:
        textout = '{} [{}] {}\r'.format(time.strftime('%Y%m%d-%H:%M:%S'),
                                        atype.rjust(7),
                                        text)
    else:
        textout = '{} [{}] {}\n'.format(time.strftime('%Y%m%d-%H:%M:%S'),
                                        atype.rjust(7),
                                        text)

    log.write('{}{}{}'.format(colors[atype], textout, reset))
    if atype == 'error':
        sys.exit(1)


def success(text, log=sys.stderr):
    alert('success', text, log)


def error(text, log=sys.stderr):
    alert('error', text, log)


def warning(text, log=sys.stderr):
    alert('warning', text, log)


def info(text, log=sys.stderr, repeat=False):
    alert('info', text, log)


def parse_cmdline():
    """ Parse command-line args. """
    # parse cmd-line ---------------------------------------------------------
    description = 'Read GTF file and convert to bed12. ' + \
                  'Only transcripts are extracted and exons will make up ' + \
                  'the bockSizes and blockStarts of bed12.'
    version = 'version {}, date {}'.format(__version__, __date__)
    epilog = 'Copyright {} ({})'.format(__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='{}'.format(version))

    parser.add_argument(
        'str_file',
        metavar='FILE',
        help='GTF-file. [use "-" or "stdin" to read from standard in]')
    parser.add_argument('--db',
                        metavar='STRING',
                        dest='annotation_db_name',
                        default="annotation.db",
                        help='Annotation database filename. [default: "annotation.db"]')
    parser.add_argument('-o',
                        '--out',
                        metavar='STRING',
                        dest='outfile_name',
                        default=None,
                        help='Out-file. [default: "stdout"]')

    # if no arguments supplied print help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args, parser


def load_file(filename):
    """ LOADING FILES """
    if filename in ['-', 'stdin']:
        filehandle = sys.stdin
    elif filename.split('.')[-1] == 'gz':
        filehandle = gzip.open(filename, 'rt')
    elif filename.split('.')[-1] == 'bz2':
        filehandle = bz2.open(filename, 'rt')
    elif filename.split('.')[-1] == 'zip':
        filehandle = zipfile.ZipFile(filename)
    else:
        filehandle = open(filename)
    return filehandle


def main():
    """ The main funtion. """
    args, parser = parse_cmdline()

    try:
        db = gffutils.create_db(args.str_file,
                                dbfn=args.annotation_db_name,
                                force=True,
                                keep_order=True,
                                merge_strategy='merge',
                                sort_attribute_values=True,
                                disable_infer_genes=True,
                                disable_infer_transcripts=True)
    except IOError:
        error('gffutils.create_db could not load file "{}". EXIT.'.format(args.str_file))

    # create outfile object
    if not args.outfile_name:
        outfileobj = sys.stdout
    elif args.outfile_name in ['-', 'stdout']:
        outfileobj = sys.stdout
    elif args.outfile_name.split('.')[-1] == 'gz':
        outfileobj = gzip.open(args.outfile_name, 'wt')
    elif args.outfile_name.split('.')[-1] == 'bz2':
        outfileobj = bz2.BZ2File(args.outfile_name, 'wt')
    else:
        outfileobj = open(args.outfile_name, 'w')

    # For printing to stdout
    # SIGPIPE is throwing exception when piping output to other tools
    # like head. => http://docs.python.org/library/signal.html
    # use a try - except clause to handle
    try:
        for tx in db.features_of_type('transcript', order_by='start'):
            bed = [s.strip() for s in db.bed12(tx).split('\t')]
            bed[3] = tx.id
            outfileobj.write('{}\n'.format('\t'.join(bed)))
        # flush output here to force SIGPIPE to be triggered
        # while inside this try block.
        sys.stdout.flush()
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shut-down
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE

    # ------------------------------------------------------
    outfileobj.close()
    return


if __name__ == '__main__':
    sys.exit(main())
