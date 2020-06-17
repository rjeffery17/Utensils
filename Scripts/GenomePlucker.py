
'''
GenomePlucker.py
=================

:Author: Rebecca
:Tags: Python

Purpose
-------

Retrieve 16S from bacteria of interest from 16S database


Usage
-----




Example::

   python GenomePlucker.py --BugIDs BugsofInterest.txt --Infile /gfs/mirror/dada2/database.fasta.gz --Outfile /gfs/work/rjeffery/MM12 

Type::

   python GenomePlucker.py --help

for command line help.

Command line options
--------------------

'''


import glob
import os
import sys
import cgatcore.experiment as E

###############################################
###############################################
###############################################
###############################################     
###### Open Bacteria of interest file,
### iterate through dada2 database and take out
# sequences based on bacteria of interest

import cgatcore.iotools as iotools
import cgat.FastaIterator as fastaiterator

def BugReader(Bacteria_IDs):
    Bugs = [x.strip("\n") for x in open(Bacteria_IDs, "r").readlines()]
    return(set(Bugs)) 

################################################
################################################

def plucker(Infile, Bugs, Outfile):
    infile = iotools.open_file(Infile)
    fastas = fastaiterator.iterate(infile)
    outfile = open(Outfile, "w")

    found = set()
    for fasta in fastas:
        if fasta.title in Bugs:
            found.add(fasta.title)
            outfile.write(">" + fasta.title + "\n" + fasta.sequence + "\n")

    notfound = Bugs.difference(found)
    if len(notfound) > 0:
        E.warn("Couldn't find %s in db" % ",".join(list(notfound)))
        
            
################################################
################################################
################################################
### Give the command line arguments which ######
######## feed into function execution ##########
################################################
################################################
################################################
            
def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--BugIDs", dest="Bacteria_IDs", type=str,
                        help= "Supply file with name of bacteria you wish to retrieve from db")

    parser.add_argument("--Db-File", dest="Infile", type=str,
                        help= "Supply file containing 16S database")

    parser.add_argument("--Outfile", dest="Outfile", type=str,
                        help= "Supply desired outfile")

    
    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    ###############################################
    ###############################################
    ############## Execute Functions ##############
    ###############################################
    ###############################################


    Bugs = BugReader(args.Bacteria_IDs)
    plucker(args.Infile, Bugs, args.Outfile)
     
    
    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
       
 
