
'''
build.a.biome.py
=================

:Author: Rebecca
:Tags: Python

Purpose
-------

Assemble Individual bacteria 16S files into a gnotobiotic database


Usage
-----

E.g. doing 16S on data obtained from MM12 model

NB. Mapping file has to be in tsv format: filename \t taxonomy

Example::

   python BuildaBiome.py --Mapping-file=mapper.tsv  

Type::

   python BuildaBiome.py --help

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
###### Open mapping file, create dictionary to link
### filename to desired taxonomy name

import cgatcore.iotools as iotools
import cgat.FastaIterator as fastaiterator

def FileMap(mappingfile):

    file2tax = {}
    for line in open(mappingfile).readlines():
        filename, taxonomy = line.strip("\n").split("\t")
        file2tax[filename] = taxonomy
    return(file2tax)

################################################
################################################
def OutFile(outfile):
    outfile = open(outfile, "w")
    return(outfile)


def RenameFastaTitle(fastafile, file2tax, outfile):
    taxonomy = file2tax[fastafile]
    suffix=1
    for fasta in fastaiterator.iterate(iotools.open_file(fastafile)):
        suffix_str = "(" + str(suffix) + ")"
        new_title = taxonomy + suffix_str
        suffix += 1
        outfile.write(">" + new_title + "\n" + fasta.sequence + "\n")

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

    parser.add_argument("--Mapping-file", dest="mappingfile", type=str,
                        help= "Supply mapping file in tsv format filename \t taxonomy")

    parser.add_argument("--Outfile", dest="outfile", type=str,
                        help= "Desired outfile name")

    
    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    ###############################################
    ###############################################
    ############## Execute Functions ##############
    ###############################################
    ###############################################


    file2tax = FileMap(args.mappingfile)
    outfile = OutFile(args.outfile)
    for fastafile, taxonomy in file2tax.items():
        RenameFastaTitle(fastafile, file2tax, outfile)
     
    
    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
       
 
