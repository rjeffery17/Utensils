'''
FileShuttle.py
=================

:Author: Rebecca
:Tags: Python

Purpose
-------


Symlink files from an original directory into new directories based on unique identifiers within file names.


--ID-file = IDs linking sequencing number to sample name containing unique identifier
--Infile-dir = Path to original files
--Outfile-ID-1 = Unique identifier within sample name that is used to symlink to a new directory
--Outfile-ID-2 = Unique identifier within sample name that is used to symlink to a new directory
--Outfile-dir-1 = Path to directory that you want to symlink files with unique ID 1 
--Outfile-dir-2 = Path to directory that you want to symlink files with unique ID 2 


Usage
-----

When you get sequencing for multiple projects, can split the files into appropriate directories based on sequencing number that corresponds to a unique ID within the sample name


Example::

   python FileShuttle.py --ID-file <WTandSampleID.txt> --Infile-dir /gfs/work/rjeffery/GFU005b/data --Outfile-ID-1 CMS --Outfile-ID-2 GFU --Outfile-dir-1 /gfs/work/rjeffery/CMSData --Outfile-dir-2 /gfs/work/rjeffery/GFUData 

Type::

   python FileShuttle.py --help

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
#### Create a dictionary that attributes  #####
######### sample name to WT number ############
###############################################
###############################################
###############################################

def CreateIDDict(ID_file):

    Identifier = open(ID_file, "r")
    IDDict = {}
    for ID in Identifier:
        IDS = ID.strip().split("\t")
        WT = IDS[0]
        Sample = IDS[1]
        IDDict[WT] = Sample

    return(IDDict)

###############################################
###############################################
###############################################
###############################################     
## Now need to go over filenames in current ###
################ directory ####################
######### and send to a new folder ############
###############################################
###############################################
###############################################

def Shuttle(IDDict, Infile_Dir, Outfile1_ID, Outfile2_ID, Outfile1_Dir, Outfile2_Dir):

    filenames = glob.glob(f'{Infile_Dir}/*.fastq.gz')
    for file in filenames:
        files = file.split("/")
        for file in files:
            for ID in IDDict:
                if ID in file:
                    if Outfile1_ID in IDDict[ID]:
                        os.symlink(f'{Infile_Dir}/{file}', f'{Outfile1_Dir}/{file}')
                    elif Outfile2_ID in IDDict[ID]: 
                        os.symlink(f'{Infile_Dir}/{file}', f'{Outfile2_Dir}/{file}')

                           
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

    parser.add_argument("--ID-file", dest="ID_file", type=str,
                        help="Supply txt file with Sample Name assigned to Sequencing Number")

    parser.add_argument("--Infile-dir", dest="Infile_Dir", type=str,
                        help="Supply path to files to be moved")

    parser.add_argument("--Outfile-ID-1", dest="Outfile1_ID", type=str,
                        help="Supply identifier that can be used to move file to a particular directory e.g. CMS")

    parser.add_argument("--Outfile-ID-2", dest="Outfile2_ID", type=str,
                        help="Supply identifier within file name that can be used to move file to a particular directory e.g. GFU")

    parser.add_argument("--Outfile-dir-1", dest="Outfile1_Dir", type=str,
                        help="Supply desired directory for file containing identifier given in --Outfile-ID-1")

    parser.add_argument("--Outfile-dir-2", dest="Outfile2_Dir", type=str,
                        help="Supply desired directory for file containing identifier given in --Outfile-ID-2")

    
    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    ###############################################
    ###############################################
    ############## Execute Functions ##############
    ###############################################
    ###############################################


    IDDict = CreateIDDict(args.ID_file)
    Shuttle(IDDict, args.Infile_Dir, args.Outfile1_ID, args.Outfile2_ID, args.Outfile1_Dir, args.Outfile2_Dir)
     
    
    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
       
