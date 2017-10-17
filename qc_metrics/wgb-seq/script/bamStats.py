#!/usr/bin/env python
# coding:utf-8
# Author:  Richard Corbett
# Purpose: Provide a python script for iterating and calculate basic stats from a bam file
# Created: 04/12/2009

from optparse import OptionParser
import sys
import os
import re
import subprocess
import fnmatch
import operator
import math
import profile
import pstats
# import collections # Not found in all pythons on the network :(

samtools = "/path/to/your/samtools/install/samtools"
WQT=10
WsnpQThresh = 40  #parameter was never used so now fixed at 40

usage = """ %prog [options] -b bamFile.bam [-q QUALITY_THRESHOLD] [-c ]

The X coverage is calculated from non-dup, non-chastity, aligned reads.  These 
reads are multiplied by their read length.

Note that secondary alignments are currently ignored in the analysis.

"""
def main():

    parser = OptionParser(usage)
    parser.add_option("-b", "--bamfile", dest="bam",
                      action="store", type="string",
                      help="path to bam file.")
    parser.add_option("-c", dest="chip",
                      action="store_true",
                      help="Don't filter chastity reads, used for chip libraries.")  
    parser.add_option("-2", dest="second",
                      action="store_true",
                      help="Use only the second read in a pair to generate results.",
                      default=False)
    parser.add_option("-g","--genomesize", dest="gsize",
                      action="store", type=int,
                      help="Number of bases in the reference genome.  If provided, will be used to calculte X coverage. Without this parameter, the X coverage will not be printed.",
                      default=-1) 

    #get at the arguments
    (options, args) = parser.parse_args()
    bam = options.bam
    gsize = options.gsize
    if (options.bam == None):
        parser.print_help()
        sys.exit()

    #Check that our bam file exists
    if (not os.path.isfile(bam)):
        sys.stderr.write(bam + " is not a valid file\n")
        sys.exit()

    if(options.chip == None):
        keepChastity = False
    else:
        keepChastity = True

    #Slurp from our bam file
    proc = subprocess.Popen(samtools+' view '+bam,
                            shell=True,
                            stdout=subprocess.PIPE)

    print "Script_path:\t%s" % (os.path.abspath(sys.argv[0]))
    print "Bam_path:\t%s" % (os.path.abspath(bam))
    #print "Quality_threshold:\t%d" % (WsnpQThresh)  #nuked to reduce confusion.  This number is reported below on the lines is affects.
    print "Reference:\t%s" % (getReference(bam))
    
    #Counter variables
    nChastPass=0
    nChastFail=0
    nReads=0
    mapScores=[0]*11
    nMapped=0
    nUnmapped=0
    nUnAmb=0
    nAmb=0
    WnUnambThresh=0
    nUnAmb0=0  #for counting the number of mismatches in unambiguous alignments
    nUnAmb1=0
    nUnAmb2=0
    nUnAmb3=0
    nQCsumNoDups=0
    WnDistinctUniqueThresh=0
    #WnSnpBases=0
    WsnpReads=0
    chr_count={}
    nPairedAligned=0
    sumInsert=0
    WgenecovReads=0
    WmapScoreQ40=0
    nDups=0
    WnDistinctUnique=0
    bitFlag=0
    xCalcSum=0
    readLengths={}

    #compile our regex's in hopes it makes this a little faster
    m0=re.compile(r'NM:i:0')
    m1=re.compile(r'NM:i:1')
    m2=re.compile(r'NM:i:2')
    m3=re.compile(r'NM:i:3')
    m0s=re.compile(r'CM:i:0')
    m1s=re.compile(r'CM:i:1')
    m2s=re.compile(r'CM:i:2')
    m3s=re.compile(r'CM:i:3')
    c_spl=re.compile(r'.*[A-Z]')

    # a dict container for counting bases
    #bases = collections.defaultdict(int)
    
    for line in proc.stdout:
       
        linep=line.split("\t")
        bitFlag=int(linep[1])

        #Skip if secondary alignment
        if(bitFlag & 256):
            continue
        
        #Do we only want the second reads of the pairs?
        if(options.second == True and (bitFlag & 64)):
            continue
        
        nReads+=1
        if( (nReads % 1000000) == 0):
            print >> sys.stderr, "line number %d" % (nReads)

        l = len(linep[9])
        if(readLengths.has_key(l)):
            readLengths[l] +=1
        else:
            readLengths[l] = 1

        #check for chastity filtering
        if((bitFlag & 512)==512 and keepChastity == False):
            nChastFail+=1
            continue
        else:
            nChastPass+=1
            

        if((bitFlag & 4) != 4): #Not unmapped
            nMapped+=1
            #Add to the mapping score histogram
            mapScore=int(linep[4])
            if(mapScore==0):
                mapScores[0]+=1
            elif(mapScore>0 and mapScore<10):
                mapScores[1]+=1
            elif(mapScore>=10 and mapScore<20):
                mapScores[2]+=1
            elif(mapScore>=20 and mapScore<30):
                mapScores[3]+=1
            elif(mapScore>=30 and mapScore<40):
                mapScores[4]+=1
            elif(mapScore>=40 and mapScore<50):
                mapScores[5]+=1
            elif(mapScore>=50 and mapScore<60):
                mapScores[6]+=1
            elif(mapScore>=60 and mapScore<70):
                mapScores[7]+=1
            elif(mapScore>=70 and mapScore<80):
                mapScores[8]+=1
            elif(mapScore>=80 and mapScore<90):
                mapScores[9]+=1
            elif(mapScore>=90):
                mapScores[10]+=1

            if(mapScore != 0):
                nUnAmb+=1
            else:
                nAmb+=1

            if(mapScore >= WQT):
                WnUnambThresh+=1

            #Count the mismatches in the unique alignments
            if(mapScore>0):   
                if(m0.search(line) or m0s.search(line)):
                    nUnAmb0+=1
                if(m1.search(line) or m1s.search(line)):
                    nUnAmb1+=1
                if(m2.search(line) or m2s.search(line)):
                    nUnAmb2+=1
                if(m3.search(line) or m3s.search(line)):
                    nUnAmb3+=1

            #if(mapScore > WQT):
            #    WgenecovReads+=1

               
            #Not dups
            if((bitFlag & 1024)==0):
                xCalcSum+=(len(linep[9]))
                nQCsumNoDups+=1
                if(mapScore>=WQT):
                    WnDistinctUniqueThresh+=1
                if(mapScore>0):
                    key=linep[2]
                    if not chr_count.has_key(key):
                        chr_count[key] = 0
                    chr_count[key]+=1
              #      if(mapScore >=40 ):
              #          WmapScoreQ40+=1
                    
        else: #unmapped
            nUnmapped+=1

        if((bitFlag&2)==2): #paired read
            nPairedAligned+=1
            sumInsert+=abs(int(linep[8]))

        if((bitFlag&1024)==1024): #dup flag
            nDups+=1
      #  elif((bitFlag&4)!=4 and mapScore >0): #not a dup, aligned, and mapScore > 0
      #      WnDistinctUnique+=1

        # count the characters in the string and create a dictionary of char:count pairs
        #if(nReads % 100) :
        #    for c in linep[9]:
        #        bases[c] += 1

    #Now spew everything to stdout
    #print "Read_length:\t%s" % ( ",".join(map(str, readLengths)) )
    print "Read_length:\t", 
    l_string = ''
    for l in readLengths.keys():
        i_string = "%d:%d" % (l, readLengths[l])
        if(len(l_string) > 0):
            l_string = l_string + "," + i_string
        else:
            l_string = i_string
    print l_string
   
    print "Total_Number_Of_Reads:\t%d" % (nReads)
    print "Number_of_Chastity_Passed_Reads:\t%d" % (nChastPass)
    print "Number_of_Chastity_Failed:\t%d" % (nChastFail)
    print "Number_of_Duplicates:\t%d" % nDups
    print "Number_Reads_Aligned:\t%d" % (nMapped)
    print "Number_Reads_Unaligned:\t%d" % (nUnmapped)
    print "Mapping_Qualities_of_Aligned:"
    print "\t0:\t%d" % mapScores[0]
    print "\t1-9:\t%d" % mapScores[1]
    print "\t10-19:\t%d" % mapScores[2]
    print "\t20-29:\t%d" % mapScores[3]
    print "\t30-39:\t%d" % mapScores[4]
    print "\t40-49:\t%d" % mapScores[5]
    print "\t50-59:\t%d" % mapScores[6]
    print "\t60-69:\t%d" % mapScores[7]
    print "\t70-79:\t%d" % mapScores[8]
    print "\t80-89:\t%d" % mapScores[9]
    print "\t>=90:\t%d" % mapScores[10]
    print "Number_Unique_Alignments:\t%d" % nUnAmb
    print "Number_Non_Unique_Alignments:\t%d" % mapScores[0]
    print "Number_Mismatches_in_Unique_Alignments:"
    print "\t0_mismatches:\t%d" % nUnAmb0
    print "\t1_mismatches:\t%d" % nUnAmb1
    print "\t2_mismatches:\t%d" % nUnAmb2
    print "\t3_mismatches:\t%d" % nUnAmb3
    print "Number_of_Paired_Alignments:\t%d" % nPairedAligned
    if(nPairedAligned>0):
        print "Average_Insert_Size:\t%d" % (sumInsert/nPairedAligned)
    else:
        print "Average_Insert_Size:\t%f" % 0
    print "Number_of_Uniquely_Aligned_Reads_with_Q_>=_%d:\t%d" % ( WQT, WnUnambThresh)
    #print "Number_of_Uniquely_Aligned_Reads_without_Dups_Q_>_0:\t%d" % WnDistinctUnique
    print "Number_of_Uniquely_Aligned_Reads_without_Dups_and_Q_>=_%d:\t%d" % ( WQT, WnDistinctUniqueThresh )

    #Removed Feb, 2012
    #print "Number_of_Reads_Used_for_SNV_calling_with_Q_>=_%d:\t%d" % (WsnpQThresh, WsnpReads)
    
    #Number of bases deprecated -> Aug18, 2010
    #print "Number_of_Bases_Used_for_SNV_calling_with_Q_>=_%d:\t%d" % (WsnpQThresh, WnSnpBases)

    #print "Number_of_Reads_Used_for_Gene_Coverage_with_Q_>=_%d:\t%d" % (WQT, WgenecovReads )
    #print "Number_of_Uniquely_Aligned_Reads_without_Dups_and_Q_>=_40:\t%d" % ( WmapScoreQ40 )
    print "Number_Uniquely_Aligned_Without_Dups_to_Each_Chr:"
    for i in sorted(chr_count.iterkeys()):
        print "\t%s:\t%d" % (i, chr_count[i])
    if(gsize > 0) :
        print "Effective_Genome_size:\t%d" % (gsize)
        print "Estimate_for_genome_X_coverage:\t%f" % (float(xCalcSum)/float(gsize))
    #print "Relative Bases:"
    #for c in sorted(bases, key=bases.get, reverse=True):
    #    print '\t%s:%6d' % (c, bases[c])
    print "END"
    

#Determine the reference to use from the header of the bam file.
#Returns None if unavailable
def getReference(bamFile):
    rtn = None
    c = samtools + " view -H " + bamFile
    for line in subprocess.Popen(c, stdout=subprocess.PIPE, shell=True).stdout.readlines():
        #print line,
        if(re.match(r'@SQ.*NCBI-Build-36.1.*', line)):
            #print bamFile + " appears to be hg18"
            return "hg18"
        if(re.match(r'@SQ.*hg19.*', line)):
            #print bamFile + " appears to be hg19"
            return "hg19"
        if(re.match(r'@SQ.*AS:NCBI-Build-37.*', line)):
            return "hg19a"
            
    #if we get this far we didn't find either of the references we are looking for
    sys.stderr.write("Cannot determine reference from " + bamFile + "\n")
    return rtn
    
    
if __name__ == "__main__":
    #p=profile.Profile()
    #p.run("main()")
    #s = pstats.Stats(p)
    #s.sort_stats("time", "name").print_stats()

    main()

