#!/usr/bin/env python3
'''
Created on Monday, May 15 2017
Last Updated Friday, Oct 20 2017

Computes molecule size and stdev (for error bounds of molecule size)

Columns: Rname Start End Size BX MI Reads Mapq_median AS_median NM_median

Version 0.0.1

@author: cjustin
'''

from optparse import OptionParser
import pysam
import statistics

class Molecule:
    def __init__(self, rname, start, end, \
                 newMolecID, barcode, \
                 interArrivals, count, \
                 mapqMedian, asMedian, nmMedian):
        
        self.rname = rname
        self.start = start
        self.end = end
        self.barcode = barcode
        self.newMolecID = newMolecID
        self.interArrivals = interArrivals
        self.count = count
        self.mapqMedian = mapqMedian
        self.asMedian = asMedian
        self.nmMedian = nmMedian

    def asTSV(self):
        return self.rname + "\t" + str(self.start+1) + "\t" + str(self.end+1) \
            + "\t" + str(self.end - self.start) + "\t" + self.barcode \
            + "\t" + str(self.newMolecID) + "\t" + str(self.count) \
            + "\t" + str(self.mapqMedian) + "\t" + str(self.asMedian) \
            + "\t" + str(self.nmMedian)

    def getLength(self):
        return self.end-self.start
        
class MolecIdentifier:
    
    def setBAM(self,bam):
        self._bam = bam
    
    def setDist(self, dist):
        self._maxDist = int(dist)
        
    def setMin(self, min):
        self._min = int(min)
    
    def setMAPQ(self, mapq):
        self._mapq = int(mapq)
    
    def setASRatio(self, asRatio):
        self._asRatio = float(asRatio)
    
    def setNM(self, nm):
        self._nm = int(nm)
    
    def setNewBam(self, filename):
        self._newBamFilename = filename
    
    def setOutput(self, filename):
        self._tsvFilename = filename
        
    def printTSV(self, molec):
        if self._tsvFilename:
            self._newMolecFH.write(molec.asTSV() + "\n")
        else:
            print(molec.asTSV())
    
    def __init__(self):
        """
        Constructor, identifies molecules based on inter-arrival time threshold
        """
        self._min = 4
        self._maxDist = 60000
        self._mapq = 1
        self._asRatio = 0.8
        self._nm = 5
        self._newBamFilename = ""
        self._tsvFilename = ""
        self._bam = ""
        
    def run(self):
        if self._bam:
            samfile = pysam.AlignmentFile(self._bam, "rb")
        else:
            samfile = pysam.AlignmentFile("-", "rb")
        
        if self._newBamFilename:
            self._outfilebam = pysam.AlignmentFile(self._newBamFilename, "wb", template=samfile)
        else:
            self._outfilebam = None
        
        header = "Rname\tStart\tEnd\tSize\tBX\tMI\tReads\tMapq_median\tAS_median\tNM_median"
        if self._tsvFilename:
            self._newMolecFH = open(self._tsvFilename, "w");
            self._newMolecFH.write(header + "\n")
        else:
            self._newMolecFH = None
            print(header)
            
        prevBarcode = ""
        prevChr = ""
        curReads = []
        trueMolecs = {}
        
        newMolecID = 0
        for read in samfile:
            barcode = ""
            if read.is_unmapped or \
            read.is_supplementary or \
            read.mapping_quality < self._mapq or \
            read.get_tag("AS") < self._asRatio*len(read.query_sequence) or \
            read.get_tag("NM") >= self._nm:
                continue
            
            # extract barcode
            barcodeList = [bc for bc in read.tags if "BX" in bc]
            if len(barcodeList) != 0:
                barcode = barcodeList[0][1]
            else:
                if self._newBamFilename:
                    self._outfilebam.write(read)
                continue
            if prevChr == "" or prevBarcode == "":
                prevBarcode = barcode
                prevChr = read.reference_id
            if prevBarcode != barcode or read.reference_id != prevChr:
                prevVal = 0
                prevRead = curReads[0]
                prevVal1 = 0
                prevVal2 = 0
                start = curReads[0].pos
                rname = curReads[0].reference_name
                interArrivals = []
                mapQs = []
                alSs = []
                noMs = []
                count = 0
                
                for curRead in curReads:                    
                    value = curRead.pos
                    absDist = value - prevVal
                    mapQs.append(curRead.mapping_quality)
                    alSs.append(curRead.get_tag("AS"))
                    noMs.append(curRead.get_tag("NM"))

                    #check if molecules should be terminated
                    if absDist > self._maxDist and prevVal > 0:
                        end = prevRead.reference_end
                        
                        #find distance from nearest read
                        molec = Molecule(rname, start, end, \
                                 newMolecID, prevBarcode, \
                                 interArrivals, count, \
                                 statistics.median(mapQs), \
                                 statistics.median(alSs), \
                                 statistics.median(noMs))
                        
                        if prevRead.is_reverse:
                            prevVal2 = value
                            prevVal1 = 0
                        else:
                            prevVal1 = value
                            prevVal2 = 0
                        start = value;
                        if count >= self._min:
                            self.printTSV(molec)
                            newMolecID += 1
                        if self._newBamFilename:
                            curRead.tags += [("MI", newMolecID)]
                            self._outfilebam.write(curRead)
                        interArrivals = []
                        mapQs = []
                        alSs = []
                        noMs = []
                        mapQs.append(curRead.mapping_quality)
                        alSs.append(curRead.get_tag("AS"))
                        noMs.append(curRead.get_tag("NM"))
                        prevVal = value
                        count = 0
                        continue
                    else:
                        if self._newBamFilename:
                            curRead.tags += [("MI", newMolecID)]
                            self._outfilebam.write(curRead)
                    
                    #inter arrival time is distance between read of the same direction
                    interArrival = 0
                    if curRead.is_reverse:
                        if prevVal2 == 0:
                            prevVal2 = value
                            prevVal = value
                            count += 1
                            continue
                        else:
                            interArrival = value - prevVal2
                            prevVal2 = value
                    else:
                        if prevVal1 == 0:
                            prevVal1 = value
                            prevVal = value
                            count += 1
                            continue
                        else:
                            interArrival = value - prevVal1
                            prevVal1 = value
                    if interArrival > 0:
                        count += 1
                        interArrivals.append(interArrival)
                    prevVal = value
                    prevRead = curRead
                end = prevRead.reference_end
                molec = Molecule(rname, start, end, \
                                 newMolecID, prevBarcode, \
                                 interArrivals, count, \
                                 statistics.median(mapQs), \
                                 statistics.median(alSs), \
                                 statistics.median(noMs))
                
                
                if count >= self._min:
                    self.printTSV(molec)
                    newMolecID += 1
                curReads = []
            curReads.append(read)
            prevBarcode = barcode
            prevChr = read.reference_id
        
        #clean up
        samfile.close()
        if self._newMolecFH != None:
            self._newMolecFH.close()
        if self._outfilebam != None:
            self._outfilebam.close()
    
if __name__ == '__main__':
    
    # specify parser options
    parser = OptionParser()
    parser.set_description("Takes a bam file via stdin and outputs molecules to a bed-like (1-based coordinates) TSV file. Read to genome bam file used must be sorted by BX tag and then by position.")
    parser.add_option("-b", "--bam", dest="bam",
                  help="Read to genome BAM file file instead of stdin (optional)", metavar="BAM")
    parser.add_option("-d", "--dist", dest="dist",
                  help="Minimum distance between reads to be considered the same molecule [60000]", metavar="DIST")
    parser.add_option("-o", "--output", dest="output",
                  help="file name of tsv file instead of stdout (optional)", metavar="OUTPUT")
    parser.add_option("-w", "--new_bam", dest="newBam",
                  help="New bam file with MI tags added (optional)", metavar="NEWBAM")
    parser.add_option("-m", "--min", dest="min",
                  help="minimum number of reads in alignment to consider (dupes are not considered) [4]", metavar="MIN")
    parser.add_option("-q", "--mapq", dest="mapq",
                  help="Reads MAPQ greater or equal to this will be kept [1]", metavar="MAPQ")
    parser.add_option("-a", "--asRatio", dest="asRatio",
                  help="Reads with an AS/Read length ratio greater or equal to this will be kept [0.8]", metavar="AS")
    parser.add_option("-n", "--nm", dest="nm",
                  help="Reads that have NM tag lower than this will be kept [5]", metavar="NM")
    
    (options, args) = parser.parse_args()  
  
    molecID = MolecIdentifier()
    if options.bam:
        molecID.setBAM(options.bam)
    if options.dist:
        molecID.setDist(options.dist)
    if options.min:
        molecID.setMin(options.min)
    if options.mapq:
        molecID.setMAPQ(options.mapq)
    if options.asRatio:
        molecID.setASRatio(options.asRatio)
    if options.nm:
        molecID.setNM(options.nm)
    if options.newBam:
        molecID.setNewBam(options.newBam)
    if options.output:
        molecID.setOutput(options.output)
    molecID.run()
