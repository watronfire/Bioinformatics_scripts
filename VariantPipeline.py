from csv import reader
import os
from fnmatch import fnmatch
from math import log

## Inputs
###########################################################################################
###########################################################################################

# Path to folder with input files
pathToFiles = ""

# Path to virus map
pathToVirus = ""

# Path to output folder.
pathToOutput = ""


## Generate Variant files.
###########################################################################################
###########################################################################################

# Data Structure for holding the variant information
class Variant:
    def __init__( self, position, ancestral, substitution, counts, frequency ):
        self.position = position
        self.ancestral = ancestral
        self.substitution = substitution
        self.counts = counts
        self.frequency = frequency

class confirmedVariant:
    def __init__(self, position, ancestral, substitution, coverageA, coverageB, freq1, freq2):
        if "/" not in freq2:
            self.freq2 = float( freq2[:-1] ) / 100.0
        else:
            self.freq2 = float( freq2.split( "/" )[0][:-1] ) / 100.0
        if "/" not in freq1:
            self.freq1 = float( freq1[:-1] ) / 100.0
        else:
            self.freq1 = float( freq1.split( "/" )[ 0 ][ :-1 ] ) / 100.0

        self.freqA = ( self.freq1 + self.freq2 ) / 2
        self.coverageB = coverageB
        self.coverageA = coverageA
        self.substitution = substitution
        self.ancestral = ancestral
        self.position = position
        if self.freqA != 1.0:
            self.Sn = -( ( ( 1 - self.freqA ) * log( 1 - self.freqA ) ) + ( self.freqA * log( self.freqA ) ) ) / log( 2 )
        else:
            self.Sn = 0.0

    def __str__(self):

        freq1Str = str( self.freq1 * 100 ) + "%"
        freq2Str = str( self.freq2 * 100 ) + "%"
        freqAStr = str( self.freqA * 100 ) + "%"

        attributes = [ str( self.position ), self.substitution, self.ancestral + " -> " + self.substitution, str( self.coverageB ), str( self.coverageA ), freq2Str, freq1Str, freqAStr ]
        return ",".join( attributes )


def variantParser( variantFile ):
    returnDict = dict()
    for ro in reader( variantFile ):
        # We only care about SNPs here, all other lines are disregarded.
        if "SNP" in ro[ 6 ]:
            # Calculate and assign the variant determinants.
            iposition = ro[ 1 ].replace( ",", "" )
            iancestral = ro[ 4 ][ 0 ]
            isubstitution = ro[ 0 ][ 0 ]
            icounts = ro[ 5 ].replace( ",", "" )
            ifrequency = ro[ 7 ]

            # Create a new dictionary entry with the variant information mapped to its bp position.
            returnDict[ iposition ] = Variant( iposition, iancestral, isubstitution, icounts, ifrequency )

    return returnDict

# Generate a translation table dictionary
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
aminoAcids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codonTable = dict( zip( codons, aminoAcids ) )

# Iterate throught the virus map, and load into list.
lines = list()
with open( pathToVirus, "r" ) as virusMap:
    rd = reader( virusMap )
    seq = "n"
    fl = True
    for row in rd:
        lines.append( row )
        if not fl:
            seq += row[ 5 ]
        else:
            fl = False

# Generates a list of files in the input folder. Only csv files are accepted.
inputFilesList = list()
for inputPair in os.listdir( pathToFiles ):
    if fnmatch( inputPair, "*.csv" ):
        inputFilesList.append( inputPair )

# Sort the input files.
inputFilesList.sort()

# Since files represent one of a technical replicates, files will need to be merged. Next loop will generation a list
# of paired replicates.
replicateList = list()
for i in range( len( inputFilesList ) - 1 ):
    generalName = inputFilesList[ i ][ :-5 ]
    if generalName == inputFilesList[ i + 1 ][ :-5 ]:
        replicateList.append( (inputFilesList[ i ], inputFilesList[ i + 1 ] ) )

print( inputFilesList )

# Iterate through each input file
for inputPair in replicateList:

    # Create a dictionary with the variants in the first file
    with open( pathToFiles + inputPair[0], "r" ) as inputFile:
        variantDictA = variantParser( inputFile )

    # Create a dictionary with the variants in the second file
    with open( pathToFiles + inputPair[1], "r" ) as inputFile:
        variantDictB = variantParser( inputFile )

    # Generate a final dictionary containing only entries shared between
    variantDict = dict()

    for pos in variantDictA.keys():
        if pos in variantDictB.keys():
            if variantDictA[pos].substitution == variantDictB[pos].substitution:

                # position, ancestral, substitution, coverageA, coverageB, freq1, freq2):
                variantDict[pos] = confirmedVariant( pos,
                                                     variantDictA[pos].ancestral,
                                                     variantDictA[pos].substitution,
                                                     variantDictA[pos].counts,
                                                     variantDictB[pos].counts,
                                                     variantDictA[pos].frequency,
                                                     variantDictB[pos].frequency )

    # First line of virus map is the header, so denote that first line is being read.
    firstLine = True

    # Name for output file is inputFile_output.csv
    fileName = inputPair[0][ :-5 ] + "_output.csv"

    # Creates, or overrides, an output file.
    with open( pathToOutput + fileName, "w" ) as outputFile:

        # List to hold lines before writing to file.
        writingBuffer = list()

        # Some statistical measures which will be calculated along the way.
        complexitySn = 0
        diversityNT = 0
        richness = 0
        #distance = 0
        distanceN = 0
        distanceS = 0

        for line in lines:

            # String which will hold the updated line in document.
            tempLine = ",".join( line )

            # If reading the first line then append additional column names.
            if firstLine:
                tempLine += ",var_site,variant,Change,Coverage_a,Coverage_b,var_freq_a,var_freq_b,var_freq_ave,var_3nt_aa,var_aa,N_S,Sn"
                writingBuffer.append( tempLine )
                firstLine = False
                continue


            # If the position is found in the variant dictionary, then we add the information we have for the variant to the line.
            variantFound = False
            if line[0] in variantDict:
                tempLine += "," + str( variantDict[line[0]] ) + ","
                variantFound = True

            # Else add the necessary spaces and the variants sequence.
            else:
                tempLine += ",," + line[5] + ",,,,,,,"

            # Next translate variant sequence.
            codon = ""
            if "UTR" not in line[1]:

                pos = int( line[ 0 ] )
                codonPos = int( line[2] )

                # If a variant was found then we use its substitution.
                var = ""
                if variantFound:
                    var = variantDict[line[0]].substitution
                else:
                    var = seq[pos]

                # Create the Codon at each position.
                if codonPos == 1:
                    codon = var + seq[pos+1] + seq[pos+2]
                elif codonPos == 2:
                    codon = seq[pos-1] + var + seq[pos+1]
                else:
                    codon = seq[pos-2] + seq[pos-1] + var

                # Determine whether mutation is synonymous or non-synonymous. Also going to calculate some statistics here.
                mutationType = ""
                if variantFound:

                    # Calculate specified statistics.
                    complexitySn += variantDict[ line[ 0 ] ].Sn
                    diversityNT += variantDict[ line[ 0 ] ].freqA
                    richness += 1

                    if line[6] == codonTable[codon]:
                        mutationType = "S" + "," + str( variantDict[line[0]].Sn )
                        distanceS += variantDict[line[0]].freqA
                    else:
                        mutationType = "N" + "," + str( variantDict[line[0]].Sn )
                        distanceN += variantDict[line[0]].freqA

                        # add the codon, and its translation to the line.
                tempLine += codon + "," + codonTable[codon] + "," + mutationType

            # write the line to file.
            #outputFile.write( tempLine  )
            writingBuffer.append( tempLine )

        distance = diversityNT
        complexitySn /= 10272
        diversityNT /= 10272
        selectionPN = distanceN / (distanceN + distanceS)

        writingBuffer[0] += ",,Test,region,result"
        writingBuffer[1] += ",,,,,Complexity_Sn,CDS," + str( complexitySn )
        writingBuffer[2] += ",,,,,Diversity_nt,CDS," + str( diversityNT )
        writingBuffer[3] += ",,,,,Richness,CDS," + str( richness )
        writingBuffer[4] += ",,,,,Distance,CDS," + str( distance )
        writingBuffer[5] += ",,,,,Distance,CDS_N," + str( distanceN )
        writingBuffer[6] += ",,,,,Distance,CDS_S," + str( distanceS )
        writingBuffer[7] += ",,,,,Selection_pN,CDS," + str( selectionPN )

        for entry in writingBuffer:
            outputFile.write( entry + "\n" )
