#!/usr/bin/env python
# Brett Youtsey 3/12/20

import argparse
import re
import cairo

parser = argparse.ArgumentParser(description = "mark and draw motifs on intron/exon map")
parser.add_argument("-f", "--gene_FASTA", type = str, help = "FASTA file of genes of interest. Intron: lower case bases. Exon: upper case bases. Gene names must be in headers")
parser.add_argument("-m", "--motif_file", type = str, help = "Text file with a single column of motif sequences. Degenerate bases can be used.")
args = parser.parse_args()

genes = args.gene_FASTA
motifs = args.motif_file

class gene:
    def __init__(self, header, length, exonSart, exonLength, motifList):
        self.length = length
        self.exonStart = exonStart
        self.exonLength = exonLength
        self.motifList = motifList
        self.header = header

class motif:
    def __init__(self, rootSequence, variantSequence, position):
        self.rootSequence = rootSequence
        self.variantSequence = variantSequence
        self.position = position
        # should generate the same color for each root sequence
        self.color = (rootSequence.count("g"), rootSequence.count("u"), rootSequence.count("y") )

### Define functions
def generateMotifs(entryMotif):
    '''Generates a list of all possible motifs given entry motif, recursively'''
    possibleList = list()

    if len(entryMotif) == 1:
    #BASE CASE
        if entryMotif == "y":
            possibleList = ["t", "c"]

        elif entryMotif == "u":
            possibleList = ["t"]
        else:
            possibleList = [entryMotif]

    else:
        first = entryMotif[0:int(len(entryMotif)/2)]
        firstList = generateMotifs(first)

        second = entryMotif[int(len(entryMotif)/2):int(len(entryMotif)+1)]
        secondList = generateMotifs(second)

        for firstCombo in firstList:
            for secondCombo in secondList:
                possibleList.append(firstCombo + secondCombo)

    return possibleList


def findExon(record):
    'Finds exons (Capitalized letter), returns two coordinates: leftmost exon position & exon length (bp)'
    for x in range(0,len(record)):
        if record[x].isupper() and record[x-1].islower():
            leftmost = x
        elif record[x].isupper() and record[x+1].islower():
            exonEnd = x
            exonLength = exonEnd - leftmost

    return leftmost, exonLength


def findMotifs(record, rootDict):
    '''Iterates through root dict (variant : root) and generates a list of motif classes'''
    motifList = list()

    record = record.lower()

    for aMotif in rootDict:
        regexMotif = "r" + aMotif + "/"
        if re.search(aMotif, record):
            #THE PROBLEM
            matchMotif = re.compile(aMotif)

            for m in matchMotif.finditer(record):
                motifList.append(motif(rootDict[aMotif], m.group(), m.start()))

    return motifList

### SAVE motifs
with open(motifs, "r") as motifFile:
    # dictionary with variable sequenes as key and root sequence as value
    rootDict = dict()
    for line in motifFile:
        rootSequence = line.strip().lower()
        newMotifs = generateMotifs(rootSequence)
        for aVariant in newMotifs:
            rootDict.setdefault(aVariant, rootSequence)

# list of gene classes
geneList = list()

### GO THROUGH GENES
with open(genes, "r") as geneFile:
    First = True
    for line in geneFile:
            if First:
                pastHeader = line[1:].strip()
                First = False
            # All sequence lines for a gene stored in record
            record = ""
            line = geneFile.readline().strip()
            while line[0] != ">":
                record += line
                line = geneFile.readline().strip()
                if line is "":
                    #end of file
                    break
            newHeader = line[1:]

            # THIS IS WHERE YOU FIND MOTIFS AND EXONS
            # only meant for one exon
            exonStart, exonLength = findExon(record)
            geneLength = len(record)

            # dictionary contains leftmost motif position as key and sequence as value
            motifList = findMotifs(record, rootDict)

            # save information gathered as a gene in geneList
            geneList.append(gene(pastHeader, geneLength, exonStart, exonLength, motifList))

            pastHeader = newHeader

### Draw with pycairo based on motifList, exonStart, exonLength, & geneLength
def dimScale(value):
    ''' Scales a value by the dimmensions of the pycairo image '''
    return value / dim

def drawGene(length):
    '''Draws main body of the gene (horizontal line)'''
    context.set_line_width(0.02)
    context.move_to(0.3, geneYStart)
    context.line_to(0.3 + dimScale(length), geneYStart)
    context.stroke()

def drawExon(distDiff, length):
    '''Draws exon from scaled start to scaled length (thick horizontal line)'''
    context.set_line_width(0.04)
    context.move_to(0.3 + distDiff, geneYStart)
    context.line_to(0.3 + distDiff + length, geneYStart)
    context.stroke()

def drawMotif(distDiff):
    '''Draws motif tic marks based on scaled distance of start coordinate'''
    context.set_line_width(0.003)
    context.move_to(distDiff + 0.3, geneYStart + .02)
    context.line_to(distDiff + 0.3, geneYStart - .02)
    context.stroke()

dim = 800
with cairo.SVGSurface("example.svg", dim, dim) as surface:
    context = cairo.Context(surface)
    context.scale(dim, dim)

    # amount of spacing between genes given the number of genes
    geneSpacing = ((dim-200)/len(geneList))/dim
    # starting y position for current gene
    geneYStart = geneSpacing

    # list of unique motifs to appear in legend
    motifLegend = list()

    for aGene in geneList:

        # main gene body
        geneLength = aGene.length
        context.set_source_rgba(0, 0, 0)
        drawGene(geneLength)

        # exon band
        context.set_source_rgba(0,0,0)
        drawExon(dimScale(aGene.exonStart), dimScale(aGene.exonLength))

        # motif ticks
        motifList = aGene.motifList
        for aMotif in motifList:
            context.set_source_rgba(aMotif.color[0], aMotif.color[1], aMotif.color[2])
            motifLoc = aMotif.position
            drawMotif(dimScale(motifLoc))

            motifSequence = aMotif.rootSequence
            if motifSequence not in motifLegend:
                # add unique motifs for legend
                motifLegend.append(motifSequence)

        # Write gene header
        context.set_source_rgba(0,0,0)
        context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL,
        cairo.FONT_WEIGHT_NORMAL)
        context.set_font_size(0.02)
        context.move_to(0.3,geneYStart - geneSpacing/4)
        context.show_text(aGene.header)

        geneYStart += geneSpacing

    # Generate Legend
    legendYSpacing = ((dim-400)/(len(motifLegend)+1))/dim
    legendYStart = 0.05
    context.move_to(0.05,legendYStart)
    context.show_text("Legend")

    for aRoot in motifLegend:
        context.set_source_rgba(0,0,0)
        legendYStart = legendYStart + legendYSpacing
        context.move_to(0.05, legendYStart)
        context.show_text(aRoot)

        context.set_source_rgba(aRoot.count("g"), aRoot.count("u"), aRoot.count("y"))
        context.set_line_width(0.04)
        context.move_to(0.05, legendYStart + legendYSpacing/3)
        context.line_to(0.1, legendYStart + legendYSpacing/3)
        context.stroke()
