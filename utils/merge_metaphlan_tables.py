#!/usr/bin/env python

# ============================================================================== 
# Merge script: from MetaPhlAn output on single sample to a joined "clades vs samples" table
# Authors: Timothy Tickle (ttickle@hsph.harvard.edu) and Curtis Huttenhower (chuttenh@hsph.harvard.edu)
# ==============================================================================

import argparse
import csv
import os
import sys


def merge( aaastrIn, astrLabels, iCol, ostm ):
	"""
	Outputs the table join of the given pre-split string collection.
	
	:param	aaastrIn:	One or more split lines from which data are read.
	:type	aaastrIn:	collection of collections of string collections
	:param	astrLabels:	File names of input data.
	:type	astrLabels:	collection of strings
	:param	iCol:		Data column in which IDs are matched (zero-indexed).
	:type	iCol:		int
	:param	ostm:		Output stream to which matched rows are written.
	:type	ostm:		output stream

	"""
	
	setstrIDs = set()
	"""The final set of all IDs in any table."""
	ahashIDs = [{} for i in range( len( aaastrIn ) )]
	"""One hash of IDs to row numbers for each input datum."""
	aaastrData = [[] for i in range( len( aaastrIn ) )]
	"""One data table for each input datum."""
	aastrHeaders = [[] for i in range( len( aaastrIn ) )]
	"""The list of non-ID headers for each input datum."""
	strHeader = "ID"
	"""The ID column header."""

	# For each input datum in each input stream...
	pos = 0

	for f in aaastrIn :
		with open(f) as csvfile :
			iIn = csv.reader(csvfile, csv.excel_tab)

			# Lines from the current file, empty list to hold data, empty hash to hold ids
			aastrData, hashIDs = (a[pos] for a in (aaastrData, ahashIDs))

			iLine = -1
			# For a line in the file
			for astrLine in iIn:
				iLine += 1

				# ID is from first column, data are everything else
				strID, astrData = astrLine[iCol], ( astrLine[:iCol] + astrLine[( iCol + 1 ):] )

				hashIDs[strID] = iLine
				aastrData.append( astrData )

			# Batch merge every new ID key set
			setstrIDs.update( hashIDs.keys( ) )

		pos += 1

	# Create writer
	csvw = csv.writer( ostm, csv.excel_tab, lineterminator='\n' )

	# Make the file names the column names
	csvw.writerow( [strHeader] + [os.path.splitext(f)[0] for f in astrLabels] )

	# Write out data
	for strID in sorted( setstrIDs ):
		astrOut = []
		for iIn in range( len( aaastrIn ) ):
			aastrData, hashIDs = (a[iIn] for a in (aaastrData, ahashIDs))
			# Look up the row number of the current ID in the current dataset, if any
			iID = hashIDs.get( strID )
			# If not, start with no data; if yes, pull out stored data row
			astrData = [0.0] if ( iID == None ) else aastrData[iID]
			# Pad output data as needed
			astrData += [None] * ( len( aastrHeaders[iIn] ) - len( astrData ) )
			astrOut += astrData
		csvw.writerow( [strID] + astrOut )


argp = argparse.ArgumentParser( prog = "merge_metaphlan_tables.py",
	description = """Performs a table join on one or more metaphlan output files.""")
argp.add_argument( "aistms",	metavar = "input.txt", nargs = "+",
	help = "One or more tab-delimited text tables to join" )

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" )

argp.usage = argp.format_usage()[7:]+"\n\n\tPlease make sure to supply file paths to the files to combine. If combining 3 files (Table1.txt, Table2.txt, and Table3.txt) the call should be:\n\n\t\tpython merge_metaphlan_tables.py Table1.txt Table2.txt Table3.txt > output.txt\n\n\tA wildcard to indicate all .txt files that start with Table can be used as follows:\n\n\t\tpython merge_metaphlan_tables.py Table*.txt > output.txt"


def _main( ):
	args = argp.parse_args( )
	merge(args.aistms, [os.path.split(os.path.basename(f))[1] for f in args.aistms], 0, sys.stdout)


if __name__ == "__main__":
	_main( )
