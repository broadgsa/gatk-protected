/*
 * Copyright (c) 2012 The Broad Institute
 * 
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.codecs.table;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.gatk.refdata.ReferenceDependentFeatureCodec;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;

/**
 * Reads tab deliminated tabular text files
 * 
 * <p>
 * <ul>
 * <li>Header: must begin with line HEADER or track (for IGV), followed by any
 * number of column names, separated by whitespace.</li>
 * <li>Comment lines starting with # are ignored</li>
 * <li>Each non-header and non-comment line is split into parts by whitespace,
 * and these parts are assigned as a map to their corresponding column name in
 * the header. Note that the first element (corresponding to the HEADER column)
 * must be a valid genome loc such as 1, 1:1 or 1:1-10, which is the position of
 * the Table element on the genome. Alternatively, the columns contig and
 * position or Chromosome, Start_Position, and End_Position must be present.
 * TableCodec requires that there be one value for each column in the header,
 * and no more, on all lines.</li>
 * </ul>
 * </p>
 * 
 * </p>
 * 
 * <h2>File format example</h2>
 * 
 * <pre>
 *     HEADER a b c
 *     1:1  1   2   3
 *     1:2  4   5   6
 *     1:3  7   8   9
 * </pre>
 * 
 * @author Mark DePristo
 * @since 2009
 */
public class TableCodec extends AsciiFeatureCodec<TableFeature> implements ReferenceDependentFeatureCodec {
    final static protected String commentDelimiter = "#";
    final static protected String delimiterRegex = "\\s+";
    final static protected String headerDelimiter = "HEADER";
    final static protected String igvHeaderDelimiter = "track";
    /**
     * The parser to use when resolving genome-wide locations.
     */
    protected GenomeLocParser genomeLocParser;

    protected ArrayList<String> header = new ArrayList<String>();

    private GenomeLocBuilder builder;

    public TableCodec() {
	super(TableFeature.class);
    }

    @Override
    public TableFeature decode(String line) {
	if (line.startsWith(headerDelimiter) || line.startsWith(commentDelimiter)
		|| line.startsWith(igvHeaderDelimiter))
	    return null;
	String[] split = line.split(delimiterRegex);
	if (split.length < 1)
	    throw new IllegalArgumentException("TableCodec line = " + line
		    + " doesn't appear to be a valid table format");
	return new TableFeature(builder.build(split), Arrays.asList(split), header);
    }

    @Override
    public Object readHeader(LineReader reader) {
	String line = "";
	try {
	    boolean isFirst = true;
	    while ((line = reader.readLine()) != null) {
		if (isFirst && !line.startsWith(headerDelimiter) && !line.startsWith(commentDelimiter)) {
		    throw new UserException.MalformedFile("TableCodec file does not have a header");
		}
		isFirst &= line.startsWith(commentDelimiter);
		if (line.startsWith(headerDelimiter)) {
		    if (header.size() > 0)
			throw new IllegalStateException(
				"Input table file seems to have two header lines.  The second is = " + line);
		    String spl[] = line.split(delimiterRegex);
		    for (String s : spl)
			header.add(s);
		    int contigIndex = header.indexOf("contig");
		    int positionIndex = header.indexOf("position");
		    builder = new FirstElementGenomeLocBuilder();
		    if (contigIndex != -1 && positionIndex != -1) {
			builder = new ContigGenomeLocBuilder(contigIndex, positionIndex);
		    } else {
			int chromosomeIndex = header.indexOf("Chromosome");
			int startPositionIndex = header.indexOf("Start_Position");
			int endPositionIndex = header.indexOf("End_Position");
			if (chromosomeIndex != -1 && startPositionIndex != -1 && endPositionIndex != -1) {
			    builder = new ContigStartEndGenomeLocBuilder(chromosomeIndex, startPositionIndex,
				    endPositionIndex);
			}
		    }

		    return header;
		} else if (!line.startsWith(commentDelimiter)) {
		    break;
		}
	    }
	} catch (IOException e) {
	    throw new UserException.MalformedFile("unable to parse header from TableCodec file", e);
	}
	return header;
    }

    /**
     * Set the parser to use when resolving genetic data.
     * 
     * @param genomeLocParser
     *            The supplied parser.
     */
    @Override
    public void setGenomeLocParser(GenomeLocParser genomeLocParser) {
	this.genomeLocParser = genomeLocParser;
    }

    private class ContigGenomeLocBuilder implements GenomeLocBuilder {
	private int contigIndex;
	private int positionIndex;

	public ContigGenomeLocBuilder(int contigIndex, int positionIndex) {
	    this.contigIndex = contigIndex;
	    this.positionIndex = positionIndex;
	}

	public GenomeLoc build(String[] split) {
	    int position = Integer.parseInt(split[positionIndex]);
	    return genomeLocParser.createGenomeLoc(split[contigIndex], position, position);
	}
    }

    private class ContigStartEndGenomeLocBuilder implements GenomeLocBuilder {
	private int contigIndex;
	private int endPositionIndex;
	private int startPositionIndex;

	public ContigStartEndGenomeLocBuilder(int contigIndex, int startPositionIndex, int endPositionIndex) {
	    this.contigIndex = contigIndex;
	    this.startPositionIndex = startPositionIndex;
	    this.endPositionIndex = endPositionIndex;
	}

	public GenomeLoc build(String[] split) {
	    return genomeLocParser.createGenomeLoc(split[contigIndex], Integer.parseInt(split[startPositionIndex]),
		    Integer.parseInt(split[endPositionIndex]));
	}
    }

    private class FirstElementGenomeLocBuilder implements GenomeLocBuilder {
	public GenomeLoc build(String[] split) {
	    return genomeLocParser.parseGenomeLoc(split[0]);
	}
    }

    /**
     * 
     * Interface for creating a GenomeLoc instance from a tokenized string.
     */
    private static interface GenomeLocBuilder {
	public GenomeLoc build(String[] split);
    }
}
