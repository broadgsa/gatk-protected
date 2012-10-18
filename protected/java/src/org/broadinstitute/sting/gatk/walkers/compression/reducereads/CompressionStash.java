package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.utils.GenomeLocComparator;

import java.util.TreeSet;

/**
 * A stash of regions that must be kept uncompressed in all samples
 *
 * In general, these are regions that were kept uncompressed by a tumor sample and we want to force
 * all other samples (normals and/or tumors) to also keep these regions uncompressed
 *
 * User: carneiro
 * Date: 10/15/12
 * Time: 4:08 PM
 */
public class CompressionStash extends TreeSet<SimpleGenomeLoc> {
    public CompressionStash() {
        super(new GenomeLocComparator());
    }
}
