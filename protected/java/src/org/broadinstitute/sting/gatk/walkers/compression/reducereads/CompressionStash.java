package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.utils.GenomeLocComparator;

import java.util.Collection;
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

    /**
     * Adds a SimpleGenomeLoc to the stash and merges it with any overlapping (and contiguous) existing loc
     * in the stash.
     *
     * @param insertLoc the new loc to be inserted
     * @return true if the loc, or it's merged version, wasn't present in the list before.
     */
    @Override
    public boolean add(SimpleGenomeLoc insertLoc) {
        TreeSet<SimpleGenomeLoc> removedLocs = new TreeSet<SimpleGenomeLoc>();
        for (SimpleGenomeLoc existingLoc : this) {
            if (existingLoc.isPast(insertLoc)) {
                break;                                          // if we're past the loc we're done looking for overlaps.
            }
            if (existingLoc.equals(insertLoc)) {
                return false;                                   // if this loc was already present in the stash, we don't need to insert it.
            }
            if (existingLoc.contiguousP(insertLoc)) {
                removedLocs.add(existingLoc);                   // list the original loc for merging
            }
        }
        for (SimpleGenomeLoc loc : removedLocs) {
            this.remove(loc);                                   // remove all locs that will be merged
        }
        removedLocs.add(insertLoc);                             // add the new loc to the list of locs that will be merged
        return super.add(SimpleGenomeLoc.merge(removedLocs));   // merge them all into one loc and add to the stash
    }

    @Override
    public boolean addAll(Collection<? extends SimpleGenomeLoc> locs) {
        boolean result = false;
        for (SimpleGenomeLoc loc : locs) {
            result |= this.add(loc);
        }
        return result;
    }
}
