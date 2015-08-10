/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 (“BROAD”) and the LICENSEE and is effective at the date the downloading is completed (“EFFECTIVE DATE”).
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system (“PHONE-HOME”) which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE’S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2015 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.engine.recalibration;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

/**
 * A tree of recal datum, where each contains a set of sub datum representing sub-states of the higher level one
 *
 * @author Mark DePristo
 * @since 07/27/12
 */
public class RecalDatumNode<T extends RecalDatum> {
    private final static double SMALLEST_CHI2_PVALUE = 1e-300;
    protected static final Logger logger = Logger.getLogger(RecalDatumNode.class);

    /**
     * fixedPenalty is this value if it's considered fixed
     */
    private final static double UNINITIALIZED = Double.NEGATIVE_INFINITY;

    private final T recalDatum;
    private double fixedPenalty = UNINITIALIZED;
    private final Set<RecalDatumNode<T>> subnodes;

    @Requires({"recalDatum != null"})
    public RecalDatumNode(final T recalDatum) {
        this(recalDatum, new HashSet<RecalDatumNode<T>>());
    }

    @Override
    public String toString() {
        return recalDatum.toString();
    }

    @Requires({"recalDatum != null", "subnodes != null"})
    public RecalDatumNode(final T recalDatum, final Set<RecalDatumNode<T>> subnodes) {
        this(recalDatum, UNINITIALIZED, subnodes);
    }

    @Requires({"recalDatum != null"})
    protected RecalDatumNode(final T recalDatum, final double fixedPenalty) {
        this(recalDatum, fixedPenalty, new HashSet<RecalDatumNode<T>>());
    }

    @Requires({"recalDatum != null", "subnodes != null"})
    protected RecalDatumNode(final T recalDatum, final double fixedPenalty, final Set<RecalDatumNode<T>> subnodes) {
        this.recalDatum = recalDatum;
        this.fixedPenalty = fixedPenalty;
        this.subnodes = new HashSet<RecalDatumNode<T>>(subnodes);
    }

    /**
     * Get the recal data associated with this node
     * @return
     */
    @Ensures("result != null")
    public T getRecalDatum() {
        return recalDatum;
    }

    /**
     * The set of all subnodes of this tree.  May be modified.
     * @return
     */
    @Ensures("result != null")
    public Set<RecalDatumNode<T>> getSubnodes() {
        return subnodes;
    }

    /**
     * Return the fixed penalty, if set, or else the the calculated penalty for this node
     * @return
     */
    public double getPenalty() {
        if ( fixedPenalty != UNINITIALIZED )
            return fixedPenalty;
        else
            return calcPenalty();
    }

    /**
     * Set the fixed penalty for this node to a fresh calculation from calcPenalty
     *
     * This is important in the case where you want to compute the penalty from a full
     * tree and then chop the tree up afterwards while considering the previous penalties.
     * If you don't call this function then manipulating the tree may result in the
     * penalty functions changing with changes in the tree.
     *
     * @param doEntireTree recurse into all subnodes?
     * @return the fixed penalty for this node
     */
    public double calcAndSetFixedPenalty(final boolean doEntireTree) {
        fixedPenalty = calcPenalty();
        if ( doEntireTree )
            for ( final RecalDatumNode<T> sub : subnodes )
                sub.calcAndSetFixedPenalty(doEntireTree);
        return fixedPenalty;
    }

    /**
     * Add node to the set of subnodes of this node
     * @param sub
     */
    @Requires("sub != null")
    public void addSubnode(final RecalDatumNode<T> sub) {
        subnodes.add(sub);
    }

    /**
     * Is this a leaf node (i.e., has no subnodes)?
     * @return
     */
    public boolean isLeaf() {
        return subnodes.isEmpty();
    }

    /**
     * Is this node immediately above only leaf nodes?
     *
     * @return
     */
    public boolean isAboveOnlyLeaves() {
        for ( final RecalDatumNode<T> sub : subnodes )
            if ( ! sub.isLeaf() )
                return false;
        return true;
    }

    /**
     * What's the immediate number of subnodes from this node?
     * @return
     */
    @Ensures("result >= 0")
    public int getNumSubnodes() {
        return subnodes.size();
    }

    /**
     * Total penalty is the sum of leaf node penalties
     *
     * This algorithm assumes that penalties have been fixed before pruning, as leaf nodes by
     * definition have 0 penalty unless they represent a pruned tree with underlying -- but now
     * pruned -- subtrees
     *
     * @return
     */
    public double totalPenalty() {
        if ( isLeaf() )
            return getPenalty();
        else {
            double sum = 0.0;
            for ( final RecalDatumNode<T> sub : subnodes )
                sum += sub.totalPenalty();
            return sum;
        }
    }

    /**
     * The maximum penalty among all nodes
     * @return
     */
    public double maxPenalty(final boolean leafOnly) {
        double max = ! leafOnly || isLeaf() ? getPenalty() : Double.MIN_VALUE;
        for ( final RecalDatumNode<T> sub : subnodes )
            max = Math.max(max, sub.maxPenalty(leafOnly));
        return max;
    }

    /**
     * The minimum penalty among all nodes
     * @return
     */
    public double minPenalty(final boolean leafOnly) {
        double min = ! leafOnly || isLeaf() ? getPenalty() : Double.MAX_VALUE;
        for ( final RecalDatumNode<T> sub : subnodes )
            min = Math.min(min, sub.minPenalty(leafOnly));
        return min;
    }

    /**
     * What's the longest branch from this node to any leaf?
     * @return
     */
    public int maxDepth() {
        int subMax = 0;
        for ( final RecalDatumNode<T> sub : subnodes )
            subMax = Math.max(subMax, sub.maxDepth());
        return subMax + 1;
    }

    /**
     * What's the shortest branch from this node to any leaf?  Includes this node
     * @return
     */
    @Ensures("result > 0")
    public int minDepth() {
        if ( isLeaf() )
            return 1;
        else {
            int subMin = Integer.MAX_VALUE;
            for ( final RecalDatumNode<T> sub : subnodes )
                subMin = Math.min(subMin, sub.minDepth());
            return subMin + 1;
        }
    }

    /**
     * Return the number of nodes, including this one, reachable from this node
     * @return
     */
    @Ensures("result > 0")
    public int size() {
        int size = 1;
        for ( final RecalDatumNode<T> sub : subnodes )
            size += sub.size();
        return size;
    }

    /**
     * Count the number of leaf nodes reachable from this node
     *
     * @return
     */
    @Ensures("result >= 0")
    public int numLeaves() {
        if ( isLeaf() )
            return 1;
        else {
            int size = 0;
            for ( final RecalDatumNode<T> sub : subnodes )
                size += sub.numLeaves();
            return size;
        }
    }

    /**
     * Calculate the phred-scaled p-value for a chi^2 test for independent among subnodes of this node.
     *
     * The chi^2 value indicates the degree of independence of the implied error rates among the
     * immediate subnodes
     *
     * @return the phred-scaled p-value for chi2 penalty, or 0.0 if it cannot be calculated
     */
    private double calcPenalty() {
        if ( isLeaf() || freeToMerge() )
            return 0.0;
        else if ( subnodes.size() == 1 )
            // only one value, so its free to merge away
            return 0.0;
        else {
            final long[][] counts = new long[subnodes.size()][2];

            int i = 0;
            for ( final RecalDatumNode<T> subnode : subnodes ) {
                // use the yates correction to help avoid all zeros => NaN
                counts[i][0] = Math.round(subnode.getRecalDatum().getNumMismatches()) + 1L;
                counts[i][1] = subnode.getRecalDatum().getNumObservations() + 2L;
                i++;
            }

            try {
                final double chi2PValue = new ChiSquareTestImpl().chiSquareTest(counts);
                final double penalty = -10.0 * Math.log10(Math.max(chi2PValue, SMALLEST_CHI2_PVALUE));

                // make sure things are reasonable and fail early if not
                if (Double.isInfinite(penalty) || Double.isNaN(penalty))
                    throw new ReviewedGATKException("chi2 value is " + chi2PValue + " at " + getRecalDatum());

                return penalty;
            } catch ( MathException e ) {
                throw new ReviewedGATKException("Failed in calculating chi2 value", e);
            }
        }
    }

    /**
     * Is this node free to merge because its rounded Q score is the same as all nodes below
     * @return
     */
    private boolean freeToMerge() {
        if ( isLeaf() ) // leaves are free to merge
            return true;
        else {
            final byte myQual = getRecalDatum().getEmpiricalQualityAsByte();
            for ( final RecalDatumNode<T> sub : subnodes )
                if ( sub.getRecalDatum().getEmpiricalQualityAsByte() != myQual )
                    return false;
            return true;
        }
    }

    /**
     * Calculate the penalty of this interval, given the overall error rate for the interval
     *
     * If the globalErrorRate is e, this value is:
     *
     * sum_i |log10(e_i) - log10(e)| * nObservations_i
     *
     * each the index i applies to all leaves of the tree accessible from this interval
     * (found recursively from subnodes as necessary)
     *
     * @param globalErrorRate overall error rate in real space against which we calculate the penalty
     * @return the cost of approximating the bins in this interval with the globalErrorRate
     */
    @Requires("globalErrorRate >= 0.0")
    @Ensures("result >= 0.0")
    private double calcPenaltyLog10(final double globalErrorRate) {
        if ( globalErrorRate == 0.0 ) // there were no observations, so there's no penalty
            return 0.0;

        if ( isLeaf() ) {
            // this is leave node
            return (Math.abs(Math.log10(recalDatum.getEmpiricalErrorRate()) - Math.log10(globalErrorRate))) * (double)recalDatum.getNumObservations();
            // TODO -- how we can generalize this calculation?
//            if ( this.qEnd <= minInterestingQual )
//                // It's free to merge up quality scores below the smallest interesting one
//                return 0;
//            else {
//                return (Math.abs(Math.log10(getEmpiricalErrorRate()) - Math.log10(globalErrorRate))) * getNumObservations();
//            }
        } else {
            double sum = 0;
            for ( final RecalDatumNode<T> hrd : subnodes)
                sum += hrd.calcPenaltyLog10(globalErrorRate);
            return sum;
        }
    }

    /**
     * Return a freshly allocated tree prunes to have no more than maxDepth from the root to any leaf
     *
     * @param maxDepth
     * @return
     */
    public RecalDatumNode<T> pruneToDepth(final int maxDepth) {
        if ( maxDepth < 1 )
            throw new IllegalArgumentException("maxDepth < 1");
        else {
            final Set<RecalDatumNode<T>> subPruned = new HashSet<RecalDatumNode<T>>(getNumSubnodes());
            if ( maxDepth > 1 )
                for ( final RecalDatumNode<T> sub : subnodes )
                    subPruned.add(sub.pruneToDepth(maxDepth - 1));
            return new RecalDatumNode<T>(getRecalDatum(), fixedPenalty, subPruned);
        }
    }

    /**
     * Return a freshly allocated tree with to no more than maxElements in order of penalty
     *
     * Note that nodes must have fixed penalties to this algorithm will fail.
     *
     * @param maxElements
     * @return
     */
    public RecalDatumNode<T> pruneByPenalty(final int maxElements) {
        RecalDatumNode<T> root = this;

        while ( root.size() > maxElements ) {
            // remove the lowest penalty element, and continue
            root = root.removeLowestPenaltyNode();
        }

        // our size is below the target, so we are good, return
        return root;
    }

    /**
     * Return a freshly allocated tree where all mergable nodes with < maxPenalty are merged
     *
     * Note that nodes must have fixed penalties to this algorithm will fail.
     *
     * @param maxPenaltyIn the maximum penalty we are allowed to incur for a merge
     * @param applyBonferroniCorrection if true, we will adjust penalty by the phred-scaled bonferroni correction
     *                                  for the size of the initial tree.  That is, if there are 10 nodes in the
     *                                  tree and maxPenalty is 20 we will actually enforce 10^-2 / 10 = 10^-3 = 30
     *                                  penalty for multiple testing
     * @return
     */
    public RecalDatumNode<T> pruneToNoMoreThanPenalty(final double maxPenaltyIn, final boolean applyBonferroniCorrection) {
        RecalDatumNode<T> root = this;

        final double bonferroniCorrection = 10 * Math.log10(this.size());
        final double maxPenalty = applyBonferroniCorrection ? maxPenaltyIn + bonferroniCorrection : maxPenaltyIn;

        if ( applyBonferroniCorrection )
        logger.info(String.format("Applying Bonferroni correction for %d nodes = %.2f to initial penalty %.2f for total " +
                "corrected max penalty of %.2f", this.size(), bonferroniCorrection, maxPenaltyIn, maxPenalty));

        while ( true ) {
            final Pair<RecalDatumNode<T>, Double> minPenaltyNode = root.getMinPenaltyAboveLeafNode();

            if ( minPenaltyNode == null || minPenaltyNode.getSecond() > maxPenalty ) {
                // nothing to merge, or the best candidate is above our max allowed
                if ( minPenaltyNode == null ) {
                    if ( logger.isDebugEnabled() ) logger.debug("Stopping because no candidates could be found");
                } else {
                    if ( logger.isDebugEnabled() ) logger.debug("Stopping because node " + minPenaltyNode.getFirst() + " has penalty " + minPenaltyNode.getSecond() + " > max " + maxPenalty);
                }
                break;
            } else {
                // remove the lowest penalty element, and continue
                if ( logger.isDebugEnabled() ) logger.debug("Removing node " + minPenaltyNode.getFirst() + " with penalty " + minPenaltyNode.getSecond());
                root = root.removeLowestPenaltyNode();
            }
        }

        // no more candidates exist with penalty < maxPenalty
        return root;
    }


    /**
     * Find the lowest penalty above leaf node in the tree, and return a tree without it
     *
     * Note this excludes the current (root) node
     *
     * @return
     */
    private RecalDatumNode<T> removeLowestPenaltyNode() {
        final Pair<RecalDatumNode<T>, Double> nodeToRemove = getMinPenaltyAboveLeafNode();
        if ( logger.isDebugEnabled() )
            logger.debug("Removing " + nodeToRemove.getFirst() + " with penalty " + nodeToRemove.getSecond());

        final Pair<RecalDatumNode<T>, Boolean> result = removeNode(nodeToRemove.getFirst());

        if ( ! result.getSecond() )
            throw new IllegalStateException("Never removed any node!");

        final RecalDatumNode<T> oneRemoved = result.getFirst();
        if ( oneRemoved == null )
            throw new IllegalStateException("Removed our root node, wow, didn't expect that");
        return oneRemoved;
    }

    /**
     * Finds in the tree the node with the lowest penalty whose subnodes are all leaves
     *
     * @return the node and its penalty, or null if no such node exists
     */
    private Pair<RecalDatumNode<T>, Double> getMinPenaltyAboveLeafNode() {
        if ( isLeaf() )
            // not allowed to remove leafs directly
            return null;
        if ( isAboveOnlyLeaves() )
            // we only consider removing nodes above all leaves
            return new Pair<RecalDatumNode<T>, Double>(this, getPenalty());
        else {
            // just recurse, taking the result with the min penalty of all subnodes
            Pair<RecalDatumNode<T>, Double> minNode = null;
            for ( final RecalDatumNode<T> sub : subnodes ) {
                final Pair<RecalDatumNode<T>, Double> subFind = sub.getMinPenaltyAboveLeafNode();
                if ( subFind != null && (minNode == null || subFind.getSecond() < minNode.getSecond()) ) {
                    minNode = subFind;
                }
            }
            return minNode;
        }
    }

    /**
     * Return a freshly allocated tree without the node nodeToRemove
     *
     * @param nodeToRemove
     * @return
     */
    private Pair<RecalDatumNode<T>, Boolean> removeNode(final RecalDatumNode<T> nodeToRemove) {
        if ( this == nodeToRemove ) {
            if ( isLeaf() )
                throw new IllegalStateException("Trying to remove a leaf node from the tree! " + this + " " + nodeToRemove);
            // node is the thing we are going to remove, but without any subnodes
            final RecalDatumNode<T> node = new RecalDatumNode<T>(getRecalDatum(), fixedPenalty);
            return new Pair<RecalDatumNode<T>, Boolean>(node, true);
        } else {
            // did we remove something in a sub branch?
            boolean removedSomething = false;

            // our sub nodes with the penalty node removed
            final Set<RecalDatumNode<T>> sub = new HashSet<RecalDatumNode<T>>(getNumSubnodes());

            for ( final RecalDatumNode<T> sub1 : subnodes ) {
                if ( removedSomething ) {
                    // already removed something, just add sub1 back to sub
                    sub.add(sub1);
                } else {
                    // haven't removed anything yet, so try
                    final Pair<RecalDatumNode<T>, Boolean> maybeRemoved = sub1.removeNode(nodeToRemove);
                    removedSomething = maybeRemoved.getSecond();
                    sub.add(maybeRemoved.getFirst());
                }
            }

            final RecalDatumNode<T> node = new RecalDatumNode<T>(getRecalDatum(), fixedPenalty, sub);
            return new Pair<RecalDatumNode<T>, Boolean>(node, removedSomething);
        }
    }

    /**
     * Return a collection of all of the data in the leaf nodes of this tree
     *
     * @return
     */
    public Collection<T> getAllLeaves() {
        final LinkedList<T> list = new LinkedList<T>();
        getAllLeavesRec(list);
        return list;
    }

    /**
     * Helpful recursive function for getAllLeaves()
     *
     * @param list the destination for the list of leaves
     */
    private void getAllLeavesRec(final LinkedList<T> list) {
        if ( isLeaf() )
            list.add(getRecalDatum());
        else {
            for ( final RecalDatumNode<T> sub : subnodes )
                sub.getAllLeavesRec(list);
        }
    }
}
