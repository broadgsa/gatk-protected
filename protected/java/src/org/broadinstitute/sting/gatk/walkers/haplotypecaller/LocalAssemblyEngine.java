package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 14, 2011
 */
public abstract class LocalAssemblyEngine {

    public enum ASSEMBLER {
        SIMPLE_DE_BRUIJN
    }

    protected LocalAssemblyEngine() {
    }

    public abstract ArrayList<Haplotype> runLocalAssembly(ActiveRegion activeRegion, Haplotype refHaplotype, byte[] fullReferenceWithPadding, GenomeLoc refLoc, int PRUNE_FACTOR, ArrayList<VariantContext> activeAllelesToGenotype);
}
