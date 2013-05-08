package org.broadinstitute.sting.utils.pairhmm;

import java.util.*;

import org.broadinstitute.sting.utils.haplotype.Haplotype;

public final class CnyPairHMM extends PairHMM implements BatchPairHMM {
    private static class HmmInput {
	public List<Haplotype> haplotypes;
	public byte[] readBases;
	public byte[] readQuals;
	public byte[] insertionGOP;
	public byte[] deletionGOP;
	public byte[] overallGCP;
    };

    private static boolean loaded = false;
    private List<HmmInput> pending = new LinkedList<HmmInput>();

    static public boolean isAvailable() {
	return true;
    }

    public void batchAdd(final List<Haplotype> haplotypes, 
			 final byte[] readBases,
			 final byte[] readQuals,
			 final byte[] insertionGOP,
			 final byte[] deletionGOP,
			 final byte[] overallGCP) {
	HmmInput test=new HmmInput();
	test.haplotypes=haplotypes;
	test.readBases=readBases;
	test.readQuals=readQuals;
	test.insertionGOP=insertionGOP;
	test.deletionGOP=deletionGOP;
	test.overallGCP=overallGCP;
	pending.add(test);
    }
    
    public double[] batchResult() {
	HmmInput test=pending.remove(0);
	double[] results=new double[test.haplotypes.size()];
	for (int i=0; i<results.length; i++) {
	    results[i]=0.0;
	}
	return results;
    }

    protected double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
								  final byte[] readBases,
								  final byte[] readQuals,
								  final byte[] insertionGOP,
								  final byte[] deletionGOP,
								  final byte[] overallGCP,
								  final int hapStartIndex,
								  final boolean recacheReadValues ) {
	return 0.0;
    }     
}
