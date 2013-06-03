package org.broadinstitute.sting.utils.pairhmm;

import java.io.File;
import java.util.*;
import java.lang.reflect.*;

import org.broadinstitute.sting.utils.haplotype.Haplotype;

public final class CnyPairHMM extends PairHMM implements BatchPairHMM {
    private static class HmmInput {
	public byte[] readBases;
	public byte[] readQuals;
	public byte[] insertionGOP;
	public byte[] deletionGOP;
	public byte[] overallGCP;
	public List<Haplotype> haplotypes;
    };

    private static class ResultQueue {
	private int offset;
	private List<double[]> batchResults;
	
	public ResultQueue() {
	    batchResults = new LinkedList<double[]>();
	    offset = 0;
	}

	public void push(double[] results) {
	    batchResults.add(results);
	}
	
	public double pop() {
	    double[] results = batchResults.get(0);
	    double top = results[offset++];
	    if (offset == results.length) {
		batchResults.remove(0);
		offset = 0;
	    }
	    return top;
	}
    }

    final static String libPath = "/opt/convey/personalities/32100.1.1.1.0";
    final static String libName = "gmvhdl_gatk_hmm";

    private static boolean loaded = false;
    private List<HmmInput> batchRequests = new LinkedList<HmmInput>();
    private ResultQueue resultQueue = new ResultQueue();

    static public boolean isAvailable() {
	if (!loaded) {
	    File library = new File(libPath + "/lib" + libName + ".so");
	    return library.exists();
	}
	return true;
    }

    private native void initFpga();
    private native int dequeueRequirement(int reflen, int readlen);
    private native int enqueue(byte[] haplotypeBases,
			       byte[] readBases,
			       byte[] readQuals,
			       byte[] insertionGOP,
			       byte[] deletionGOP,
			       byte[] overallGCP,
			       int hapStartIndex,
			       boolean recacheReadValues);
    private native int flushQueue();
    private native int dequeue(double[] results);
    private native double softHmm(byte[] haplotypeBases,
				  byte[] readBases,
				  byte[] readQuals,
				  byte[] insertionGOP,
				  byte[] deletionGOP,
				  byte[] overallGCP,
				  int hapStartIndex,
				  boolean recacheReadValues);
    
    public native void reportStats();

    public void initialize( final int READ_MAX_LENGTH, final int HAPLOTYPE_MAX_LENGTH ) {
	if (!loaded) {
	    addLibraryPath(libPath);
	    System.loadLibrary(libName);
	    initFpga();
	    loaded = true;
	    System.out.println("FPGA HMM Initialized");
	}
    }

    public void batchAdd(final List<Haplotype> haplotypes, 
			 final byte[] readBases,
			 final byte[] readQuals,
			 final byte[] insertionGOP,
			 final byte[] deletionGOP,
			 final byte[] overallGCP) {
        final int numHaplotypes = haplotypes.size();
	HmmInput test = new HmmInput();
	test.readBases = readBases;
	test.readQuals = readQuals;
	test.insertionGOP = insertionGOP;
	test.deletionGOP = deletionGOP;
	test.overallGCP = overallGCP;
	test.haplotypes = haplotypes;
	batchRequests.add(test);
	for (int jjj = 0; jjj < numHaplotypes; jjj++) {
	    final boolean recacheReadValues = (jjj == 0);
	    final Haplotype haplotype = haplotypes.get(jjj);
	    enqueuePrepare(haplotype.getBases(), readBases);
	    if (enqueue(haplotype.getBases(), readBases, readQuals, insertionGOP, deletionGOP, overallGCP, 0, recacheReadValues) == 0)
		throw new RuntimeException("FPGA queue overflow in batchAdd");
	}
    }
    
    public double[] batchGetResult() {
	double[] results;

	int n = flushQueue();
	if (n > 0) {
	    results = new double[n];
	    if (dequeue(results) != n)
		System.out.println("queue underflow in enqueuePrepare");
	    resultQueue.push(results);
	}

	final HmmInput test = batchRequests.remove(0);
        final int numHaplotypes = test.haplotypes.size();
	results = new double[numHaplotypes];
	for (int jjj = 0; jjj < numHaplotypes; jjj++) {
	    results[jjj] = resultQueue.pop();
	    if (results[jjj]<-60.0) {
		final Haplotype haplotype = test.haplotypes.get(jjj);
		results[jjj]=softHmm(haplotype.getBases(), test.readBases, test.readQuals, test.insertionGOP, test.deletionGOP, test.overallGCP, 0, true);
	    }
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

    private void enqueuePrepare(byte[] haplotypeBases, byte[] readBases) {
	double[] results = null;
	int n = dequeueRequirement(haplotypeBases.length, readBases.length);
	if (n>0) {
	    results = new double[n];
	    if (dequeue(results)!=n)
		System.out.println("queue underflow in enqueuePrepare");
	} else if (n<0) {
	    n = flushQueue();
	    if (n > 0) {
		results = new double[n];
		if (dequeue(results) != n)
		    System.out.println("queue underflow in enqueuePrepare");
	    }
	}
	
	if (results != null)
	    resultQueue.push(results);
    }

    public static void addLibraryPath(String pathToAdd) {
	try {
	    final Field usrPathsField = ClassLoader.class.getDeclaredField("usr_paths");
	    usrPathsField.setAccessible(true);
 
	    //get array of paths
	    final String[] paths = (String[])usrPathsField.get(null);
	    
	    //check if the path to add is already present
	    for(String path : paths) {
		if(path.equals(pathToAdd)) {
		    return;
		}
	    }
	    
	    //add the new path
	    final String[] newPaths = Arrays.copyOf(paths, paths.length + 1);
	    newPaths[newPaths.length-1] = pathToAdd;
	    usrPathsField.set(null, newPaths);
	} catch (Exception ex) {
	}
    }
}
