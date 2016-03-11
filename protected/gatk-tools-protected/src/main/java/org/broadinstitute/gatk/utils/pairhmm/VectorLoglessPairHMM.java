/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
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
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system ("PHONE-HOME") which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE'S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.pairhmm;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

//For loading library from jar

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin, carneiro
 * Date: 10/16/12
 */
public class VectorLoglessPairHMM extends JNILoglessPairHMM {

    protected final static Logger logger = Logger.getLogger(VectorLoglessPairHMM.class);

    //Used to copy references to byteArrays to JNI from reads
    protected class JNIReadDataHolderClass {
        public byte[] readBases = null;
        public byte[] readQuals = null;
        public byte[] insertionGOP = null;
        public byte[] deletionGOP = null;
        public byte[] overallGCP = null;
    }

    //Used to copy references to byteArrays to JNI from haplotypes
    protected class JNIHaplotypeDataHolderClass {
        public byte[] haplotypeBases = null;
    }

    /**
     * Return 64-bit mask representing machine capabilities
     * Bit 0 is LSB, bit 63 MSB
     * Bit 0 represents sse4.1 availability
     * Bit 1 represents sse4.2 availability
     * Bit 2 represents AVX availability
     */
    public native long jniGetMachineType();

    /**
     * Function to initialize the fields of JNIReadDataHolderClass and JNIHaplotypeDataHolderClass from JVM.
     * C++ codegets FieldIDs for these classes once and re-uses these IDs for the remainder of the program. Field IDs do not
     * change per JVM session
     *
     * @param readDataHolderClass      class type of JNIReadDataHolderClass
     * @param haplotypeDataHolderClass class type of JNIHaplotypeDataHolderClass
     * @param mask                     a 64 bit integer identical to the one received from jniGetMachineType(). Users can disable usage of some hardware features by zeroing bits in the mask
     */
    private native void jniInitializeClassFieldsAndMachineMask(Class<?> readDataHolderClass, Class<?> haplotypeDataHolderClass, long mask);

    private static Boolean isVectorLoglessPairHMMLibraryLoaded = false;

    //The constructor is called only once inside PairHMMLikelihoodCalculationEngine
    public VectorLoglessPairHMM(final PairHMM.HMM_SUB_IMPLEMENTATION pairHMMSub, final boolean alwaysLoadVectorLoglessPairHMMLib) throws UserException.HardwareFeatureException {
        super();

        synchronized (isVectorLoglessPairHMMLibraryLoaded) {
            // Get the mask for the requested hardware sub-implementation
            // If a specifically requested hardware feature can not be supported, throw an exception
            long mask = pairHMMSub.getMask();
            throwIfHardwareFeatureNotSupported(mask, pairHMMSub);

            // Load the library and initialize the FieldIDs
            // Load if not loaded or if the the always load flag is true
            if (!isVectorLoglessPairHMMLibraryLoaded || alwaysLoadVectorLoglessPairHMMLib) {
                try
                {
                    //Try loading from Java's library path first
                    //Useful if someone builds his/her own library and wants to override the bundled
                    //implementation without modifying the Java code
                    System.loadLibrary("VectorLoglessPairHMM");
                    logger.info("libVectorLoglessPairHMM found in JVM library path");
                } catch (UnsatisfiedLinkError ule) {
                    //Could not load from Java's library path - try unpacking from jar
                    try
                    {
                        logger.debug("libVectorLoglessPairHMM not found in JVM library path - trying to unpack from GATK jar file");
                        loadLibraryFromJar("/org/broadinstitute/gatk/utils/pairhmm/libVectorLoglessPairHMM.so");
                        logger.info("libVectorLoglessPairHMM unpacked successfully from GATK jar file");
                    } catch (IOException ioe) {
                        //Throw the UnsatisfiedLinkError to make it clear to the user what failed
                        throw ule;
                    }
                }
                logger.info("Using vectorized implementation of PairHMM");
                isVectorLoglessPairHMMLibraryLoaded = true;

                //need to do this only once
                jniInitializeClassFieldsAndMachineMask(JNIReadDataHolderClass.class, JNIHaplotypeDataHolderClass.class, mask);
            }
        }
    }

    private native void jniInitializeHaplotypes(final int numHaplotypes, JNIHaplotypeDataHolderClass[] haplotypeDataArray);

    //Hold the mapping between haplotype and index in the list of Haplotypes passed to initialize
    //Use this mapping in computeLikelihoods to find the likelihood value corresponding to a given Haplotype
    private HashMap<Haplotype, Integer> haplotypeToHaplotypeListIdxMap = new HashMap<>();
    private JNIHaplotypeDataHolderClass[] mHaplotypeDataArray;

    @Override
    public HashMap<Haplotype, Integer> getHaplotypeToHaplotypeListIdxMap() {
        return haplotypeToHaplotypeListIdxMap;
    }

    //Used to transfer data to JNI
    //Since the haplotypes are the same for all calls to computeLikelihoods within a region, transfer the haplotypes only once to the JNI per region

    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize(final List<Haplotype> haplotypes, final Map<String, List<GATKSAMRecord>> perSampleReadList,
                           final int readMaxLength, final int haplotypeMaxLength) {
        int numHaplotypes = haplotypes.size();
        mHaplotypeDataArray = new JNIHaplotypeDataHolderClass[numHaplotypes];
        int idx = 0;
        haplotypeToHaplotypeListIdxMap.clear();
        for (final Haplotype currHaplotype : haplotypes) {
            mHaplotypeDataArray[idx] = new JNIHaplotypeDataHolderClass();
            mHaplotypeDataArray[idx].haplotypeBases = currHaplotype.getBases();
            haplotypeToHaplotypeListIdxMap.put(currHaplotype, idx);
            ++idx;
        }
        jniInitializeHaplotypes(numHaplotypes, mHaplotypeDataArray);
    }

    /**
     * Tell JNI to release arrays - really important if native code is directly accessing Java memory, if not
     * accessing Java memory directly, still important to release memory from C++
     */
    private native void jniFinalizeRegion();

    /**
     * {@inheritDoc}
     */
    @Override
    public void finalizeRegion()
    {
        jniFinalizeRegion();
    }

    /**
     * Real compute kernel
     */
    private native void jniComputeLikelihoods(int numReads, int numHaplotypes, JNIReadDataHolderClass[] readDataArray,
                                              JNIHaplotypeDataHolderClass[] haplotypeDataArray, double[] likelihoodArray, int maxNumThreadsToUse);

    /**
     * {@inheritDoc}
     */
    @Override
    public void computeLikelihoods(final ReadLikelihoods.Matrix<Haplotype> likelihoods, final List<GATKSAMRecord> processedReads, final Map<GATKSAMRecord, byte[]> gcp) {
        if (processedReads.isEmpty())
            return;
        if (doProfiling)
            startTime = System.nanoTime();
        int readListSize = processedReads.size();
        int numHaplotypes = likelihoods.alleleCount();
        JNIReadDataHolderClass[] readDataArray = new JNIReadDataHolderClass[readListSize];
        int idx = 0;
        for (GATKSAMRecord read : processedReads) {
            readDataArray[idx] = new JNIReadDataHolderClass();
            readDataArray[idx].readBases = read.getReadBases();
            readDataArray[idx].readQuals = read.getBaseQualities();
            readDataArray[idx].insertionGOP = read.getBaseInsertionQualities();
            readDataArray[idx].deletionGOP = read.getBaseDeletionQualities();
            readDataArray[idx].overallGCP = gcp.get(read);
            ++idx;
        }

        mLikelihoodArray = new double[readListSize * numHaplotypes];      //to store results
        if (doProfiling)
            threadLocalSetupTimeDiff = (System.nanoTime() - startTime);
        //for(reads)
        //   for(haplotypes)
        //       compute_full_prob()
        jniComputeLikelihoods(readListSize, numHaplotypes, readDataArray, mHaplotypeDataArray, mLikelihoodArray, 12);

        int readIdx = 0;
        for (int r = 0; r < readListSize; r++) {
            int hapIdx = 0;
            for (final Haplotype haplotype : likelihoods.alleles()) {

                //Since the order of haplotypes in the List<Haplotype> and alleleHaplotypeMap is different,
                //get idx of current haplotype in the list and use this idx to get the right likelihoodValue
                final int idxInsideHaplotypeList = haplotypeToHaplotypeListIdxMap.get(haplotype);
                likelihoods.set(hapIdx, r, mLikelihoodArray[readIdx + idxInsideHaplotypeList]);
                ++hapIdx;
            }
            readIdx += numHaplotypes;
        }
        if (doProfiling) {
            threadLocalPairHMMComputeTimeDiff = (System.nanoTime() - startTime);
            //synchronized(doProfiling)
            {
                pairHMMComputeTime += threadLocalPairHMMComputeTimeDiff;
                pairHMMSetupTime += threadLocalSetupTimeDiff;
            }
        }
    }

    /**
     * Print final profiling information from native code
     */
    public native void jniClose();

    @Override
    public void close() {
        if (doProfiling)
            logger.info("Time spent in setup for JNI call : " + (pairHMMSetupTime * 1e-9));
        super.close();
        jniClose();
    }

    //Copied from http://frommyplayground.com/how-to-load-native-jni-library-from-jar 

    /**
     * Loads library from current JAR archive
     * <p/>
     * The file from JAR is copied into system temporary directory and then loaded. The temporary file is deleted after exiting.
     * Method uses String as filename because the pathname is "abstract", not system-dependent.
     *
     * @param path The filename inside JAR as absolute path (beginning with '/'), e.g. /package/File.ext
     * @throws IOException  If temporary file creation or read/write operation fails
     * @throws IllegalArgumentException If source file (param path) does not exist
     * @throws IllegalArgumentException If the path is not absolute or if the filename is shorter than three characters (restriction of {@see File#createTempFile(java.lang.String, java.lang.String)}).
     */
    public static void loadLibraryFromJar(String path) throws IOException {

        if (!path.startsWith("/")) {
            throw new IllegalArgumentException("The path to be absolute (start with '/').");
        }

        // Obtain filename from path
        String[] parts = path.split("/");
        String filename = (parts.length > 1) ? parts[parts.length - 1] : null;

        // Split filename to prexif and suffix (extension)
        String prefix = "";
        String suffix = null;
        if (filename != null) {
            parts = filename.split("\\.", 2);
            prefix = parts[0];
            suffix = (parts.length > 1) ? "." + parts[parts.length - 1] : null; // Thanks, davs! :-)
        }

        // Check if the filename is okay
        if (filename == null || prefix.length() < 3) {
            throw new IllegalArgumentException("The filename has to be at least 3 characters long.");
        }

        // Prepare temporary file
        File temp = File.createTempFile(prefix, suffix);
        //System.out.println("Temp lib file "+temp.getAbsolutePath());
        temp.deleteOnExit();

        if (!temp.exists()) {
            throw new FileNotFoundException("File " + temp.getAbsolutePath() + " does not exist.");
        }

        // Prepare buffer for data copying
        byte[] buffer = new byte[1024];
        int readBytes;

        // Open and check input stream
        InputStream is = VectorLoglessPairHMM.class.getResourceAsStream(path);
        if (is == null) {
            throw new FileNotFoundException("File " + path + " was not found inside JAR.");
        }

        // Open output stream and copy data between source file in JAR and the temporary file
        OutputStream os = new FileOutputStream(temp);
        try {
            while ((readBytes = is.read(buffer)) != -1) {
                os.write(buffer, 0, readBytes);
            }
        } finally {
            // If read/write fails, close streams safely before throwing an exception
            os.close();
            is.close();
        }

        // Finally, load the library
        System.load(temp.getAbsolutePath());
    }

    /**
     * If the machine does not support the requested hardware feature, throw an exception
     * <p/>
     * If requesting a specific hardware feature, check if the machine supports this feature.
     * If it does not, throw an exception.
     *
     * @param mask a 64 bit integer identical to the one received from jniGetMachineType(). Users can disable usage of some hardware features by zeroing some bits in the mask.
     * @param pairHMMSub the PairHMM machine dependent sub-implementation to use for genotype likelihood calculations
     * @throws UserException.HardwareFeatureException if the hardware feature is not supported
     */
    private void throwIfHardwareFeatureNotSupported(long mask, PairHMM.HMM_SUB_IMPLEMENTATION pairHMMSub) throws UserException.HardwareFeatureException
    {
        if ( pairHMMSub.getIsSpecificHardwareRequest() ) {
            if ( !isHardwareFeatureSupported(mask) )
                throw new UserException.HardwareFeatureException("Machine does not support pairHMM hardware dependent sub-type = " + pairHMMSub);
        }
    }

    /**
     * Check if the machine supports the requested hardware feature
     * <p/>
     * Mask the bits for the hardware feature and check if they are set by the machine
     * If the bits are set, the machine supports this feature
     *
     * @param mask a 64 bit integer identical to the one received from jniGetMachineType(). Users can disable usage of some hardware features by zeroing some bits in the mask.
     * @return true of machine supports the requested hardware feature, false otherwise
     */
    private boolean isHardwareFeatureSupported(long mask)
    {
        return (mask & jniGetMachineType()) != 0x0;
    }
}
