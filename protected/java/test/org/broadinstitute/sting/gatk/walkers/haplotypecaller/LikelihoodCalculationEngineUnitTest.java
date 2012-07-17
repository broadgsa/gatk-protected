package org.broadinstitute.sting.gatk.walkers.haplotypecaller;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/14/12
 */

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.MathUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for LikelihoodCalculationEngine
 */
public class LikelihoodCalculationEngineUnitTest extends BaseTest {

    @Test
    public void testNormalizeDiploidLikelihoodMatrixFromLog10() {
        double[][] likelihoodMatrix = {
            {-90.2,     0,      0},
            {-190.1, -2.1,      0},
            {-7.0,  -17.5,  -35.9}
        };
        double[][] normalizedMatrix = {
            {-88.1,     0,      0},
            {-188.0,  0.0,      0},
            {-4.9,  -15.4,  -33.8}
        };


        Assert.assertTrue(compareDoubleArrays(LikelihoodCalculationEngine.normalizeDiploidLikelihoodMatrixFromLog10(likelihoodMatrix), normalizedMatrix));

        double[][] likelihoodMatrix2 = {
                {-90.2,     0,      0,        0},
                {-190.1, -2.1,      0,        0},
                {-7.0,  -17.5,  -35.9,        0},
                {-7.0,  -17.5,  -35.9,  -1000.0},
        };
        double[][] normalizedMatrix2 = {
                {-88.1,     0,      0,        0},
                {-188.0,  0.0,      0,        0},
                {-4.9,  -15.4,  -33.8,        0},
                {-4.9,  -15.4,  -33.8,   -997.9},
        };
        Assert.assertTrue(compareDoubleArrays(LikelihoodCalculationEngine.normalizeDiploidLikelihoodMatrixFromLog10(likelihoodMatrix2), normalizedMatrix2));
    }

    private class BasicLikelihoodTestProvider extends TestDataProvider {
        public Double readLikelihoodForHaplotype1;
        public Double readLikelihoodForHaplotype2;
        public Double readLikelihoodForHaplotype3;
        
        public BasicLikelihoodTestProvider(double a, double b) {
            super(BasicLikelihoodTestProvider.class, String.format("Diploid haplotype likelihoods for reads %f / %f",a,b));
            readLikelihoodForHaplotype1 = a;
            readLikelihoodForHaplotype2 = b;
            readLikelihoodForHaplotype3 = null;
        }

        public BasicLikelihoodTestProvider(double a, double b, double c) {
            super(BasicLikelihoodTestProvider.class, String.format("Diploid haplotype likelihoods for reads %f / %f / %f",a,b,c));
            readLikelihoodForHaplotype1 = a;
            readLikelihoodForHaplotype2 = b;
            readLikelihoodForHaplotype3 = c;
        }
        
        public double[][] expectedDiploidHaplotypeMatrix() {
            if( readLikelihoodForHaplotype3 == null ) {
                double maxValue = Math.max(readLikelihoodForHaplotype1,readLikelihoodForHaplotype2);
                double[][] normalizedMatrix = {
                        {readLikelihoodForHaplotype1 - maxValue, Double.NEGATIVE_INFINITY},
                        {Math.log10(0.5*Math.pow(10,readLikelihoodForHaplotype1) + 0.5*Math.pow(10,readLikelihoodForHaplotype2)) - maxValue, readLikelihoodForHaplotype2 - maxValue}
                };
                return normalizedMatrix;
            } else {
                double maxValue = MathUtils.max(readLikelihoodForHaplotype1,readLikelihoodForHaplotype2,readLikelihoodForHaplotype3);
                double[][] normalizedMatrix = {
                        {readLikelihoodForHaplotype1 - maxValue, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
                        {Math.log10(0.5*Math.pow(10,readLikelihoodForHaplotype1) + 0.5*Math.pow(10,readLikelihoodForHaplotype2)) - maxValue, readLikelihoodForHaplotype2 - maxValue, Double.NEGATIVE_INFINITY},
                        {Math.log10(0.5*Math.pow(10,readLikelihoodForHaplotype1) + 0.5*Math.pow(10,readLikelihoodForHaplotype3)) - maxValue,
                         Math.log10(0.5*Math.pow(10,readLikelihoodForHaplotype2) + 0.5*Math.pow(10,readLikelihoodForHaplotype3)) - maxValue, readLikelihoodForHaplotype3 - maxValue}
                };
                return normalizedMatrix;
            }
        }
        
        public double[][] calcDiploidHaplotypeMatrix() {
            ArrayList<Haplotype> haplotypes = new ArrayList<Haplotype>();
            for( int iii = 1; iii <= 3; iii++) {
                Double readLikelihood = ( iii == 1 ? readLikelihoodForHaplotype1 : ( iii == 2 ? readLikelihoodForHaplotype2 : readLikelihoodForHaplotype3) );
                if( readLikelihood != null ) {
                    Haplotype haplotype = new Haplotype( (iii == 1 ? "AAAA" : (iii == 2 ? "CCCC" : "TTTT")).getBytes() );
                    haplotype.addReadLikelihoods("myTestSample", new double[]{readLikelihood});
                    haplotypes.add(haplotype);
                }
            }
            return LikelihoodCalculationEngine.computeDiploidHaplotypeLikelihoods(haplotypes, "myTestSample");
        }
    }

    @DataProvider(name = "BasicLikelihoodTestProvider")
    public Object[][] makeBasicLikelihoodTests() {
        new BasicLikelihoodTestProvider(-1.1, -2.2);
        new BasicLikelihoodTestProvider(-2.2, -1.1);
        new BasicLikelihoodTestProvider(-1.1, -1.1);
        new BasicLikelihoodTestProvider(-9.7, -15.0);
        new BasicLikelihoodTestProvider(-1.1, -2000.2);
        new BasicLikelihoodTestProvider(-1000.1, -2.2);
        new BasicLikelihoodTestProvider(0, 0);
        new BasicLikelihoodTestProvider(-1.1, 0);
        new BasicLikelihoodTestProvider(0, -2.2);
        new BasicLikelihoodTestProvider(-100.1, -200.2);

        new BasicLikelihoodTestProvider(-1.1, -2.2, 0);
        new BasicLikelihoodTestProvider(-2.2, -1.1, 0);
        new BasicLikelihoodTestProvider(-1.1, -1.1, 0);
        new BasicLikelihoodTestProvider(-9.7, -15.0, 0);
        new BasicLikelihoodTestProvider(-1.1, -2000.2, 0);
        new BasicLikelihoodTestProvider(-1000.1, -2.2, 0);
        new BasicLikelihoodTestProvider(0, 0, 0);
        new BasicLikelihoodTestProvider(-1.1, 0, 0);
        new BasicLikelihoodTestProvider(0, -2.2, 0);
        new BasicLikelihoodTestProvider(-100.1, -200.2, 0);

        new BasicLikelihoodTestProvider(-1.1, -2.2, -12.121);
        new BasicLikelihoodTestProvider(-2.2, -1.1, -12.121);
        new BasicLikelihoodTestProvider(-1.1, -1.1, -12.121);
        new BasicLikelihoodTestProvider(-9.7, -15.0, -12.121);
        new BasicLikelihoodTestProvider(-1.1, -2000.2, -12.121);
        new BasicLikelihoodTestProvider(-1000.1, -2.2, -12.121);
        new BasicLikelihoodTestProvider(0, 0, -12.121);
        new BasicLikelihoodTestProvider(-1.1, 0, -12.121);
        new BasicLikelihoodTestProvider(0, -2.2, -12.121);
        new BasicLikelihoodTestProvider(-100.1, -200.2, -12.121);

        return BasicLikelihoodTestProvider.getTests(BasicLikelihoodTestProvider.class);
    }

    @Test(dataProvider = "BasicLikelihoodTestProvider", enabled = true)
    public void testOneReadWithTwoOrThreeHaplotypes(BasicLikelihoodTestProvider cfg) {
        double[][] calculatedMatrix = cfg.calcDiploidHaplotypeMatrix();
        double[][] expectedMatrix = cfg.expectedDiploidHaplotypeMatrix();
        logger.warn(String.format("Test: %s", cfg.toString()));
        Assert.assertTrue(compareDoubleArrays(calculatedMatrix, expectedMatrix));
    }

    /**
     * Private function to compare 2d arrays
     */
    private boolean compareDoubleArrays(double[][] b1, double[][] b2) {
        if( b1.length != b2.length ) {
            return false; // sanity check
        }

        for( int i=0; i < b1.length; i++ ){
            if( b1[i].length != b2[i].length) {
                return false; // sanity check
            }
            for( int j=0; j < b1.length; j++ ){
                if ( MathUtils.compareDoubles(b1[i][j], b2[i][j]) != 0 && !Double.isInfinite(b1[i][j]) && !Double.isInfinite(b2[i][j]))
                    return false;
            }
        }
        return true;
    }
}
