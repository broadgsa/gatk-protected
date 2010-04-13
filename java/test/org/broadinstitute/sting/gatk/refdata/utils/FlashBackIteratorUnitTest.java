package org.broadinstitute.sting.gatk.refdata.utils;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class FlashBackIteratorUnitTest
 *         <p/>
 *         just like a greatful dead show...this will be prone to flashbacks
 */
public class FlashBackIteratorUnitTest extends BaseTest {
    private SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;

    @Before
    public void setup() {
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
    }

    @Test
    public void testBasicIteration() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 0, 0);
        FlashBackIterator iter = new FlashBackIterator(new FakeSeekableRODIterator(loc));
        GenomeLoc lastLocation = null;
        for (int x = 0; x < 10; x++) {
            iter.next();
            GenomeLoc cur = iter.position();
            if (lastLocation != null) {
                Assert.assertTrue(lastLocation.isBefore(cur));
            }
            lastLocation = cur;
        }
    }

    @Test
    public void testBasicIterationThenFlashBack() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 0, 0);
        FlashBackIterator iter = new FlashBackIterator(new FakeSeekableRODIterator(loc));
        GenomeLoc lastLocation = null;
        for (int x = 0; x < 10; x++) {
            iter.next();
            GenomeLoc cur = iter.position();
            if (lastLocation != null) {
                Assert.assertTrue(lastLocation.isBefore(cur));
            }
            lastLocation = cur;
        }
        iter.flashBackTo(GenomeLocParser.createGenomeLoc(0, 2));
    }

    @Test
    public void testBasicIterationThenFlashBackThenIterate() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 0, 0);
        FlashBackIterator iter = new FlashBackIterator(new FakeSeekableRODIterator(loc));
        GenomeLoc lastLocation = null;
        for (int x = 0; x < 10; x++) {
            iter.next();
            GenomeLoc cur = iter.position();
            if (lastLocation != null) {
                Assert.assertTrue(lastLocation.isBefore(cur));
            }
            lastLocation = cur;
        }
        iter.flashBackTo(GenomeLocParser.createGenomeLoc(0, 1));
        int count = 0;
        while (iter.hasNext()) {
            count++;
            iter.next();
        }
        Assert.assertEquals(10, count);
    }


    @Test
    public void testFlashBackTruth() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 0, 0);
        LocationAwareSeekableRODIterator backIter = new FakeSeekableRODIterator(loc);
        // remove the first three records
        backIter.next();
        backIter.next();
        backIter.next();
        FlashBackIterator iter = new FlashBackIterator(backIter);
        GenomeLoc lastLocation = null;
        for (int x = 0; x < 10; x++) {
            iter.next();
            GenomeLoc cur = iter.position();
            if (lastLocation != null) {
                Assert.assertTrue(lastLocation.isBefore(cur));
            }
            lastLocation = cur;
        }
        Assert.assertTrue(iter.canFlashBackTo(GenomeLocParser.createGenomeLoc(0, 5)));
        Assert.assertTrue(iter.canFlashBackTo(GenomeLocParser.createGenomeLoc(0, 15)));
        Assert.assertTrue(!iter.canFlashBackTo(GenomeLocParser.createGenomeLoc(0, 2)));
        Assert.assertTrue(!iter.canFlashBackTo(GenomeLocParser.createGenomeLoc(0, 1)));
    }

    @Test
    public void testBasicIterationThenFlashBackHalfWayThenIterate() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 0, 0);
        FlashBackIterator iter = new FlashBackIterator(new FakeSeekableRODIterator(loc));
        GenomeLoc lastLocation = null;
        for (int x = 0; x < 10; x++) {
            iter.next();
            GenomeLoc cur = iter.position();
            if (lastLocation != null) {
                Assert.assertTrue(lastLocation.isBefore(cur));
            }
            lastLocation = cur;
        }
        iter.flashBackTo(GenomeLocParser.createGenomeLoc(0, 5));
        int count = 0;
        while (iter.hasNext()) {
            count++;
            iter.next();
        }
        Assert.assertEquals(6, count); // chr1:5, 6, 7, 8, 9, and 10
    }
}


class FakeSeekableRODIterator implements LocationAwareSeekableRODIterator {

    // current location
    private GenomeLoc location;
    private FakeRODatum curROD;
    private int recordCount = 10;

    public FakeSeekableRODIterator(GenomeLoc startingLoc) {
        this.location = GenomeLocParser.createGenomeLoc(startingLoc.getContigIndex(), startingLoc.getStart() + 1, startingLoc.getStop() + 1);
        ;
    }

    @Override
    public GenomeLoc peekNextLocation() {
        System.err.println("Peek Next -> " + location);
        return location;
    }

    @Override
    public GenomeLoc position() {
        return location;
    }

    @Override
    public RODRecordList seekForward(GenomeLoc interval) {
        this.location = interval;
        return next();
    }

    @Override
    public boolean hasNext() {
        return (recordCount > 0);
    }

    @Override
    public RODRecordList next() {
        RODRecordList list = new FakeRODRecordList();
        curROD = new FakeRODatum("STUPIDNAME", location);
        location = GenomeLocParser.createGenomeLoc(location.getContigIndex(), location.getStart() + 1, location.getStop() + 1);
        list.add(curROD);
        recordCount--;
        return list;
    }

    @Override
    public void remove() {
        throw new IllegalStateException("GRRR");
    }
}


/** for testing only */
class FakeRODatum extends GATKFeature implements ReferenceOrderedDatum {

    final GenomeLoc location;

    public FakeRODatum(String name, GenomeLoc location) {
        super(name);
        this.location = location;
    }

    @Override
    public String getName() {
        return "false";
    }

    @Override
    public boolean parseLine(Object header, String[] parts) throws IOException {
        return false;
    }

    @Override
    public String toSimpleString() {
        return "";
    }

    @Override
    public String repl() {
        return "";
    }

    /**
     * Used by the ROD system to determine how to split input lines
     *
     * @return Regex string delimiter separating fields
     */
    @Override
    public String delimiterRegex() {
        return "";
    }

    @Override
    public GenomeLoc getLocation() {
        return location;
    }

    @Override
    public Object getUnderlyingObject() {
        return this;
    }

    @Override
    public int compareTo(ReferenceOrderedDatum that) {
        return location.compareTo(that.getLocation());
    }

    /**
     * Backdoor hook to read header, meta-data, etc. associated with the file.  Will be
     * called by the ROD system before streaming starts
     *
     * @param source source data file on disk from which this rod stream will be pulled
     *
     * @return a header object that will be passed to parseLine command
     */
    @Override
    public Object initialize(File source) throws FileNotFoundException {
        return null;
    }

    @Override
    public String getChr() {
        return location.getContig();
    }

    @Override
    public int getStart() {
        return (int)location.getStart();
    }

    @Override
    public int getEnd() {
        return (int)location.getStop();
    }
}

class FakeRODRecordList extends AbstractList<GATKFeature> implements RODRecordList {
    private final List<GATKFeature> list = new ArrayList<GATKFeature>();

    public boolean add(GATKFeature data) {
        return list.add(data);
    }

    @Override
    public GATKFeature get(int i) {
        return list.get(i);
    }

    @Override
    public int size() {
        return list.size();
    }

    @Override
    public GenomeLoc getLocation() {
        return list.get(0).getLocation();
    }

    @Override
    public String getName() {
        return "test";
    }

    @Override
    public int compareTo(RODRecordList rodRecordList) {
        return this.list.get(0).getLocation().compareTo(rodRecordList.getLocation());
    }
}