package seqsounder;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import seqsounder.depthresponder.*;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class DepthWorker implements Runnable {
    private final int minimumMapScore;
    private final int minimumBaseQuality;
    private final boolean keepDupes;
    private final List<Interval> intervalsOfInterest;
    private final DepthResponder responder;
    private long totalReadsPaired = 0;
    private long totalPairedReadsWithMappedMates = 0;
    private long totalAlignedBases = 0;
    private long duplicateReads = 0;
    private long totalBases = -1;
    private long basesSoFar = -1;

    private final SamReader reader;

    public DepthWorker(int minimumMapScore, int minimumBaseQuality, boolean keepDupes, int histogramMaxDepth,
                       String bamFileName, List<Interval> intervalsOfInterest, DepthResponder responder) throws IOException {
        this.minimumMapScore = minimumMapScore;
        this.minimumBaseQuality = minimumBaseQuality;
        this.keepDupes = keepDupes;

        reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFileName));
        if (intervalsOfInterest.isEmpty()) {
            intervalsOfInterest = wholeGenomeIntervalsForBamFile(reader);
        }
        this.intervalsOfInterest = intervalsOfInterest; // to be really safe, could copy...
        this.responder = responder;
    }

    private List<Interval> wholeGenomeIntervalsForBamFile(SamReader reader) {
        ArrayList<Interval> intervals = new ArrayList<Interval>();
        for (SAMSequenceRecord seq : reader.getFileHeader().getSequenceDictionary().getSequences()) {
            intervals.add(new Interval(seq.getSequenceName(), 1, seq.getSequenceLength()));
        }
        return intervals;
    }

    private static final long totalBasesInIntervals(List<Interval> intervals) {
        long total = 0;
        for (Interval i : intervals) {
            total += i.length();
        }
        return total;
    }

    private int flushPrior(CoverageMap coverages, int alreadyFlushed, int barrier, Interval region) {
        if (alreadyFlushed + 1 >= barrier) return alreadyFlushed;
        for (int i = alreadyFlushed + 1; i < barrier; i++) {
            SiteDepth mutCoverage = coverages.remove(i);
            if (mutCoverage == null) {
                mutCoverage = new SiteDepth(region.getContig(), i, 0);
            }
            responder.markDepth(mutCoverage);
            basesSoFar++;
        }
        alreadyFlushed = barrier - 1;
        if (barrier >= region.getEnd() + 1) { // finishing region
            alreadyFlushed = -1;
        } else {
        }
        return alreadyFlushed;
    }

    @Override
    public void run() {
        totalBases = totalBasesInIntervals(intervalsOfInterest);
        for (Interval interval : intervalsOfInterest) {
            SAMRecordIterator query = reader.query(interval.getContig(), interval.getStart(), interval.getEnd(), false);
            responder.startRegion(interval);
            CoverageMap coverages = new CoverageMap(interval.getContig());
            int flushed = interval.getStart() - 1;
            while (query.hasNext()) {
                SAMRecord rec = query.next();
                if (rec.getMappingQuality() < minimumMapScore || rec.getReadFailsVendorQualityCheckFlag() || rec.getNotPrimaryAlignmentFlag() || rec.getReadUnmappedFlag()) {
                    continue;
                }
                if (rec.getReadPairedFlag()) {
                    totalReadsPaired++;
                    if (!rec.getMateUnmappedFlag()) totalPairedReadsWithMappedMates++;
                }
                if (rec.getDuplicateReadFlag()) {
                    duplicateReads++;
                    if (!keepDupes) continue;
                }

                byte[] baseQualities = rec.getBaseQualities();
                for (int i = 0; i < rec.getReadLength(); i++) {
                    int refPos = rec.getReferencePositionAtReadPosition(i + 1);  // htsjdk expects ONE-based read offset
                    if (refPos >= interval.getStart() && refPos <= interval.getEnd()) { // 0 is used for 'no corresponding position', e.g. insertion
                        totalAlignedBases++;
                        if (baseQualities[i] >= minimumBaseQuality) coverages.getDefault(refPos).increment();
                    }
                }
                flushed = flushPrior(coverages, flushed, rec.getAlignmentStart(), interval);
            }
            flushed = flushPrior(coverages, flushed, interval.getEnd() + 1, interval);
            responder.finishRegion(interval);
            query.close();
        }
        responder.finishAll();
        try {
            reader.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }

    public String summary() {
        StringBuilder sb = new StringBuilder();
        Formatter reportOut = new Formatter(sb, Locale.US);

        reportOut.format("Total Reads Paired:\t%d\n", totalReadsPaired);
        reportOut.format("Total Paired Reads With Mapped Mates:\t%d\n", totalPairedReadsWithMappedMates);
        reportOut.format("Duplicate Reads:\t%d\n", duplicateReads);
        reportOut.format("Total aligned bases:\t%d\n", totalAlignedBases);
/*
        long[] coverageHistogram = histogramCreator.getHistogram();
        long hist_at_least[] = new long[coverageHistogram.length];
        long cumulative = 0;
        for (int i = coverageHistogram.length - 1; i >= 0; i--) {
            cumulative += coverageHistogram[i];
            hist_at_least[i] = cumulative;
        }
        cumulative = 0;
        reportOut.println("#-----------------------------------------------------------------------------------------#");
        reportOut.println("Coverage\tCount\tCumulative Below\tCumulative Above");
        for (int i = 0; i < coverageHistogram.length; i++) {
            cumulative += coverageHistogram[i];
            reportOut.printf("%d\t%d\t%d\t%d\n", i, coverageHistogram[i], cumulative, hist_at_least[i]);
        }
        reportOut.println("#-----------------------------------------------------------------------------------------#");
        reportOut.close();
*/

        return sb.toString();
    }

    private class CoverageMap extends HashMap<Integer, SiteDepth> {
        public final String contig;
        CoverageMap(String contig) { this.contig = contig; }
        public SiteDepth getDefault(Integer key) {
            SiteDepth mi = super.get(key);
            if (null == mi) {
                mi = new SiteDepth(contig, key, 0);
                super.put(key, mi);
            }
            return mi;
        }
    }


}
