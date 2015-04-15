package seqsounder.depthresponder;

import htsjdk.samtools.util.Interval;
import seqsounder.SiteDepth;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

public class CovFastaResponder extends DepthResponder {
    private final PrintStream output;
    private int regionOffset = 0;

    public CovFastaResponder(PrintStream output) {
        this.output = output;
    }

    public CovFastaResponder(String filename) throws FileNotFoundException {
        this.output = new PrintStream(new BufferedOutputStream(new FileOutputStream(filename)));
    }

    @Override
    public void startRegion(Interval region) {
        output.print(">" + region.getContig() + ":" + region.getStart() + "-" + region.getEnd());
    }

    @Override
    public void markDepth(SiteDepth siteDepth) {
        for (int i = siteDepth.start; i < siteDepth.end; i++) {
            output.print(regionOffset++ % 100 == 0 ? "\n" : " ");
            output.print(siteDepth.depth);
        }
    }

    @Override
    public void finishRegion(Interval region) {
        output.print("\n");
    }

    @Override
    public void finishAll() {
        output.flush();
        output.close();
    }
}
