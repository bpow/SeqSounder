package seqsounder.depthresponder;

import htsjdk.samtools.util.Interval;
import seqsounder.SiteDepth;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

public class BedGraphResponder extends DepthResponder {
    private final PrintStream output;

    public BedGraphResponder(PrintStream output) {
        this.output = output;
    }

    public BedGraphResponder(String outputFile) throws FileNotFoundException {
        PrintStream ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
        this.output = ps;
    }

    @Override
    public void markDepth(SiteDepth siteDepth) {
        output.println(siteDepth.contig + "\t" + siteDepth.start + "\t" + siteDepth.end + "\t" + siteDepth.depth);
    }

    @Override
    public void finishAll() {
        output.flush();
        output.close();
    }
}
