package seqsounder;

public class SiteDepth {
    public final String contig;
    public final int start;
    public final int end;
    public int depth;

    public SiteDepth(String contig, int start, int end, int depth) {
        this.contig = contig;
        this.start = start;
        this.end = end;
        this.depth = depth;
    }

    public SiteDepth(String contig, int position, int depth) {
        this.contig = contig;
        this.start = position - 1;
        this.end = position;
        this.depth = depth;
    }

    @Override
    public int hashCode() { return depth; }

    public final int increment() { return ++depth; }

    @Override
    public String toString() {
        return String.format("%s\t%d\t%s\t%d", contig, start, end, depth);
    }
}

