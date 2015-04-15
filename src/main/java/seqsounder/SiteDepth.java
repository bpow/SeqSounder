package seqsounder;

public class SiteDepth {
    public final int start;
    public final int end;
    public int depth;

    public SiteDepth(int start, int end, int depth) {
        this.start = start;
        this.end = end;
        this.depth = depth;
    }

    public SiteDepth(int position, int depth) {
        this.start = position - 1;
        this.end = position;
        this.depth = depth;
    }

    @Override
    public int hashCode() { return depth; }

    public final int increment() { return ++depth; }
}

