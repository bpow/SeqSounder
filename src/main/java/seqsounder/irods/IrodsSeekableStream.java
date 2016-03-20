package seqsounder.irods;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.irods.jargon.core.pub.io.FileIOOperations;
import org.irods.jargon.core.pub.io.IRODSRandomAccessFile;

import java.io.IOException;

public class IrodsSeekableStream extends SeekableStream{
    private final IRODSRandomAccessFile raf;

    public IrodsSeekableStream(IRODSRandomAccessFile raf) { this.raf = raf; }

    @Override
    public long length() {
        try {
            return raf.length();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public long position() throws IOException {
        return raf.getFilePointer();
    }

    @Override
    public void seek(long position) throws IOException {
        raf.seek(position, FileIOOperations.SeekWhenceType.SEEK_START);
    }

    @Override
    public int read() throws IOException {
        return raf.read();
    }

    @Override
    public int read(byte[] buffer, int offset, int length) throws IOException {
        return raf.read(buffer, offset, length);
    }

    @Override
    public void close() throws IOException {
        raf.close();
    }

    @Override
    public boolean eof() throws IOException {
        return raf.getFilePointer() == raf.length();
    }

    @Override
    public String getSource() {
        return null;
    }
}
