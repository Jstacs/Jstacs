package projects.gemoma;

import java.io.File;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class ReadStats {

	private double meanSplit;
	private double sdSplit;
	private double meanReadLen;
	private double factor;
	
	public ReadStats(int minIntronLength, double factor, String... bams) {
		this.factor = factor;
		SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency( ValidationStringency.SILENT );
		
		double sum = 0;
		double sumSq = 0;
		double n = 0;
		
		double lenSum = 0;
		double lenN = 0;
		
		for(int i=0;i<bams.length;i++) {
			SamReader reader = srf.open(new File(bams[i]));

			Iterator<SAMRecord> recIt = reader.iterator();

			while(recIt.hasNext()) {
				SAMRecord curr = recIt.next();

				lenSum += curr.getReadLength();
				lenN++;

				List<AlignmentBlock> li = curr.getAlignmentBlocks();
				if(li.size() > 1) {
					Iterator<AlignmentBlock> it = li.iterator();
					AlignmentBlock al = it.next();
					while(it.hasNext()) {
						AlignmentBlock ne = it.next();
						double diff = Math.abs( ne.getReferenceStart()-al.getReferenceStart()-al.getLength())+1;
						if(diff>minIntronLength) {
							sum += diff;
							sumSq += diff*diff;
							n++;
						}
						al = ne;
					}
				}
			}
		}
		
		meanSplit = sum/n;
		sdSplit = Math.sqrt( sumSq/n - meanSplit*meanSplit );
		meanReadLen = lenSum/lenN;
		//System.out.println("stats: "+meanSplit+" "+sdSplit);
		
	}

	public double getSdSplitLength() {
		return sdSplit;
	}

	public double getMeanSplitLength() {
		return meanSplit;
	}

	public double getMeanReadLength() {
		return meanReadLen;
	}
	
	public boolean isOK(int len, int num) {
		double min = factor*(len - getMeanSplitLength())/getSdSplitLength();
		
		return num >= min;
	}
	
}
