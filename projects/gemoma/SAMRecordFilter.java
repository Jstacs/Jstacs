package projects.gemoma;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;

public class SAMRecordFilter {

	private int minQuality;
	private int positionsAroundSpliceSite;
	private int maxMismatches;
	
	private HashMap<String,String> genome;

	
	public SAMRecordFilter(int minQuality, int positionsAroundSpliceSite, int maxMismatches, String targetGenome) throws Exception {
		this.minQuality = minQuality;
		this.positionsAroundSpliceSite = positionsAroundSpliceSite;
		this.maxMismatches = maxMismatches;
		if(targetGenome != null) {
			genome = Tools.getFasta(targetGenome,20,".*");
		}

	}

	
	public boolean accept(SAMRecord read) {
		
		if(read.getMappingQuality() < minQuality) {
			return false;
		}
		

		List<AlignmentBlock> blockLi = read.getAlignmentBlocks();

		Iterator<AlignmentBlock> blockIt = blockLi.iterator();

		if(genome != null && blockLi.size()>1) {

			String chr = genome.get( read.getReferenceName() );

			String rs = read.getReadString();

			int nmm = 0;

			while(blockIt.hasNext()) {
				AlignmentBlock block = blockIt.next();

				int l = block.getLength();
				int c = block.getReferenceStart()-1;
				int r = block.getReadStart()-1;

				for(int i=0;i<l/2&&i<positionsAroundSpliceSite;i++) {
					if(chr.charAt(c+i) != rs.charAt(r+i)) {
						nmm++;
					}
				}

				for(int i=l-1;i>=l/2 && i>l-positionsAroundSpliceSite-1;i--) {
					if(chr.charAt(c+i) != rs.charAt(r+i)) {
						nmm++;
					}
				}
			}
			if(nmm>maxMismatches) {
				return false;
			}
		}

		return true;

		
	}
	
	
}
