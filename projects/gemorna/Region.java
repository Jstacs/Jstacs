package projects.gemorna;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.function.Predicate;

import de.jstacs.utils.IntList;
import de.jstacs.utils.ToolBox;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import projects.gemoma.ExtractRNAseqEvidence.Stranded;

public class Region {

	private static Random r = new Random(113);
	private double maxcov;
	private double sample;
	
	private LinkedList<SAMRecord> list;
	private HashMap<Integer,Character> strands;
	private HashMap<String,Integer> idMap;
	private Integer regionStart;
	private Integer regionEnd;
	private String chrom;
	private double sampleProb;
	private double nOut;
	
	private int[] coverageByBlocks;
	
	private Region revRegion;
	
	private char strand;
	
	public Region(double maxcov, double sample, char strand) {
		this.list = new LinkedList<SAMRecord>();
		this.chrom = null;
		this.regionStart = null;
		this.regionEnd = null;
		this.sampleProb = 1.0;
		this.sample = sample;
		this.maxcov = maxcov;
		this.nOut = 0;
		this.strand = strand;
	}
	
	public char getStrand() {
		return strand;
	}
	
	public void setRevRegion(Region region) {
		this.revRegion = region;
	}
	
	public Region getRevRegion() {
		return revRegion;
	}
	
	public HashMap<Integer,Character> getStrands() {
		return strands;
	}
	
	public String getChrom() {
		return chrom;
	}
	
	public LinkedList<SAMRecord> getReads(){
		return list;
	}
	
	public int getTheoreticalNumberOfReads() {
		return (int)Math.round( (double)list.size()/sampleProb);
	}
	
	public HashMap<String,Integer> getIDMap(){
		return idMap;
	}
	
	public int[] getCoverageByBlocks() {
		if(coverageByBlocks == null) {
			this.computeCoverageByBlocks();
		}
		return coverageByBlocks;
	}
	
	private int[] computeCoverageByReads() {
		int[] cov = new int[regionEnd-regionStart+1];
		for(SAMRecord rec: list) {
			for(int i=rec.getAlignmentStart();i<rec.getAlignmentEnd();i++) {
				cov[i-regionStart]++;
			}
		}
		return cov;
	}
	
	private void computeCoverageByBlocks() {
		int[] cov = new int[regionEnd-regionStart+1];
		for(SAMRecord rec: list) {
			List<AlignmentBlock> blocks = rec.getAlignmentBlocks();
			for(AlignmentBlock block : blocks) {
				for(int i=block.getReferenceStart();i<block.getReferenceStart()+block.getLength();i++) {
					cov[i-regionStart]++;
				}
			}
		}
		this.coverageByBlocks = cov;
	}
	
	
	public Region[] splitByCov(int maxLen) {
		int[] cov = computeCoverageByReads();
		
		IntList splitPoints = new IntList();
		splitPoints(splitPoints, 0, cov.length, cov, maxLen);
		
		splitPoints.add(0);
		splitPoints.add(cov.length);
		
		splitPoints.sort();
		
		int[] sar = splitPoints.toArray();
		
		Region[] regions = new Region[splitPoints.length()-1];
		for(int i=0;i<regions.length;i++) {
			regions[i] = new Region(maxcov, sample,strand);
		}
		
		for(SAMRecord rec : list) {
			
			int astart = rec.getAlignmentStart()-regionStart;
			int aend = rec.getAlignmentEnd()-regionStart;
			
			for(int j=1;j<sar.length;j++) {
				if(astart >= sar[j-1] && aend < sar[j]) {
					regions[j-1].addRead(rec);
					break;
				}
			}
			
		}
		
		return regions;
	}
	
	private void splitPoints(IntList splitPoints, int start, int end, int[] cov, int maxLen) {
		int len = end-start+1;		
		
		int idx = ToolBox.getMinIndex(start+(int)(len*0.2), start+(int)(len*0.8), cov);
		splitPoints.add(idx);
		
		if(idx - start > maxLen) {
			splitPoints(splitPoints, start, idx, cov, maxLen);
		}
		if(end - idx > maxLen) {
			splitPoints(splitPoints, idx, end, cov, maxLen);
		}
		
	}
	
	
	public void addRead(SAMRecord read) {
		if(sampleProb == 1.0 || r.nextDouble()<sampleProb) {
			list.add(read);
			if(chrom == null) {
				chrom = read.getReferenceName();
			}else {
				if(!chrom.equals(read.getReferenceName())) {
					throw new RuntimeException();
				}
			}
			if(regionStart == null || read.getAlignmentStart()< regionStart) {
				regionStart = read.getAlignmentStart();
			}
			if(regionEnd == null || read.getAlignmentEnd()> regionEnd) {
				regionEnd = read.getAlignmentEnd();
			}
		}else {
			nOut++;
		}
		
		if(list.size() > maxcov*(regionEnd - regionStart+1)) {
			int before = list.size();
			System.out.print("#"+list.size());
			
			sampleProb *= sample;
			
			list.removeIf(new Predicate<SAMRecord>() {

				@Override
				public boolean test(SAMRecord t) {
					return r.nextDouble() >= sample;
				}
				
			});
			nOut += list.size() - before;
			System.out.println("->"+list.size());
		}
		
	}
	
	public double getScale() {
		return (double)(list.size()+nOut)/(double)list.size();
	}


	public Integer getReferenceIndex() {
		if(list.size() == 0) {
			return null;
		}else {
			return list.getLast().getReferenceIndex();
		}
	}

	public Integer getRegionStart() {
		return regionStart;
	}
	
	public Integer getRegionEnd() {
		return regionEnd;
	}

	public String toString() {
		return chrom+" "+regionStart+"-"+regionEnd+": "+list.size();
	}
	
	public ReadGraph buildGraph(int minIntronLength, int maxGapFilled, Stranded stranded, int maxMM, boolean longReads) {
		
		ReadGraph sg = new ReadGraph(chrom, regionStart, regionEnd, minIntronLength);
		
		Iterator<SAMRecord> it = list.iterator();
		
		this.strands = new HashMap<Integer, Character>();
		this.idMap = new HashMap<String, Integer>();
		
		while(it.hasNext()) {
			SAMRecord sr = it.next();
			
			char strand = '.';
			
			if(stranded == Stranded.FR_FIRST_STRAND) {
				strand = '-';
				if(sr.getReadPairedFlag() && sr.getFirstOfPairFlag() && sr.getReadNegativeStrandFlag()) {
					strand = '+';
				}else if(sr.getReadPairedFlag() && sr.getSecondOfPairFlag() && sr.getMateNegativeStrandFlag()) {
					strand = '+';
				}else if(!sr.getReadPairedFlag()) {
					if(sr.getReadNegativeStrandFlag()) {//TODO was -/+
						strand = '+';
					}else {
						strand = '-';
					}
				}
			}else if(stranded == Stranded.FR_SECOND_STRAND) {
				strand = '+';
				if(sr.getReadPairedFlag() && sr.getFirstOfPairFlag() && sr.getReadNegativeStrandFlag()) {
					strand = '-';
				}else if(sr.getReadPairedFlag() && sr.getSecondOfPairFlag() && sr.getMateNegativeStrandFlag()) {
					strand = '-';
				}else if(!sr.getReadPairedFlag()) {
					if(sr.getReadNegativeStrandFlag()) {//TODO was -/+
						strand = '-';
					}else {
						strand = '+';
					}
				}
			}
			
			int i=-1;
			String name = sr.getReadName();
			if(idMap.containsKey(name)) {
				i = idMap.get(name);
			}else {
				i = idMap.size();
				idMap.put(name, i);
			}
			

			strands.put(i, strand);
			
			sg.addRead(sr,maxMM,longReads);
			
		}
		
		sg.finalizeGaps(maxGapFilled);
		
		return sg;
	}

	public void join(Region revTemp) {
		if(! this.chrom.equals(revTemp.chrom ) ) {
			throw new RuntimeException("chromosomes do not match");
		}
		
		if(revTemp.getReads().size() == 0) {
			return;
		}
		
		if(this.idMap != null || revTemp.idMap != null) {
			throw new RuntimeException("too late");
		}
		
		if(this.sampleProb == revTemp.sampleProb) {
			this.regionStart = Math.min(this.regionStart, revTemp.regionStart);
			this.regionEnd = Math.max(this.regionEnd, revTemp.regionEnd);
			this.list.addAll(revTemp.list);
			this.nOut += revTemp.nOut;
			
		}else if(this.sampleProb <= revTemp.sampleProb) {
			
			double temp = this.sampleProb;
			this.sampleProb /= revTemp.sampleProb;
			for(SAMRecord sr : revTemp.list) {
				this.addRead(sr);
			}
			this.sampleProb = temp;
			
		}else {
			
			double temp = revTemp.sampleProb;
			revTemp.sampleProb /= this.sampleProb;
			for(SAMRecord sr : this.list) {
				revTemp.addRead(sr);
			}
			revTemp.sampleProb = temp;
			
			this.sampleProb = revTemp.sampleProb;
			this.nOut = revTemp.nOut;
			this.regionStart = revTemp.regionStart;
			this.regionEnd = revTemp.regionEnd;
			this.list = revTemp.list;
			
			
		}
		
		
	}
	
	
}
