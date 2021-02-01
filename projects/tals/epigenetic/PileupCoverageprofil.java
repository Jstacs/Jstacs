package projects.tals.epigenetic;

public class PileupCoverageprofil {

	private String chrom;
	private int chromLength;
	private int[][] numberOfCoveragePositions=new int[2][];// peakBool[0] entspricht "+" Strang und peakBool[1] entspricht "-" Strang
	private int before;
	private int after;
	private boolean calculateAlwaysOnCompleteSeq;
	
	public PileupCoverageprofil(String chrom,int chromLength,int[][] numberOfCoveragePositions,int before,int after, boolean calculateAlwaysOnCompleteSeq){
		this.chrom=chrom;
		this.chromLength=chromLength;
		this.numberOfCoveragePositions=numberOfCoveragePositions;
		this.before=before;
		this.after=after;
		this.calculateAlwaysOnCompleteSeq=calculateAlwaysOnCompleteSeq;
		
	}
	
	public int getnumberOfCoveragePositionsSurroundPos(int Pos,boolean strand) {
		if(strand){//+  true
			return numberOfCoveragePositions[0][Pos];
		}else{//-  false
			return numberOfCoveragePositions[1][Pos];
		}
	} 

}
