package projects.tals.epigenetic;

public class Methylationprofil {

	private int startPos;
	private String chrom;
	int seqLength;
	private float[] MethylationProb;
	private boolean origStrand;
	private float pseudoCounts;
	private float probCMethylated;
	
	public Methylationprofil(int startPos, String chrom,int chromLength,float[] MethylationProb, boolean strand, float pseudoCounts, float probCMethylated){
		this.startPos = startPos;
		this.chrom=chrom;
		this.seqLength=chromLength;
		this.MethylationProb=MethylationProb;
		this.origStrand=strand;
		this.pseudoCounts=pseudoCounts;
		this.probCMethylated=probCMethylated;
	}
	
	public double getMethylPropAtPos() {
		return MethylationProb[this.startPos];
	} 
	
	public double getMethylPropAtPos(int offset) {
		if(this.origStrand==false){
			return MethylationProb[this.startPos-offset-1];
		}else{
			return MethylationProb[this.startPos+offset];
		}
		
	} 
	
	public void setStartPos ( int i ) {
		
			this.startPos = i;
		
		
	}
	
	public int getStartPos () {
		return this.startPos;
	}
	
	public boolean getStrand(){
		return this.origStrand;
	}
	
	public void setStrand(boolean strand){
		this.origStrand=strand;
	}
	
	public double getPseudoProb(){
		return this.pseudoCounts*this.probCMethylated;
	}
	public void setSeqLength(int seqLength){
		this.seqLength=seqLength;
	}
}
