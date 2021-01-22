package projects.tals.epigenetic;

public class NarrowPeak {
	private String chrom;
	private int startPos;
	private int endPos;
	private float score;
	private float peakValue;
	
	public NarrowPeak(String chrom, int startPos, int endPos, float score,float peakValue){
		this.chrom=chrom;
		this.startPos=startPos;
		this.endPos=endPos;
		this.score=score;
		this.peakValue=peakValue;
	}
	public String getChrom() {
		return this.chrom;
	}
	public int getPeakStartPos() {
		return this.startPos;
	}
	public int getPeakEndPos() {
		return this.endPos;
	}
	public float getPeakScore() {
		return this.score;
	}
	public float getPeakValue() {
		return this.peakValue;
	}
}
