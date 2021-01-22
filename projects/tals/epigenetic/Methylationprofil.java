package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.training.SamplingHMMTrainingParameterSet;
import sun.security.util.Length;

public class Methylationprofil {

	private int startPos;
	private String chrom;
	int chromLength;
	private float[] MethylationProb;
	private boolean origStrand;
	private float pseudoCounts;
	private float probCMethylated;
	
	public Methylationprofil(int startPos, String chrom,int chromLength,float[] MethylationProb, boolean strand, float pseudoCounts, float probCMethylated){
		this.startPos = startPos;
		this.chrom=chrom;
		this.chromLength=chromLength;
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
		if(this.origStrand){
			this.startPos = i;
		}else{
			this.startPos = this.chromLength-i;
		}
		
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
	
}
