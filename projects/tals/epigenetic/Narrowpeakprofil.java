package projects.tals.epigenetic;

public class Narrowpeakprofil {

	private boolean[][] peakBool=new boolean[2][];// peakBool[0] = "+" strand and peakBool[1] = "-" strand
	public Narrowpeakprofil(String chrom,int chromLength,boolean[][] peakBool,int before,int after){
		this.peakBool=peakBool;
	}
	
	public boolean isPeakSurroundPos(int Pos,boolean strand) {
		if(strand){//+  true
			return peakBool[0][Pos];
		}else{//-  false
			return peakBool[1][Pos];
		}
	} 
}