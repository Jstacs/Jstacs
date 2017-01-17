package projects.dream2016;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;

public class RegionTractScan2 {

	
	public static void main(String[] args) throws Exception {
		
		Pair<int[][], DataSet> seqs = null;
		
		StringBuffer lastHeader = new StringBuffer();
		BufferedReader read = new BufferedReader(new FileReader(args[0]));
		
		DecimalFormat nf = new DecimalFormat("#.####");
		DecimalFormatSymbols syms = nf.getDecimalFormatSymbols();
		syms.setInfinity("Inf");
		syms.setNaN("NaN");
		nf.setDecimalFormatSymbols(syms);
		
		PrintWriter wr = new PrintWriter(args[0]+"_tracts.txt");
		
		double[] scoresSing = new double[50];
		Arrays.fill(scoresSing, 0);
		double[] scoresBoth = new double[50];
		Arrays.fill(scoresBoth, 0);
		
		
		int currIdx = 0;
		int lastOff = 0;
		String lastId = null;
		
		while( (seqs = readNextSequences(read, lastHeader)) != null ){
			DataSet ds = seqs.getSecondElement();
			int[][] idxs = seqs.getFirstElement();
			for(int i=0;i<ds.getNumberOfElements();i++){
				Sequence seq = ds.getElementAt(i);
				String id = seq.getSequenceAnnotationByType("id", 0).getIdentifier();
				
				if(!id.equals(lastId)){
					currIdx = 0;
					lastOff = 0;
					Arrays.fill(scoresSing, 0);
					Arrays.fill(scoresBoth, 0);
					wr.println("["+id+"]");
				}else{
					//wr.println("new part");
				}
				
				int off = idxs[0][i];
				//wr.println("["+id+" "+off+" "+lastOff+" "+currIdx+"]");
				
				for(int j=lastOff;j<off;j++){
					scoresSing[currIdx] = 0;
					scoresBoth[currIdx] = 0;
					currIdx++;
					if(currIdx == scoresSing.length){

						print(wr,id,scoresSing,scoresBoth, lastOff,j,currIdx);
						
						currIdx = 0;
						
						Arrays.fill(scoresSing, 0);
						Arrays.fill(scoresBoth, 0);
						
					}
				}
				
				
				for(int j=0;j<seq.getLength();j++){
					
					scoresSing[currIdx] = getTractLength(seq,j);
					scoresBoth[currIdx] = getGCATTractLength(seq,j);
					
					currIdx++;
					if(currIdx == scoresSing.length){
					
						
						print(wr,id,scoresSing,scoresBoth,off,j,currIdx);
						
						currIdx = 0;
						
						Arrays.fill(scoresSing, 0);
						Arrays.fill(scoresBoth, 0);
						
					}
					lastOff = off+j+1;
				}
				lastId = id;
			}
		}
		
		read.close();
		wr.close();
		
		
	}
	
	private static int getTractLength(Sequence seq, int j) {
		int len = 0;
		int curr = seq.discreteVal(j);
		int i=j;
		while( i>= 0 && seq.discreteVal(i)==curr ){
			len++;
			i--;
		}
		i=j+1;
		while( i<seq.getLength() && seq.discreteVal(i) == curr ){
			len++;
			i++;
		}
		if(curr == 0 || curr == 3){
			return len;
		}else{
			return -len;
		}
	}
	
	private static int getGCATTractLength(Sequence seq, int j){
		int len = 0;
		int curr = seq.discreteVal(j);
		int i=j;
		while( i>= 0 && (seq.discreteVal(i)==curr||seq.discreteVal(i)==3-curr) ){
			len++;
			i--;
		}
		i=j+1;
		while( i<seq.getLength() && (seq.discreteVal(i)==curr||seq.discreteVal(i)==3-curr) ){
			len++;
			i++;
		}
		if(curr == 0 || curr == 3){
			return len;
		}else{
			return -len;
		}
	}
	

	private static void print(PrintWriter wr, String id, double[] scoresSing, double[] scoresBoth, int off, int j, int end){
		int min, max, min1, max1;
		if(scoresSing == null){
			min = max = min1 = max1 = 0;
		}else{
			
			max = (int)ToolBox.max(scoresSing);
			min = (int)ToolBox.min(scoresSing);
			max1 = (int)ToolBox.max(scoresBoth);
			min1 = (int)ToolBox.min(scoresBoth);
			//System.out.println(Arrays.toString(scores1));
		}
		
		if(max < 0){
			max = 0;
		}
		
		if(min > 0){
			min = 0;
		}
		
		if(max1 < 0){
			max1 = 0;
		}
		
		if(min1 > 0){
			min1 = 0;
		}
		
	
		
		//wr.println(id+"\t"+(off+j-50+1)+"\t"+(off+j+1)+"\t"+nf.format(mean)+"\t"+nf.format(max)+"\t"+nabove1+"\t"+nabove2+"\t"+strand+"\t"+idx+"\t"+ns);
		//wr.println(nf.format(-mean)+"\t"+nf.format(-max)+"\t"+nabove1+"\t"+nabove2+"\t"+strand);
		
	    wr.println((max)+"\t"+(-min)+"\t"+(max1)+"\t"+(-min1));
		
	}

	private static int aboveThreshold(double[] scores1, double t, int end) {
		int n = 0;
		for(int i=0;i<end;i++){
			if(scores1[i] > t){
				n++;
			}
		}
		return n;
	}

	public static double getConsensusScore(double[][] pwm){
		double score = 0.0;
		for(int i=0;i<pwm.length;i++){
			score += ToolBox.max(pwm[i]);
		}
		return score;
	}
	
	public static Pair<int[][],DataSet> readNextSequences(BufferedReader read, StringBuffer lastHeader) throws Exception {
		//System.out.println("started reading");
		String str = null;
		
		StringBuffer line = new StringBuffer();
		
		IntList starts = new IntList();
	//	IntList ends = new IntList();
		
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();

		Pattern acgt = Pattern.compile( "[ACGT]+", Pattern.CASE_INSENSITIVE );
		
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		
		int size = 0;
		
		while( (str = read.readLine()) != null || line.length() > 0 ){
			if(str != null){
				str = str.trim();
			}
			if(str == null || str.startsWith( ">" )){//next sequence
				String header = lastHeader.toString();
				if(str != null){
					lastHeader.delete( 0, lastHeader.length() );
					lastHeader.append( str.substring( 1 ).trim() );
				}
				if(line.length() > 0){//we have a sequence
					String seqStr = line.toString();
					line.delete( 0, line.length() );
					Matcher match = acgt.matcher( seqStr );
					while(match.find()){
						int start = match.start();
						int end = match.end();
						
						SequenceAnnotation annotation = new SequenceAnnotation( "id",header);
						Sequence seq = Sequence.create( DNAAlphabetContainer.SINGLETON, seqStr.substring( start, end ) );
						seq = seq.annotate( false, annotation );
						seqs.add( seq );
						size += end-start;
						starts.add( start );
					//	ends.add( seqStr.length()-end );
					}
					if(size > 1E7 || str == null){
						int s = 0;
						return new Pair<int[][],DataSet>(new int[][]{starts.toArray()/*,ends.toArray()*/},new DataSet( "", seqs ));
					}
				}
			}else{
				line.append( str );
			}	
		}
		return null;
	}
	
}
