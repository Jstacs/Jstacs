package projects.dream2016;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.FileManager;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM.InitMethod;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import projects.dimont.AbstractSingleMotifChIPper;

public class RegionModelScan2 {

	
	/**
	 * Arguments:
	 * 0: Slim xml file
	 * 1: path to hg19.genome.fa file
	 * 
	 * output will be generated as args[0]_winscores2.txt
	 */
	public static void main(String[] args) throws Exception {
		
		GenDisMixClassifier cl = new GenDisMixClassifier(FileManager.readFile(args[0]));
		
		DifferentiableStatisticalModel model = ((AbstractSingleMotifChIPper) cl.getDifferentiableSequenceScore(0)).getFunction(0);
		
		StrandDiffSM model2 = new StrandDiffSM(model, 1, true, InitMethod.INIT_FORWARD_STRAND, 0.5);
		
		
		
		double t1 = model2.getLength()*Math.log(0.25);
		double t2 = model2.getLength()*Math.log(1.0/3.0);
		
		
		Pair<int[][], DataSet> seqs = null;
		
		StringBuffer lastHeader = new StringBuffer();
		BufferedReader read = new BufferedReader(new FileReader(args[1]));
		
		DecimalFormat nf = new DecimalFormat("#.####");
		DecimalFormatSymbols syms = nf.getDecimalFormatSymbols();
		syms.setInfinity("Inf");
		syms.setNaN("NaN");
		nf.setDecimalFormatSymbols(syms);
		
		PrintWriter wr = new PrintWriter(args[0]+"_winscores2.txt");
		
		double[] scores1 = new double[50];
		Arrays.fill(scores1, Double.NEGATIVE_INFINITY);
		
		byte[] strand1 = new byte[50];
		
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
					Arrays.fill(scores1, Double.NEGATIVE_INFINITY);
					wr.println("["+id+"]");
				}else{
					//wr.println("new part");
				}
				
				int off = idxs[0][i];
				//wr.println("["+id+" "+off+" "+lastOff+" "+currIdx+"]");
				
				for(int j=lastOff;j<off;j++){
					scores1[currIdx] = Double.NEGATIVE_INFINITY;
					strand1[currIdx] = 0;
					currIdx++;
					if(currIdx == scores1.length){

						print(wr,nf,id,scores1,strand1,t1,t2,lastOff,j,currIdx);
						
						currIdx = 0;
						
						Arrays.fill(scores1, Double.NEGATIVE_INFINITY);
						Arrays.fill(strand1, (byte)0);
						
					}
				}
				
				
				for(int j=0;j<seq.getLength()-model2.getLength()+1;j++){
					double[] compScore = model2.getComponentScores(seq, j);
					
					double score = Normalisation.getLogSum(compScore);
					byte strand = compScore[0] == compScore[1] ? 0 : (compScore[0] > compScore[1] ? (byte)1 : (byte)-1);
					
					scores1[currIdx] = score;
					strand1[currIdx] = strand;
					
					currIdx++;
					if(currIdx == scores1.length){
					
						
						print(wr,nf,id,scores1,strand1,t1,t2,off,j,currIdx);
						
						currIdx = 0;
						
						Arrays.fill(scores1, Double.NEGATIVE_INFINITY);
						Arrays.fill(strand1, (byte)0);
						
					}
					lastOff = off+j+1;
				}
				lastId = id;
			}
		}
		
		read.close();
		wr.close();
		
		
	}
	
	private static void print(PrintWriter wr, NumberFormat nf, String id, double[] scores1, byte[] strands, double t1, double t2, int off, int j, int end){
		double mean, max;
		int nabove1, nabove2;
		int idx = -1;
		String strand;
		double ns = 0;
		if(scores1 == null){
			mean = max = Double.NEGATIVE_INFINITY;
			strand = "..";
			nabove1 = 0;
			nabove2 = 0;
		}else{
			mean = Normalisation.getLogSum(0,end,scores1);
			ns = 0;
			for(int i=0;i<end;i++){
				ns += strands[i]*Math.exp(scores1[i]-mean);
			}
			idx = ToolBox.getMaxIndex(0,end,scores1);
			max = scores1[idx];//ToolBox.max(0,end,scores1);
			strand = (strands[idx] == 0 ? "." : (strands[idx]>0 ? "+" : "-"))+(ns==0 || Double.isNaN(ns) ? "." : (ns>0 ? "+" : "-"));
			nabove1 = aboveThreshold(scores1,t1,end);
			nabove2 = aboveThreshold(scores1,t2,end);
		}
		
		
		//wr.println(id+"\t"+(off+j-50+1)+"\t"+(off+j+1)+"\t"+nf.format(mean)+"\t"+nf.format(max)+"\t"+nabove1+"\t"+nabove2+"\t"+strand+"\t"+idx+"\t"+ns);
		//wr.println(nf.format(-mean)+"\t"+nf.format(-max)+"\t"+nabove1+"\t"+nabove2+"\t"+strand);
		
	    wr.println(nf.format(-mean)+"\t"+nf.format(-max)+"\t"+nabove1+"\t"+nabove2+"\t"+strand+idx);
		
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
