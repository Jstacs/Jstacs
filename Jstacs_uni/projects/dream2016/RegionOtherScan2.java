package projects.dream2016;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.DiscreteSequenceEnumerator;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;

public class RegionOtherScan2 {

	
	public static void main(String[] args) throws Exception {
		
		Pair<int[][], DataSet> seqs = null;
		
		StringBuffer lastHeader = new StringBuffer();
		BufferedReader read = new BufferedReader(new FileReader(args[0]));
		
		DecimalFormat nf = new DecimalFormat("#.####");
		DecimalFormatSymbols syms = nf.getDecimalFormatSymbols();
		syms.setInfinity("Inf");
		syms.setNaN("NaN");
		nf.setDecimalFormatSymbols(syms);
		
		PrintWriter wr = new PrintWriter(args[0]+"_other.txt");
		
		double[] score = null;
		int sl = 50;
		
		int currIdx = 0;
		int lastOff = 0;
		String lastId = null;
		
		Pair<HashMap<Sequence,Integer>,double[]> temp = getIndex(1);
		HashMap<Sequence,Integer> i1 = temp.getFirstElement();
		double[] c1 = temp.getSecondElement();
		double[] e1 = c1.clone();
		
		temp = getIndex(2);
		HashMap<Sequence,Integer> i2 = temp.getFirstElement();
		double[] c2 = temp.getSecondElement();
		double[] e2 = c2.clone();
		
		temp = getIndex(3);
		HashMap<Sequence,Integer> i3 = temp.getFirstElement();
		double[] c3 = temp.getSecondElement();
		double[] e3 = c3.clone();
		
		
		
		while( (seqs = readNextSequences(read, lastHeader)) != null ){
			DataSet ds = seqs.getSecondElement();
			int[][] idxs = seqs.getFirstElement();
			for(int i=0;i<ds.getNumberOfElements();i++){
				Sequence seq = ds.getElementAt(i);
				String id = seq.getSequenceAnnotationByType("id", 0).getIdentifier();
				
				if(!id.equals(lastId)){
					currIdx = 0;
					lastOff = 0;
					Arrays.fill(c1, 0);Arrays.fill(c2, 0);Arrays.fill(c3, 0);
					wr.println("["+id+"]");
				}else{
					//wr.println("new part");
				}
				
				int off = idxs[0][i];
				//wr.println("["+id+" "+off+" "+lastOff+" "+currIdx+"]");
				
				for(int j=lastOff;j<off;j++){
					
					currIdx++;
					if(currIdx == sl){

						score = getFeatures(null, j,i1,c1,i2,c2,i3,c3,e1,e2,e3);
						
						print(wr,id,score, lastOff,j,currIdx);
						
						Arrays.fill(c1, 0);Arrays.fill(c2, 0);Arrays.fill(c3, 0);
						
						currIdx = 0;
						
					}
				}
				
				
				for(int j=0;j<seq.getLength();j++){
					
					addCounts(seq, j, i1, c1);
					addCounts(seq, j, i2, c2);
					addCounts(seq, j, i3, c3);
					
					currIdx++;
					if(currIdx == sl){
					
						score = getFeatures(seq,j,i1,c1,i2,c2,i3,c3,e1,e2,e3);
						
						print(wr,id,score,off,j,currIdx);
						
						Arrays.fill(c1, 0);Arrays.fill(c2, 0);Arrays.fill(c3, 0);
						
						currIdx = 0;
						
						
					}
					lastOff = off+j+1;
				}
				lastId = id;
			}
		}
		
		read.close();
		wr.close();
		
		
	}
	
	private static Pair<HashMap<Sequence,Integer>,double[]> getIndex(int k) throws OperationNotSupportedException{
		DiscreteSequenceEnumerator en1 = new DiscreteSequenceEnumerator(DNAAlphabetContainer.SINGLETON, k, true);
		HashMap<Sequence, Integer> i1 = new HashMap<Sequence, Integer>();
		int idx = 0;
		while(en1.hasMoreElements()){
			Sequence s = en1.nextElement();
			i1.put(s, idx);
			i1.put(s.reverseComplement(), idx);
			idx++;
		}
		double[] c1 = new double[idx];
		en1.reset();
		while(en1.hasMoreElements()){
			Sequence s = en1.nextElement();
			c1[i1.get(s)]++;
			c1[i1.get(s.reverseComplement())]++;
		}
		Normalisation.sumNormalisation(c1);
		return new Pair<HashMap<Sequence,Integer>, double[]>(i1, c1);
	}
	
	
	private static double[] getFeatures(Sequence seq, int j, HashMap<Sequence, Integer> i1, double[] c1,
			HashMap<Sequence, Integer> i2, double[] c2,
			HashMap<Sequence, Integer> i3, double[] c3,
			double[] e1, double[] e2, double[] e3) throws IllegalArgumentException, WrongAlphabetException {
		if(seq == null){
			return new double[]{0.5,2.0/16.0,0,0,0};
		}
		
		double gc = c1[i1.get(Sequence.create(DNAAlphabetContainer.SINGLETON, "C"))];
		double at = c1[i1.get(Sequence.create(DNAAlphabetContainer.SINGLETON, "A"))];
		gc /= (gc+at);
		if(gc < 0 || gc > 1){
			System.out.println(gc);
			System.out.println(Arrays.toString(c1));
			System.exit(1);
		}
		
		double cpg = c2[i2.get(Sequence.create(DNAAlphabetContainer.SINGLETON, "CG"))];
		//System.out.println(cpg+" "+ToolBox.sum(c2)+" "+seq.toString(Math.max(0, j-50), j));
		cpg /= ToolBox.sum(c2);
		
		
		double en1 = getKL(c1,e1);
		double en2 = getKL(c2,e2);
		double en3 = getKL(c3,e3);
		return new double[]{gc,cpg,en1,en2,en3};
	}
	
	private static double getKL(double[] counts, double[] ref){
		double n = ToolBox.sum(counts);
		
		double en = 0;
		for(int j=0;j<counts.length;j++){
			if(counts[j] > 0){
				en += counts[j]/n * (Math.log(counts[j]/n) - Math.log(ref[j]));
			}
		}
		
		return en;
	}
	
	private static void addCounts(Sequence seq, int j, HashMap<Sequence, Integer> indexes, double[] counts){
		int sl = indexes.keySet().iterator().next().getLength();
		if(seq.getLength()-j>=sl){
			Sequence sub = seq.getSubSequence(j, sl);
			if(j==0 && sl == seq.getLength()){
				sub = sub.annotate(false,null);
			}
			
			try{
				counts[ indexes.get(sub) ]++;
			}catch(NullPointerException ex){
				System.out.println(j+" "+sl+" "+seq.getLength());
				System.out.println(seq.getSubSequence(j, sl));
				System.out.println(indexes);
				System.out.println(indexes.get(seq.getSubSequence(j, sl)));
				Sequence temp = seq.getSubSequence(j, sl);
				Sequence temp2 = indexes.keySet().toArray(new Sequence[0])[3];
				System.out.println(temp+" "+temp2+" "+indexes.get(temp)+" "+indexes.get(temp2)+" "+temp.hashCode()+" "+temp2.hashCode()+" "+seq.hashCode());
				System.out.println(sub+" "+temp2+" "+indexes.get(sub)+" "+indexes.get(temp2)+" "+sub.hashCode()+" "+temp2.hashCode()+" "+seq.hashCode());
				throw ex;
			}
		}
	}
	
	
	private static void print(PrintWriter wr, String id, double[] scores1, int off, int j, int end){
		double en1, en2, en3, gc, cpg;
		

		en1 = scores1[2];
		en2 = scores1[3];
		en3 = scores1[4];
		gc = scores1[0];
		cpg = scores1[1];
		
		//wr.println(id+"\t"+(off+j-50+1)+"\t"+(off+j+1)+"\t"+nf.format(mean)+"\t"+nf.format(max)+"\t"+nabove1+"\t"+nabove2+"\t"+strand+"\t"+idx+"\t"+ns);
		//wr.println(nf.format(-mean)+"\t"+nf.format(-max)+"\t"+nabove1+"\t"+nabove2+"\t"+strand);
		
	    wr.println(gc+"\t"+cpg+"\t"+en1+"\t"+en2+"\t"+en3);
		
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
