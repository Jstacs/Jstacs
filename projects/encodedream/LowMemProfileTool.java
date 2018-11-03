package projects.encodedream;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.sequenceScores.QuickScanningSequenceScore;
import de.jstacs.utils.IntList;
import de.jstacs.utils.LargeSequenceReader;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;

public class LowMemProfileTool {

	public void run(QuickScanningSequenceScore model2, String genome, PrintStream out2, int bin, String faiFile) throws OperationNotSupportedException, NumberFormatException, IOException {
		
		BufferedReader faidx = new BufferedReader(new InputStreamReader(new FileInputStream(faiFile)));
		String str = null;
		
		HashMap<String,Integer> sizes = new HashMap<>();
		
		while( (str = faidx.readLine()) != null ){
			String[] parts = str.split("\t");
			String chr = parts[0];
			int len = Integer.parseInt(parts[1]);
			sizes.put(chr,len);
		}
		
		faidx.close();
		
		
		double[] scores = new double[2*bin];
		Arrays.fill(scores, Double.NEGATIVE_INFINITY);
		
		
		int currIdx = 0;
		int lastOff = 0;
		String lastId = null;
		
		StringBuffer lastHeader = new StringBuffer();
		BufferedReader read = new BufferedReader(new FileReader(genome));
		
		Pair<IntList,ArrayList<Sequence>> pair = null;
		
		while( (pair = LargeSequenceReader.readNextSequences(read, lastHeader, model2.getLength())) != null ){
			ArrayList<Sequence> seqs = pair.getSecondElement();
			IntList starts = pair.getFirstElement();
			
			for(int i=0;i<seqs.size();i++){
				Sequence seq = seqs.get(i);
				int off = starts.get(i);
				
				String id = seq.getSequenceAnnotationByType("id", 0).getIdentifier();
				
				if(lastId != null && !id.equals(lastId)){
					
					int totLen = sizes.get(lastId);
					for(int j=lastOff;j<totLen;j++){
						scores[2*currIdx] = Double.NEGATIVE_INFINITY;
						scores[2*currIdx+1] = Double.NEGATIVE_INFINITY;
						currIdx++;
						if(currIdx == bin){
							print(lastId,j-bin,scores,out2);
							
							currIdx = 0;
							
							Arrays.fill(scores, Double.NEGATIVE_INFINITY);
						}
					}
					
					currIdx = 0;
					lastOff = 0;
					Arrays.fill(scores, Double.NEGATIVE_INFINITY);
				}
				
				
				
				for(int j=lastOff;j<off;j++){
					/*scores[2*currIdx] = Double.NEGATIVE_INFINITY;
					scores[2*currIdx+1] = Double.NEGATIVE_INFINITY;
					currIdx++;
					if(currIdx == bin){

						print(id,j-bin,scores,out2);
						
						currIdx = 0;
						
						Arrays.fill(scores, Double.NEGATIVE_INFINITY);
						
					}*/
					currIdx++;
					if(currIdx == bin){
						out2.println(id+"\t"+(j-bin+1)+"\tNA\tNA");
						currIdx=0;
					}
				}
				
				
				for(int j=0;j<seq.getLength()-model2.getLength()+1;j++){
					float[] compScore = new float[]{ (float) model2.getLogScoreFor(seq, j), (float) model2.getLogScoreFor(seq.reverseComplement(), seq.getLength()-j-model2.getLength()) };
					
					
					scores[2*currIdx] = compScore[0];
					scores[2*currIdx+1] = compScore[1];
					
					currIdx++;
					if(currIdx == bin){
					
						
						print(id,off+j-bin,scores,out2);
						
						currIdx = 0;
						
						Arrays.fill(scores, Double.NEGATIVE_INFINITY);
						
					}
					lastOff = off+j+1;
				}
				lastId = id;
			}
		}
		
		int totLen = sizes.get(lastId);
		for(int j=lastOff;j<totLen;j++){
			scores[2*currIdx] = Double.NEGATIVE_INFINITY;
			scores[2*currIdx+1] = Double.NEGATIVE_INFINITY;
			currIdx++;
			if(currIdx == bin){
				print(lastId,j-bin,scores,out2);
				
				currIdx = 0;
				
				Arrays.fill(scores, Double.NEGATIVE_INFINITY);
			}
		}
		
	}
	
	
	private void print(String chr, int position, double[] scores, PrintStream out){
		float max = (float) ToolBox.max(scores);
		double ls = Normalisation.getLogSum(scores);
		
		String maxS = Double.isInfinite(max) ? "NA" : max+"";
		String lsS = Double.isInfinite(ls) ? "NA" : ls+"";
		
		out.println(chr+"\t"+(position+1)+"\t"+maxS+"\t"+lsS);
		
		
	}

}
