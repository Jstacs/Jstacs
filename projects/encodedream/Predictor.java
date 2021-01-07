package projects.encodedream;
/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.zip.GZIPOutputStream;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;

public class Predictor {

	private GenDisMixClassifier[] cls;
	private FeatureReader reader;
	private int binsBefore;
	private int binsAfter;
	
	public Predictor(GenDisMixClassifier[] cls, FeatureReader reader, int binsBefore, int binsAfter){
		this.cls = cls;
		this.reader = reader;
		this.binsBefore = binsBefore;
		this.binsAfter = binsAfter;
	}
	
	
	public File predict(HashMap<String,Integer> sizes,LinkedList<String> chroms) throws Exception{
		
		if(chroms == null){
			chroms = new LinkedList<>(sizes.keySet());
			Collections.sort(chroms);
		}
		
		File temp = File.createTempFile("preds", ".tsv.gz");
		temp.deleteOnExit();
		PrintWriter wr = new PrintWriter(new GZIPOutputStream(new FileOutputStream(temp)));
		
		for(String chr : chroms){
			double[] preds = predict(chr,sizes.get(chr));
			reader.reset();
			reader.findChr(chr);
			int i=0;
			do{
				String chr2 = reader.getCurrentChromosome();
				int start = reader.getCurrentStart();
				if(!chr.equals(chr2)){
					break;
				}
				wr.println(chr+"\t"+start+"\t"+preds[i]);
				i++;
			}while(reader.readNextFeatureVector());
		}
		
		wr.close();
		return temp;
	}
	
	
	public double[] predict(String chr, int size) throws Exception{
		reader.reset();
		boolean found = reader.findChr(chr);
		if(!found) {
			throw new RuntimeException("Did not find chromosome "+chr+" in feature files.");
		}
		
//PrintWriter wr = new PrintWriter("/Users/dev/Downloads/pred+"+chr+".tsv");
		
		
		double[][] scores = new double[cls.length][size];
		
		int i=0;
		int lastStart = 0;
		do{
			if(i<size){
				Sequence seq = reader.getCurrentSequence();
				
				for(int j=0;j<cls.length;j++){
					scores[j][i] = cls[j].getScore(seq, 0) - cls[j].getScore(seq, 1);
				}
/*wr.print(reader.getCurrentChromosome()+" "+reader.getCurrentStart()+" "+i+" "+reader.getCurrentLabel());
for(int j=0;j<cls.length;j++){
	wr.print(" "+scores[j][i]);
}
wr.println();*/
				
				i++;
			}else{
				break;
			}
		}while(reader.readNextFeatureVector());

//wr.close();
		
		double[] pred = aggregate(scores);
		
		
/*wr = new PrintWriter("/Users/dev/Downloads/agg+"+chr+".tsv");
for(int j=0;j<pred.length;j++){
	wr.println(j+" "+pred[j]);
}
wr.close();*/
		
		return pred;
	}
	
	
	/*public File predict(String[] chrs) throws Exception{
		HashMap<String,File> used = new HashMap<String,File>();
		for(int i=0;i<chrs.length;i++){
			File file = File.createTempFile(chrs[i], ".itt");
			file.deleteOnExit();
			used.put(chrs[i],file);
		}
		
		
		HashMap<String,PrintWriter> writers = new HashMap<>();
		for(String key : used.keySet()){
			writers.put(key, new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(used.get(key))))));
		}
		
		reader.reset();
		while(reader.readNextFeatureVector()){
			String chr = reader.getCurrentChromosome();
			PrintWriter wr = writers.get(chr);
			if(used.containsKey(chr)){
				Sequence seq = reader.getCurrentSequence();
				wr.print(reader.getCurrentLabel());
				for(GenDisMixClassifier curr : cls){
					double sc = curr.getScore(seq, 0) - curr.getScore(seq, 1);
					wr.print("\t"+sc);
				}
				wr.println();
			}
		}
		for(String key : used.keySet()){
			writers.get(key).close();
		}
		
		
		LinkedList<double[]> buffer = new LinkedList<>();
		
		for(String chr : used.keySet()){
			BufferedReader read = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(used.get(chr)))));
			int i=0;
			String str = null;
			while( (str = read.readLine()) != null ){
				String[] parts = str.split("\t");
				double[] curr = new double[parts.length-1];
				for(int i=0;i<curr.)
				
				i++;
			}
		}
		
		
	}*/
	
	
	public double[] aggregate(double[][] scs){
		
		double[] preds = new double[scs[0].length];
		
		for(int j=0;j<scs[0].length;j++){
			double all = 0.0;
			for(int i=0;i<cls.length;i++){
				int start = Math.max(0, j-binsBefore);
				int end = Math.min(scs[i].length, j+binsAfter+1);
				
				double sum = 0.0;
				for(int k=start;k<end;k++){
					sum += Math.log1p(-1.0/(1.0+Math.exp(-scs[i][k])));
				}
				if(Double.isNaN(sum)){
					sum = 0.0;
				}
				all += 1.0 - Math.exp(sum);
				//System.out.println(binsBefore+" "+binsAfter+" "+start+" "+end+" "+j+" "+i+" "+(1.0 - Math.exp(sum))+" "+all);
			}
			all /= cls.length;
			preds[j] = all;
		}
		
		return preds;
	}
	
	
	public double[] predict(DataSet data) throws Exception{
		double[][] scs = new double[cls.length][];
		
		for(int i=0;i<cls.length;i++){
			scs[i] = cls[i].getScores(data);
		}
		
		return aggregate(scs);
	}
	
	
}
