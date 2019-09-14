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
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Random;
import java.util.zip.GZIPInputStream;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.sequences.ArbitrarySequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;

public class FeatureReader {

	private int numBins;
	private String labelsFile;
	private String dnaseFile;
	private String[] motifFiles;
	
	private BufferedReader labelsReader;
	private BufferedReader dnaseReader;
	private BufferedReader[] motifReaders;
	
	private LinkedList<Character> currLabels;
	private LinkedList<double[][]> currFeatures;
	private LinkedList<String[]> currLines;
	
	private Integer sequenceLength;
	private Integer motifsLength;
	private Integer dnaseLength;
	
	private AlphabetContainer alphabetContainer = new AlphabetContainer(new ContinuousAlphabet(true));
	
	
	public static HashMap<String,Integer> getSizes(String faiFile, int bin) throws NumberFormatException, IOException{
		BufferedReader faidx = new BufferedReader(new InputStreamReader(new FileInputStream(faiFile)));
		String str = null;
		
		HashMap<String,Integer> sizes = new HashMap<>();
		
		while( (str = faidx.readLine()) != null ){
			String[] parts = str.split("\t");
			String chr = parts[0];
			int len = Integer.parseInt(parts[1]);
			sizes.put(chr,len/bin);
		}
		
		
		return sizes;
	}
	
	
	public FeatureReader(int numBins, String labelsFile, String dnaseFile, String... motifFiles){
		this.numBins = numBins;
		this.labelsFile = labelsFile;
		this.dnaseFile = dnaseFile;
		this.motifFiles = motifFiles;
		this.motifReaders = new BufferedReader[motifFiles.length];
	}
	
	public void reset() throws FileNotFoundException, IOException{
		close();
		if(labelsFile != null){
			this.labelsReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(labelsFile))));
		}else{
			this.labelsReader = null;
		}

		this.dnaseReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(dnaseFile))));
		for(int i=0;i<motifFiles.length;i++){
			motifReaders[i] = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(motifFiles[i])))); 
		}
		currLabels = new LinkedList<>();
		currFeatures = new LinkedList<>();
		currLines = new LinkedList<>();
	}
	
	public void close() throws IOException{
		if(labelsReader != null){
			labelsReader.close();
		}
		if(dnaseReader != null){
			dnaseReader.close();
		}
		for(int i=0;i<motifReaders.length;i++){
			if(motifReaders[i] != null){
				motifReaders[i].close();
			}
		}
	}
	
	public String[] readNextLines() throws IOException{
		String[] temp = new String[2+motifReaders.length];
		if(labelsReader != null){
			temp[0] = labelsReader.readLine();
		}else{
			temp[0] = "";
		}
		temp[1] = dnaseReader.readLine();
		for(int i=0;i<motifReaders.length;i++){
			temp[2+i] = motifReaders[i].readLine();
		}
		if(!check(temp)){//TODO efficiency
			System.err.println("Chromosomes do not match between input files.");
			System.err.println(Arrays.toString(temp));
			System.exit(1);
		}
		return temp;
	}
	
	public double[][] getFeatureValues(String[] lines) throws NumberFormatException {
		int num = 0;
		int numDnase = 0;
		int numMotifs = 0;
		double[][] feats = new double[lines.length-1][];
		for(int i=1;i<lines.length;i++){
			String[] parts = lines[i].split("\t");
			feats[i-1] = new double[parts.length-2];
			if(i-1==0){
				numDnase += feats[i-1].length;
			}else{
				numMotifs += feats[i-1].length;
			}
			num += feats[i-1].length;
			for(int j=2;j<parts.length;j++){
				if(parts[j].charAt(0)=='N'){
					feats[i-1][j-2] = Double.NaN;
				}else{
					feats[i-1][j-2] = Double.parseDouble(parts[j]);
				}
			}
		}
		if(sequenceLength == null){
			sequenceLength = num*numBins;
			motifsLength = numMotifs*numBins;
			dnaseLength = numDnase*numBins;
		}else{
			if(sequenceLength != num*numBins){
				throw new NumberFormatException("Problem on input line "+Arrays.toString(lines));
			}
		}
		return feats;
	}
	
	public boolean check(String[] lines){
		if(lines[0] == null){
			return true;
		}
		if (lines[1] == null) {
            return false;
        }
		String[] parts = lines[1].split("\t");
		String ref = parts[0]+"\t"+parts[1]+"\t";
		boolean check = true;
		int i = (lines[0].length()>0 ? 0 : 1);
		while(check && i < lines.length){
			check &= lines[i] != null && lines[i].startsWith(ref);
			i++;
		}
		return check;
	}
	
	public double getCurrentDNaseMedian(){
		return currFeatures.get((numBins+1)/2)[0][1];
	}
	
	public double getCurrentDNaseMin(){
		return currFeatures.get((numBins+1)/2)[0][0];
	}
	
	public double getCurrentMotifMax(int numMotif){
		return currFeatures.get((numBins+1)/2)[numMotif+1][0];
	}
	
	public double getDNaseMedian(String[] lines){
		return Double.parseDouble(lines[1].split("\t")[3]);
	}
	
	public Character getLabel(String[] lines){
		if(lines[0].length() == 0){
			return null;
			//throw new RuntimeException("No labels");
		}
		return lines[0].split("\t")[2].charAt(0);
	}
	
	
	public double[][] getPositiveHistogram() throws FileNotFoundException, IOException{
		reset();
		DoubleList values = new DoubleList();
		String[] temp = null;
		while( (temp = readNextLines())[0] != null ){
			char lab = getLabel(temp);
			if(lab == 'S'){
				values.add(getDNaseMedian(temp));
			}
		}
		
		if(values.length()==0){
			throw new RuntimeException("Not a single bin with label 'S'. Please check input.");
		}
		
		double min = values.min(0, values.length());
		double max = values.max(0, values.length());
		double[] breaks = new double[20];//number of histogram bins
		double step = (max-min)/breaks.length;
		breaks[0] = min+step;
		for(int i=1;i<breaks.length;i++){
			breaks[i] = breaks[i-1]+step;
		}
		
		double[] counts = new double[breaks.length];
		for(int i=0;i<values.length();i++){
			double curr = values.get(i);
			int idx = Arrays.binarySearch(breaks, curr);
			if(idx < 0){
				idx = -idx-1;
			}
			if(idx >= breaks.length){//numerical accuracy
				idx = breaks.length-1;
			}
			counts[idx]++;
		}
		return new double[][]{breaks,counts};
	}
	
	public double[] getNegativeHistogram(double[] breaks) throws FileNotFoundException, IOException{
		double[] counts = new double[breaks.length];
		
		reset();
		String[] temp = null;
		while( (temp = readNextLines())[0] != null ){
			char lab = getLabel(temp);
			if(lab == 'U'){
				double curr = getDNaseMedian(temp);
				int idx = Arrays.binarySearch(breaks, curr);
				if(idx < 0){
					idx = -idx-1;
				}
				if(idx >= breaks.length){//fixedBins
					idx = breaks.length-1;
				}
				counts[idx]++;
				
			}
		}
		
		return counts;
		
	}
	
	public double[][] getSamplingProbsAndWeights(double[] posHist, double[] negHist){
		double[] probs = new double[posHist.length];
		double[] weights = new double[posHist.length];
		for(int i=0;i<posHist.length;i++){
			probs[i] = Math.min(posHist[i]*4.0/negHist[i],1.0);
			weights[i] = probs[0]/probs[i];
		}
		return new double[][]{probs,weights};
	}
	
	
	
	public char getCurrentLabel(){
		return currLabels.get((numBins+1)/2);
	}
	
	
	public Sequence getCurrentSequence() throws WrongAlphabetException, WrongSequenceTypeException{
		double[] temp = new double[sequenceLength];
		int k=0;
		for(int i=0;i<currFeatures.size();i++){
			double[][] feat = currFeatures.get(i);
			for(int j=0;j<feat.length;j++){
				for(int l=0;l<feat[j].length;l++,k++){
					temp[k] = feat[j][l];
				}
			}
		}
		return new ArbitrarySequence(alphabetContainer, temp);
	}
	
	public Sequence getCurrentDNaseSequence() throws WrongAlphabetException, WrongSequenceTypeException{
		double[] temp = new double[dnaseLength];
		int k=0;
		for(int i=0;i<currFeatures.size();i++){
			double[][] feat = currFeatures.get(i);
			for(int l=0;l<feat[0].length;l++,k++){
					temp[k] = feat[0][l];
			}
		}
		return new ArbitrarySequence(alphabetContainer, temp);
	}
	
	public Sequence getCurrentMotifsSequence() throws WrongAlphabetException, WrongSequenceTypeException{
		double[] temp = new double[motifsLength];
		int k=0;
		for(int i=0;i<currFeatures.size();i++){
			double[][] feat = currFeatures.get(i);
			for(int j=1;j<feat.length;j++){
				for(int l=0;l<feat[j].length;l++,k++){
					temp[k] = feat[j][l];
				}
			}
		}
		return new ArbitrarySequence(alphabetContainer, temp);
	}
	
	
	public boolean readNextFeatureVector() throws IOException{
		
		while(currLabels.size()<(numBins-1)/2-1){
			String[] temp = readNextLines();
			if(temp[0] == null){
				return false;
			}
			currLabels.add(getLabel(temp));
			currFeatures.add(getFeatureValues(temp));
			currLines.add(temp);
			
		}
		
		Character lab = null;
		double[][] vals = null;
		String[] temp = null;
		while(vals==null){
			temp = readNextLines();
			if(temp[0] == null){
				return false;
			}
			lab = getLabel(temp);
			vals = getFeatureValues(temp);
		}
		
		
		currLabels.add(lab);
		currFeatures.add(vals);
		currLines.add(temp);
		
		while(currLabels.size()<numBins){
			currLabels.addFirst(currLabels.getFirst());
			currFeatures.addFirst(currFeatures.getFirst());
			currLines.addFirst(currLines.getFirst());
		}
		
		if(currLabels.size()>numBins){
			currLabels.removeFirst();
			currFeatures.removeFirst();
			currLines.removeFirst();
		}
		return true;
	}
	
	
	public boolean readNextFeatureVector2() throws IOException{
		while(currLabels.size()<(numBins-1)/2){
			String[] temp = readNextLines();
			if(temp[0] == null){
				return false;
			}
			currLabels.add(getLabel(temp));
			currFeatures.add(getFeatureValues(temp));
			currLines.add(temp);
			
		}
		
		Character lab = null;
		double[][] vals = null;
		String[] temp = null;
		while(vals==null){
			temp = readNextLines();
			if(temp[0] == null){
				return false;
			}
			lab = getLabel(temp);
			vals = getFeatureValues(temp);
		}
		
		while(currLabels.size()<numBins-1){
			currLabels.add(lab);
			currFeatures.add(vals);
			currLines.add(temp);
		}
		
		currLabels.add(lab);
		currFeatures.add(vals);
		currLines.add(temp);
		
		if(currLabels.size()>numBins){
			currLabels.removeFirst();
			currFeatures.removeFirst();
			currLines.removeFirst();
		}
		
		return true;
	}
	
	
	public static DataSet replaceNaN(DataSet data) throws WrongAlphabetException, WrongSequenceTypeException, EmptyDataSetException{
		double[] means = new double[data.getElementLength()];
		Arrays.fill(means, Double.MAX_VALUE);
		for(int i=0;i<data.getNumberOfElements();i++){
			Sequence seq = data.getElementAt(i);
			for(int j=0;j<seq.getLength();j++){
				if(!Double.isNaN(seq.continuousVal(j))){
					if(seq.continuousVal(j)<means[j]){
						means[j] = seq.continuousVal(j);
					}
					//means[j] += seq.continuousVal(j)/data.getNumberOfElements();
				}
			}
		}
		
		double[] temp = new double[data.getElementLength()];
		LinkedList<Sequence> seqs = new LinkedList<>();
		for(int i=0;i<data.getNumberOfElements();i++){
			Sequence seq = data.getElementAt(i);
			boolean replace = false;
			for(int j=0;j<seq.getLength();j++){
				if(Double.isNaN(seq.continuousVal(j))){
					replace = true;
				}
			}
			
			if(replace){
				for(int j=0;j<seq.getLength();j++){
					if(Double.isNaN(seq.continuousVal(j))){
						temp[j] = means[j];
					}else{
						temp[j] = seq.continuousVal(j);
					}
				}
				seqs.add(new ArbitrarySequence(seq.getAlphabetContainer(), temp));
			}else{
				seqs.add(seq);
			}
		}
		
		return new DataSet("",seqs);
	}
	
	public Pair<DataSet[],double[][]> getInitialData(HashSet<String> chromosomes) throws FileNotFoundException, IOException, WrongAlphabetException, WrongSequenceTypeException, EmptyDataSetException{
		
		double[][] hist = getPositiveHistogram();
		double[] negh = getNegativeHistogram(hist[0]);
		double[][] sw = getSamplingProbsAndWeights(hist[1], negh);
		double[] breaks = hist[0];
		
		
		double posSum = ToolBox.sum(hist[1]);
		double negSum = ToolBox.sum(negh);
		double pUnif = posSum*10/negSum;
		
		ArrayList<Sequence> pos = new ArrayList<>();
		DoubleList posw = new DoubleList();
		ArrayList<Sequence> neg = new ArrayList<>();
		DoubleList negw = new DoubleList();
		
		Random r = new Random(131);
		
		
		reset();
		while(readNextFeatureVector()){
			char label = getCurrentLabel();
			if(chromosomes == null || chromosomes.contains(getCurrentChromosome())){
				if(label=='S'){
					Sequence seq = getCurrentSequence();
					pos.add(seq);
					posw.add(1.0);
				}else if(label=='U'){
					double p = r.nextDouble();
					if(p<pUnif){
						Sequence seq = getCurrentSequence();
						neg.add(seq);
						negw.add(1.0);
					}

					double dnase = getCurrentDNaseMedian();

					int idx = Arrays.binarySearch(breaks, dnase);
					if(idx < 0){
						idx = -idx-1;
					}
					if(idx >= breaks.length){//fixedBins
						idx = breaks.length-1;
					}

					if(p<sw[0][idx]){
						Sequence seq = getCurrentSequence();
						neg.add(seq);
						negw.add(sw[1][idx]);
					}
				}
			}
		}
		
		
		return new Pair<DataSet[], double[][]>(new DataSet[]{
				replaceNaN(new DataSet("",pos)),
				replaceNaN(new DataSet("",neg))
		}, new double[][]{
			posw.toArray(),
			negw.toArray()
		});
		
	}
	
	
	
	
	/*public static void main(String[] args) throws Exception{
		FeatureReader reader = new FeatureReader(5, args[0], args[1], args[2]);
		
	//	Pair<DataSet[],double[][]> data = reader.getInitialData(null);
		
		IterativeTraining training = new IterativeTraining(reader, 4, getSizes(args[3], 50));
		
		
		HashSet<String> trainChroms = new HashSet<>();
		trainChroms.add("chr1");trainChroms.add("chr2"); trainChroms.add("chr3");
		LinkedList<String> itChroms = new LinkedList<>();
		itChroms.add("chr1");itChroms.add("chr2");
		
		GenDisMixClassifier[] cls = training.iterativeTraining(4, trainChroms, itChroms,0.9);
		
		
		LinkedList<String> predChroms = new LinkedList<>();
		predChroms.add("chr4"); predChroms.add("chr5");
		
		Predictor pred = new Predictor(cls,reader);
		
	//	pred.predict(data.getFirstElement()[0]);
		
		File f = pred.predict(getSizes(args[3], 50),predChroms);
		System.out.println(f.getAbsolutePath());
		
	}*/

	public int getNumBins() {
		return numBins;
	}

	public String getCurrentChromosome() {
		String temp = currLines.get((numBins+1)/2)[1];
		return temp.substring(0, temp.indexOf("\t"));
	}

	public boolean findChr(String chr) throws IOException {
		while(currLines.size()==0||!chr.equals(getCurrentChromosome())){
			if(!readNextFeatureVector()){
				break;
			}
		}
		if(!chr.equals(getCurrentChromosome())){
			reset();
			while(currLines.size()==0||!chr.equals(getCurrentChromosome())){
				if(!readNextFeatureVector()){
					break;
				}
			}
		}
		return chr.equals(getCurrentChromosome());
	}


	public int getCurrentStart() {
		String temp = currLines.get((numBins+1)/2)[1];
		temp = temp.substring(temp.indexOf("\t")+1);
		temp = temp.substring(0, temp.indexOf("\t"));
		return Integer.parseInt(temp);
	}
	
	
	
}
