package projects.dream2016;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.Random;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.Alphabet;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.sequences.ArbitrarySequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.io.ArrayHandler;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.ToolBox.TiedRanks;

/**
 * This class allows for reading {@link DataSet}s in a comfortable way.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class DataParser {

	
	public enum FeatureType{
		Region,
		Distance,
		Length,
		Percent,
		Entropy,
		Score,
		Coverage,
		StrandStrandPos,
		Seq
	}
	
	public enum Collection{
		Each,
		Min,
		Max,
		Mean,
		Median,
		Mix
	}
	
	public enum Bins{
		Each,
		Center3,
		Center,
		Min,
		Max,
		Mean,
		Median,
		MinCenter3,
		MinCenter5,
		MeanCenter3,
		MaxCenter3,
		MaxCenter5,
		LogSum,
		LogSumCenter3
	}
	
	private static class ParserEntry{
		
		private FeatureType ft;
		private Collection coll;
		private Bins bins;
		private int[] columns;
		
		public ParserEntry(FeatureType ft, Collection coll, Bins bins, int... columns){
			this.ft = ft;
			this.coll = coll;
			this.bins = bins;
			this.columns = columns.clone();
		}
		
		public String toString(){
			return ft+" "+coll+" "+bins+" "+Arrays.toString(columns);
		}
		
	}
	
	
	private int skip;
	private int columns;
	private ParserEntry[] entries;
	private Alphabet[] regionAlphabets;
	private AlphabetContainer regionContainer;
	private DiscreteAlphabet strandAlph;
	private DiscreteAlphabet chrAlph = null;
	
	public static ArrayList<String> getConf( String file ) throws IOException {
		ArrayList<String> conf = new ArrayList<String>();
		BufferedReader read = new BufferedReader(new FileReader(file));
		String str;
		while( (str = read.readLine()) != null ){
			if( str.charAt(0) != '#' ) {
				conf.add(str);
			}
		}
		read.close();
		return conf;
	}
	
	public DataParser(String formatFile) throws IOException, IllegalArgumentException, DoubleSymbolException{
		this( getConf(formatFile).toArray(new String[0]) );
	}
	
	public DataParser(String[] lines) throws IOException, IllegalArgumentException, DoubleSymbolException{
		String[] parts = lines[0].split("\t");
		skip = Integer.parseInt(parts[0]);
		columns = Integer.parseInt(parts[1]);
		
		ArrayList<ParserEntry> temp = new ArrayList<DataParser.ParserEntry>();
		for( int i = 1; i < lines.length; i++ ) {
			parts = lines[i].split("\t");
			FeatureType ft = FeatureType.valueOf(parts[0]);
			Collection coll = Collection.valueOf(parts[2]);
			Bins bins = Bins.valueOf(parts[3]);
			parts = parts[1].split(",");
			int[] cols = new int[parts.length];
			for(int j=0;j<cols.length;j++){
				cols[j] = Integer.parseInt(parts[j]);
			}
			temp.add(new ParserEntry(ft,coll,bins,cols));
		}
		
		entries = temp.toArray(new ParserEntry[0]);
		
		regionAlphabets = new Alphabet[]{
				new DiscreteAlphabet(false, "S","-"),
				new DiscreteAlphabet(false, "C","-"),
				new DiscreteAlphabet(false, "U","-"),
				new DiscreteAlphabet(false, "E","-"),
				new DiscreteAlphabet(false, "T","-"),
				new DiscreteAlphabet(false, "s","-"),
				new DiscreteAlphabet(false, "c","-"),
				new DiscreteAlphabet(false, "u","-"),
				new DiscreteAlphabet(false, "e","-"),
				new DiscreteAlphabet(false, "t","-")
		};
		regionContainer = new AlphabetContainer(regionAlphabets);
		
		strandAlph = new DiscreteAlphabet(false, "+","-",".");
		/*
		chrAlph = new DiscreteAlphabet(true, 
				"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
				"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
				"chr21", "chr22", "chr23", "chrX", "chrY");
		/**/
	}
	
	public String toString(){
		return Arrays.toString(entries)+"\n"+columns+" "+skip;
	}
	
	public int getNumberOfColumns() {
		return columns;
	}
	
	public DifferentiableSequenceScore[] getScores(){
		return null;
	}

	public DataSet parseData(String filename) throws Exception {
		return parseData(filename, -1);
	}
	
	public DataSet parseData(String filename, double acceptProb) throws Exception {
		BufferedReader read = new BufferedReader(new FileReader(filename));
		
		String str = null;
		AlphabetContainer cont = null;
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		Random r = new Random();
		while( (str = read.readLine()) != null ){
			if( r.nextDouble() <= acceptProb ) {
				try{
					seqs.add(parse(cont,str.split("\t")));
				}catch(StringIndexOutOfBoundsException ex){
					System.out.println(str);
					//System.out.println(Arrays.toString(parts));
					System.exit(1);
				}	
				if(cont == null){
					cont = seqs.get(0).getAlphabetContainer();
				}
			}
		}	
		
		read.close();
		
		return new DataSet(filename+(acceptProb==1d ? "" : (" randomly selected ~" + (acceptProb*100) + "%")),seqs);
	}
	
	private DataSet parseData(String filename, int max) throws Exception {
		
		BufferedReader read = new BufferedReader(new FileReader(filename));
		
		String str = null;
		AlphabetContainer cont = null;
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		while( (str = read.readLine()) != null && (max<0 || seqs.size()<max) ){
			try{
				seqs.add(parse(cont,str.split("\t")));
			}catch(StringIndexOutOfBoundsException ex){
				System.out.println(str);
				//System.out.println(Arrays.toString(parts));
				System.exit(1);
			}	
			if(cont == null){
				cont = seqs.get(0).getAlphabetContainer();
			}
		}	
		
		read.close();
		
		return new DataSet(filename+(max>0?": first " + max + " sequences":""),seqs);
	}
	
	public ArbitrarySequence parse(AlphabetContainer cont, String[] parts) throws Exception {
		if( cont == null ) {
			cont = getAlphabets(parts.length);
		}
		int numCols = (parts.length-skip)/columns;
		
		DoubleList vals = new DoubleList();
		if( chrAlph != null ) {
			vals.add(chrAlph.getCode(parts[0]));
		}
		
		for(int i=0;i<entries.length;i++){
			int numInOne = entries[i].ft == FeatureType.Region ? regionAlphabets.length : entries[i].columns.length;
			numInOne = entries[i].ft == FeatureType.StrandStrandPos ? 3 : numInOne;
			numInOne = entries[i].ft == FeatureType.Seq ? 50 : numInOne;
			double[][] temp = new double[numInOne][numCols];
			
			for(int j=0;j<entries[i].columns.length;j++){
				
				for(int k=0;k<numCols;k++){
					String str = parts[skip + k*columns + entries[i].columns[j]];
					switch(entries[i].ft){
						case Seq:
							Sequence dnaSeq = Sequence.create(DNAAlphabetContainer.SINGLETON, str);
							for(int l=0;l<dnaSeq.getLength();l++){
								temp[l][k] = dnaSeq.discreteVal(l);
							}
							break;
						case Region:
							Sequence tempSeq = Sequence.create(regionContainer, str);
							for(int l=0;l<tempSeq.getLength();l++){
								temp[l][k] = tempSeq.discreteVal(l);
							}
							break;
						case Score:
							if( str.charAt(0)=='I' ) {//.startsWith("Inf") ) {
								temp[j][k] = Double.NEGATIVE_INFINITY;
							} else {
								temp[j][k] = -Double.parseDouble(str);
							}
							break;
						case Percent:
							if( str.charAt(0)=='I' ) {//.startsWith("Inf") ) {
								temp[j][k] = Double.POSITIVE_INFINITY;
							} else {
								double a = Double.parseDouble(str);
								double b = 1.0-a+1E-10;
								a += 1E-10;
								temp[j][k] = a/(a+b);
							}
							break;
						case Distance://TODO Diskussion
							if(str.length() == 0){
								temp[j][k] = 1E6;
							} else {
								temp[j][k] = Double.parseDouble(str);
								if( temp[j][k]<0 ) {
									temp[j][k]*= -1;
								}
							}
							temp[j][k] = Math.log1p(temp[j][k]);
							break;
						case StrandStrandPos:
							temp[0][k] = strandAlph.getCode( str.substring(0, 1) );
							temp[1][k] = strandAlph.getCode( str.substring(1, 2) );
							temp[2][k] = Double.parseDouble( str.substring(2) );
							break;
						case Coverage:
							try{
								if( str.charAt(0)=='I' ){
									temp[j][k] = Double.POSITIVE_INFINITY;
								}else{
									temp[j][k] = Math.log1p( Double.parseDouble(str) );
								}
							}catch(NumberFormatException ex){
								System.out.println(Arrays.toString(parts));
								throw ex;
							}
							break;
						default:
							temp[j][k] = Double.parseDouble(str);
					}
				}
				
			}
			
			for(int j=0;j<temp.length;j++){
				int m = (temp[j].length-1)/2;
				switch(entries[i].bins){
					case Center: temp[j] = new double[]{temp[j][m]}; break;
					case Center3: temp[j] = new double[]{temp[j][m-1],temp[j][m],temp[j][m+1]}; break;
					case MeanCenter3: temp[j] = new double[]{ToolBox.mean(m-1,m+2,temp[j])}; break;
					case MinCenter3: temp[j] = new double[]{ToolBox.min(m-1,m+2,temp[j])}; break;
					case MinCenter5: temp[j] = new double[]{ToolBox.min(m-2,m+3,temp[j])}; break;
					case MaxCenter3: temp[j] = new double[]{ToolBox.max(m-1,m+2,temp[j])}; break;
					case MaxCenter5: temp[j] = new double[]{ToolBox.max(m-2,m+3,temp[j])}; break;
					case Max: temp[j] = new double[]{ToolBox.max(temp[j])}; break;
					case Mean: temp[j] = new double[]{ToolBox.mean(0,temp[j].length,temp[j])}; break;
					case Median: temp[j] = new double[]{ToolBox.median(temp[j])}; break;
					case Min: temp[j] = new double[]{ToolBox.min(temp[j])}; break;
					case LogSum: temp[j] = new double[]{Normalisation.getLogSum(temp[j])}; break;
					case LogSumCenter3: temp[j] = new double[]{Normalisation.getLogSum(m-1,m+2,temp[j])}; break;
					//case Each: do nothing
				}
			}
			
			
			double[][] temp2 = ToolBox.transpose(temp);
			
			switch(entries[i].coll){
				case Each:
				case Mix:
					for(int j=0;j<temp2.length;j++){
						for(int k=0;k<temp2[j].length;k++){
							vals.add(temp2[j][k]);
						}
					}
					break;
				case Max:
					for(int j=0;j<temp2.length;j++){
						vals.add(ToolBox.max(temp2[j]));
					}
					break;
				case Mean:
					for(int j=0;j<temp2.length;j++){
						vals.add(ToolBox.mean(0,temp2[j].length,temp2[j]));
					}
					break;
				case Median:
					for(int j=0;j<temp2.length;j++){
						vals.add(ToolBox.median(temp2[j]));
					}
					break;
				case Min:
					for(int j=0;j<temp2.length;j++){
						vals.add(ToolBox.min(temp2[j]));
					}
					break;
			}
			
		}
		double[] d = vals.toArray();
		for( int i = 0; i < d.length; i++ ) {
			if( Double.isInfinite(d[i]) ) {
				if( cont.getAlphabetLengthAt(i) == 1d ) {
					//System.out.println("PERCENT " + i + "\t" + d[i]);
					d[i] = 0.5; //percent;
				} else {
					//System.out.println("SCORE " + i + "\t" + d[i]);
					d[i] = 1E2;
				}
			}
		}
		
		return new ArbitrarySequence(cont, d);
	}
	
	private static String toString( String[] parts ) {
		StringBuffer sb = new StringBuffer();
		for( int i = 0; i < parts.length; i++ ) {
			sb.append( (i>0?"\t":"") + parts[i] ); 
		}
		return sb.toString();
	}	

	public AlphabetContainer getAlphabets(int length) throws IllegalArgumentException, DoubleSymbolException {
		length -= skip;
		int num = length/columns;
		ArrayList<Alphabet> alphabets = new ArrayList<Alphabet>();
		if( chrAlph!= null ) {
			alphabets.add(chrAlph);
		}
		for(int i=0;i<entries.length;i++){
			int numCols = entries[i].columns.length;
			
			int numSingle;
			switch(entries[i].coll){
				case Each:
				case Mix: numSingle = numCols; break;
				default: numSingle = 1;
			}
			
			int numTotal;
			switch(entries[i].bins){
				case Each: numTotal = numSingle*num; break;
				default: numTotal = numSingle;
			}
			
			Alphabet[] template;
			switch(entries[i].ft){
				case Seq: template = new Alphabet[entries[i].bins == Bins.Center3 ? 150 : 50];
					for(int j=0;j<template.length;j++){
						template[j] = DNAAlphabet.SINGLETON;
					}
					break;
				case Coverage:
				case Entropy: template = new Alphabet[]{new ContinuousAlphabet(0, Double.MAX_VALUE)}; break;
				case Length: //TODO?
				case Distance: template = new Alphabet[]{new ContinuousAlphabet(0, Integer.MAX_VALUE)}; break;//PROBLEM: Memory
				case Percent: template = new Alphabet[]{new ContinuousAlphabet(0, 1)}; break;
				case Region: template = regionAlphabets; break;
				case Score: template = new Alphabet[]{new ContinuousAlphabet()}; break;//TODO?
				case StrandStrandPos: template = new Alphabet[]{strandAlph,
						strandAlph,
						new ContinuousAlphabet(0, Integer.MAX_VALUE)}; break;
				default: throw new RuntimeException("Not implemented");
			}
			
			for(int j=0;j<numTotal;j++){
				for(int k=0;k<template.length;k++){
					alphabets.add(template[k]);
				}
			}
		}
		Alphabet[] abc = alphabets.toArray(new Alphabet[0]);
		//old return new AlphabetContainer(abc);
		
		HashMap<Alphabet, Integer> nonRed = new HashMap<Alphabet, Integer>();
		int[] assignment = new int[abc.length];
		for( int i = 0, ref = 0; i < abc.length; i++ ) {
			Iterator<Alphabet> it = nonRed.keySet().iterator();
			Integer a = null;
			while(it.hasNext()){
				Alphabet temp = it.next();
				if(temp.compareTo(abc[i])==0){
					a = nonRed.get(temp);
				}
			}
			if( a == null ) {
				a = ref++;
				nonRed.put(abc[i], a);
			}
			assignment[i]=a;
		}
		Alphabet[] abc2 = new Alphabet[nonRed.size()];
		Iterator<Entry<Alphabet, Integer>> it = nonRed.entrySet().iterator();
		while( it.hasNext() ) {
			Entry<Alphabet,Integer> e = it.next();
			abc2[e.getValue()]=e.getKey();
		}
		//System.out.println(Arrays.toString(abc));
		//System.out.println(Arrays.toString(abc2));
		//System.out.println(abc.length  + " vs. " +abc2.length);
		return new AlphabetContainer( abc2, assignment );
	}
	
	public static double[][] getStats( DataSet... d ) throws WrongAlphabetException {
		AlphabetContainer con = d[0].getAlphabetContainer();
		int n = d[0].getNumberOfElements();
		for( int i = 1; i < d.length; i++ ) {
			if( !con.checkConsistency(d[i].getAlphabetContainer()) ) {
				throw new WrongAlphabetException();
			}
			n+=d[i].getNumberOfElements();
		}
			
		double[][] stats = new double[d[0].getElementLength()][];
		double[] values = null;
		for( int p = 0; p < stats.length; p++ ) {
			if( !con.isDiscreteAt(p) ) {
				int c = 0;
				if( values == null ) {
					 values = new double[n];
				}
				for( int i = 0; i < d.length; i++ ) {
					for( int j = 0; j < d[i].getNumberOfElements(); j++, c++ ) {
						values[c] = d[i].getElementAt(j).continuousVal(p);
					}
				}
				
				stats[p] = new double[2];
				stats[p][0] = ToolBox.median(values); //median
				double[] v = new double[values.length];
				for( int k = 0; k < values.length; k++ ) {
					v[k] = Math.abs(values[k]-stats[p][0]);
				}
				stats[p][1] = ToolBox.median(v)*1.4826; //sd ~ mad*1.4826
				if( stats[p][1] == 0 ) {
					System.out.println("SD" + p);
					stats[p][1] = ToolBox.sd(0, values.length, values);
				}
				System.out.println( p + "\t" + Arrays.toString(stats[p]));
			}
		}
		return stats;
	}
	
	public static DataSet zTransform( DataSet d, double[][] stats ) throws Exception {
		AlphabetContainer con = zTransform( d.getAlphabetContainer(), stats );
		ArrayList<Sequence> seqs = new ArrayList<Sequence>();
		for( int j = 0; j < d.getNumberOfElements(); j++ ) {
			seqs.add( zTransform(d.getElementAt(j), stats, con) );
		}
		return new DataSet("normalized " + d.getAnnotation(), seqs);
	}
	
	public static AlphabetContainer zTransform( AlphabetContainer con, double[][] stats ) {
		int[] assign = new int[stats.length];
		HashMap<Integer,Integer> a = new HashMap<Integer, Integer>();
		ArrayList<Alphabet> list = new ArrayList<Alphabet>();
		int set = -1, v=0;
		for( int i = 0; i < stats.length; i++ ) {
			if( stats[i]!=null ) {
				if( set < 0 ) {
					set = v++;
					list.add( new ContinuousAlphabet() );
				}
				assign[i] = set;
			} else {
				Integer k = a.get(con.getAlphabetIndexForPosition(i));
				if( k == null ) {
					k = v++;
					a.put(con.getAlphabetIndexForPosition(i), k);
					list.add( con.getAlphabetAt(i) );
				}
				assign[i] = k;
			}
		}
		return new AlphabetContainer(list.toArray(new Alphabet[0]), assign);
	}
	
	public static Sequence zTransform( Sequence s, double[][] stats, AlphabetContainer con ) throws WrongAlphabetException, WrongSequenceTypeException {
		double[] values = new double[s.getLength()];
		for( int p = 0; p < stats.length; p++ ) {
			values[p] = s.continuousVal(p);
			if( stats[p] != null ) {
				values[p] = (values[p]-stats[p][0]) / stats[p][1];
			}
		}
		return new ArbitrarySequence(con, values);
	}
	

	//Only for tests
	public static void main(String[] args) throws Exception {
		DataParser pars = new DataParser(args[0]);
		DataSet set = pars.parseData(args[1]);
		System.out.println(set);
		System.out.println();
		
		double[][] stats = getStats(set);
		for( int i = 0; i < stats.length; i++ ) {
			System.out.println(i + "\t" + (stats[i]==null?null:Arrays.toString(stats[i])) );
		}
		set = zTransform(set, stats);
		System.out.println();
		System.out.println(set);
		
	}	
	
	
	
	public static enum Weighting {
		ONE,
		DIRECT,
		SIGNAL,
		RANK;
	}
	
	public static double[] getWeights(String fName, Weighting w ) throws IOException {
		DoubleList vals = new DoubleList();
		BufferedReader r = new BufferedReader( new FileReader(fName) );
		String line;
		while( (line=r.readLine()) != null ) {
			vals.add(Double.parseDouble(line));
		}
		r.close();
		double[] res = vals.toArray();
		switch ( w ) {
			case DIRECT:
				break;
			case ONE:
				Arrays.fill(res,1d);
				break;
			case SIGNAL:
				double min = ToolBox.min(res);
				double max = ToolBox.max(res);
				double delta = max-min;
				for( int i = 0; i < res.length; i++ ) {
					res[i] = (2d + (res[i]-min)/delta)/3d;
				}
				break;
			case RANK:
				int[] rank = ToolBox.rank(res, TiedRanks.SPORTS );
				double m = res.length;
				for( int i = 0; i < res.length; i++ ) {
					res[i] = (2d + (m-rank[i])/m)/3d;
					
				}
				break;
		}
		return res;
	}
	
}