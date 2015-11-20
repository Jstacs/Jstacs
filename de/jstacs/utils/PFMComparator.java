/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */
package de.jstacs.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.Random;
import java.util.AbstractMap.SimpleEntry;

import de.jstacs.classifiers.utils.PValueComputation;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.ComplementableDiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;

/**
 * This class implements a number of methods for the comparison of position frequency matrices (PFMs) as described in the
 * <a href="http://genome.cshlp.org/content/18/7/1180.short">Amadeus paper</a> 
 * 
 * @author Jens Keilwagen
 */
//* and
//* <a href="http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000010">motif comparison paper</a>. 
public class PFMComparator {

	/**
	 * Returns a string representation of the matrix, where each row of the matrix
	 * is printed on one line and columns are separated by tabstops.
	 * @param matrix the matrix
	 * @return the string representation
	 */
	public static String matrixToString( double[][] matrix ) {
		StringBuffer sb = new StringBuffer();
		for( int j, i = 0; i < matrix.length; i++ ) {
			sb.append( matrix[i][0] );
			for( j = 1; j < matrix[i].length; j++ ) {
				sb.append( '\t' );
				sb.append( matrix[i][j] );
			}
			sb.append( '\n' );
		}
		return sb.toString();
	}
	
	/**
	 * This method extracts the 
	 * The method does not use any degenerated IUPAC code.
	 * 
	 * @param con the {@link AlphabetContainer} that was used to create the PFM
	 * @param pfm the position frequency matrix (or also possible the position weight matrix)
	 * 
	 * @return the consensus which is the string with most probable symbol at each position 
	 */
	public static String getConsensus( AlphabetContainer con, double[][] pfm ) {
		String c = "";
		for( int m, p, l = 0; l < pfm.length; l++ ) {
			m = 0;
			for( p = 1; p < pfm[l].length; p++ ) {
				if( pfm[l][m] < pfm[l][p] ) {
					m = p;
				}
			}
			c += con.getSymbol( l, m );
		}
		return c;
	}
	
	/**
	 * This method counts the occurrences of symbols in the given data sets. This method can therefore be used to create a
	 * background distribution.
	 * 
	 * @param data the data sets
	 * 
	 * @return an array containing the counts
	 */
	public static double[] getCounts( DataSet ... data ) {
		double[] counts = new double[(int)data[0].getAlphabetContainer().getAlphabetLengthAt( 0 )];
		Sequence seq;
		for( int d = 0; d < data.length; d++ ) {
			for( int n = 0; n < data[d].getNumberOfElements(); n++ ) {
				seq = data[d].getElementAt( n );
				for( int l = 0; l < seq.getLength(); l++ ) {
					counts[seq.discreteVal( l )]++;
				}
			}
		}
		return counts;
	}
	
	/**
	 * This method enables the user to normalize a array containing counts. After using this method the array contain
	 * probabilities.
	 * 
	 * @param counts the array of counts
	 */
	public static void normalize( double[] counts ) {
		double sum = 0;
		for( int l = 0; l < counts.length; l++ ) {
			sum += counts[l];
		}
		for( int l = 0; l < counts.length; l++ ) {
			counts[l] /= sum;
		}
	}
	
	/**
	 * This method creates a PFM from a {@link DataSet} of {@link Sequence}s.
	 * 
	 * @param data the {@link DataSet}
	 * 
	 * @return the PFM
	 */
	public static double[][] getPFM( DataSet data ) {
		return getPFM( data, 0, data.getNumberOfElements() );
	}
	
	/**
	 * Returns a position frequency matrix (PFM, rows=positions, columns=symbols) for the given subset of {@link DataSet}.
	 * @param data the data set
	 * @param start the first sequence to consider
	 * @param end the first sequence not not consider
	 * @return the PFM
	 */
	public static double[][] getPFM( DataSet data, int start, int end ) {
		if( data == null ) {
			return null;
		} else {
			double[][] pfm = new double[data.getElementLength()][];
			AlphabetContainer con = data.getAlphabetContainer();
			for( int l = 0; l < pfm.length; l++ ) {
				pfm[l] = new double[(int)con.getAlphabetLengthAt( l )];
			}
			Sequence seq;
			for( int n = start; n < end; n++ ) {
				seq = data.getElementAt( n );
				for( int l = 0; l < pfm.length; l++ ) {
					pfm[l][seq.discreteVal( l )]++;
				}
			}
			return pfm;
		}
	}
	
	/**
	 * Returns a position weight matrix (PWM, rows=positions, columns=symbols, containing probabilities) for the given subset of {@link DataSet}.
	 * @param data the data set
	 * @param start the first sequence to consider
	 * @param end the first sequence not not consider
	 * @return the PWM
	 */
	public static double[][] getPWM( DataSet data, int start, int end ){
		double[][] pwm = getPFM( data, start, end );
		for(int i=0;i<pwm.length;i++){
			normalize( pwm[i] );
		}
		return pwm;
	}
	
	/**
	 * This method creates a PFM from a {@link DataSet} of {@link Sequence}s.
	 * 
	 * @param data the {@link DataSet}
	 * @param weights the weights on the sequences in <code>data</code>
	 * 
	 * @return the PFM
	 */
	public static double[][] getPFM( DataSet data, double[] weights ) {
		if( data == null ) {
			return null;
		} else {
			double[][] pfm = new double[data.getElementLength()][];
			AlphabetContainer con = data.getAlphabetContainer();
			for( int l = 0; l < pfm.length; l++ ) {
				pfm[l] = new double[(int)con.getAlphabetLengthAt( l )];
			}
			Sequence seq;
			for( int n = 0; n < data.getNumberOfElements(); n++ ) {
				seq = data.getElementAt( n );
				for( int l = 0; l < pfm.length; l++ ) {
					pfm[l][seq.discreteVal( l )]+=weights[n];
				}
			}
			return pfm;
		}
	}
	
	/**
	 * This method returns the PFM that is the reverse complement of the given PFM.
	 * 
	 * @param abc the {@link de.jstacs.data.alphabets.ComplementableDiscreteAlphabet} that is used to create the reverse complement
	 * @param pfm the absolute frequencies for each position and symbol; {@latex.inline $pfm[l][a] := N_{X_\\ell=a}$}
	 * 
	 * @return the PFM that is the reverse complement of the given PFM
	 */
	public static double[][] getReverseComplement( ComplementableDiscreteAlphabet abc, double[][] pfm ) {
		int length = pfm.length, symbols = (int) abc.length();
		double[][] pfmRC = new double[length][symbols];
		for( int l = 0; l < length; l++ ) {
			for( int a = 0; a < symbols; a++ ) {
				pfmRC[l][a] = pfm[length-1-l][abc.getComplementaryCode(a)];
			}
		}
		return pfmRC;
	}
	
	/**
	 * This method reads a number of PFMs from a directory as downloaded from the Uniprobe data base and 
	 * returns them in an {@link ArrayList} together with some annotation.
	 * 
	 * @param startdir the directory containing the PFMs in Uniprobe format
	 * 
	 * @return an {@link ArrayList} containing {@link SimpleEntry}s with annotation and PFM
	 * 
	 * @throws IOException if something went wrong while reading a file
	 */
	public static ArrayList<SimpleEntry<String, double[][]>> readPFMsFromUniprobe( String startdir ) throws IOException {
		
		LinkedList<File> dirs = new LinkedList<File>();
		dirs.add( new File(startdir) );
		
		ArrayList<SimpleEntry<String, double[][]>> res = new ArrayList<SimpleEntry<String, double[][]>>();
		
		ArrayList<String> strs = new ArrayList<String>();
		
		while(dirs.size() > 0){
			File dir = dirs.pop();
			File[] files = dir.listFiles();
			for(int i=0;i<files.length;i++){
				if(files[i].isDirectory() && !files[i].isHidden()){
					dirs.add( files[i] );
				}else if(files[i].getName().endsWith( ".txt" ) || files[i].getName().endsWith( ".pwm" )){
					res.add( readPFMFromUniprobe( dir.getName(), files[i] ) );
				}
			}
		}
		return res;
	}
	
	
	/**
	 * Reads a PFM from the UniProbe format.
	 * @param annotPrefix a prefix on the PFM's annotation 
	 * @param file the {@link File} to read
	 * @return the PFM and annotation
	 * @throws IOException
	 */
	public static SimpleEntry<String, double[][]> readPFMFromUniprobe(String annotPrefix, File file) throws IOException {
		String annot = annotPrefix+file.getName().substring( 0, file.getName().lastIndexOf( '.' ) );
		BufferedReader read = new BufferedReader( new FileReader( file ) );
		String temp = read.readLine();
		annot += ": "+(temp == null ? "" : temp.trim() );
		
		String line = null;
		ArrayList<String> strs = new ArrayList<String>();
		while( (line = read.readLine()) != null ){
			line = line.trim();
			if(line.length() > 0){
				strs.add( line );
			}
		}
		
		read.close();
		
		String[] strings = new String[4];
		for(int j=strs.size()-4,k=0;j<strs.size();j++,k++){
			strings[k] = strs.get( j );
			if(strings[k].indexOf( ':' ) >= 0){
				strings[k] = strings[k].substring( strings[k].indexOf( ':' )+1 ).trim();
			}
		}
		
		String[] parts = strings[0].split( "\\s+" );
		double[][] pwm = new double[parts.length][4];
		for(int j=0;j<strings.length;j++){
			parts = strings[j].split( "\\s+" );
			if(parts.length != pwm.length){
				throw new RuntimeException();
			}
			for(int k=0;k<parts.length;k++){
				pwm[k][j] = Double.parseDouble( parts[k] );
			}
		}
		return new SimpleEntry<String, double[][]>( annot, pwm );
	}
	
	/**
	 * This method reads a number of PFMs from a file in the Jaspar FastA-like format and 
	 * returns them in an {@link ArrayList} together with some annotation.
	 * 
	 * @param filename the path to the file containing the Jaspar PFMs
	 * 
	 * @return an {@link ArrayList} containing {@link SimpleEntry}s with annotation and PFM
	 * 
	 * @throws IOException if something went wrong while reading the file
	 */
	public static ArrayList<SimpleEntry<String,double[][]>> readPFMsFromJasparFastA( String filename ) throws IOException {
		
		BufferedReader r = new BufferedReader( new FileReader( new File( filename ) ) );
		ArrayList<SimpleEntry<String,double[][]>> res = readPFMsFromJasparFastA(r);
		r.close();
		return res;
	}
	/**
	 * Reads a list of PFMs from the Jaspar format
	 * @param r the {@link BufferedReader} on the Jaspar file
	 * @return the PFMs
	 * @throws IOException if the contents of <code>r</code> could not be read
	 */
	public static  ArrayList<SimpleEntry<String,double[][]>> readPFMsFromJasparFastA( BufferedReader r ) throws IOException {
		ArrayList<SimpleEntry<String, double[][]>> res = new ArrayList<SimpleEntry<String, double[][]>>();
		String line, annot = null;
		double[][] mat = null;
		int i=0;
		while( (line = r.readLine()) != null ){
			line = line.trim();
			if(line.length() > 0){
				if(line.startsWith( ">" )){
					if(mat != null){
						res.add( new SimpleEntry<String, double[][]>( annot, mat ) );
						annot = null;
						mat = null;
					}

					annot = line.substring( 1 ).trim();
					i=0;
				}else{
					line = line.substring( line.indexOf( '[' )+1, line.lastIndexOf( ']' ) ).trim();
					String[] parts = line.split( "\\s+" );
					if(mat == null){
						mat = new double[parts.length][4];
					}
					for(int j=0;j<parts.length;j++){
						mat[j][i] = Double.parseDouble( parts[j] );
					}
					i++;
				}
			}
		}
		if(mat != null){
			res.add( new SimpleEntry<String, double[][]>( annot, mat ) );
		}
		return res;
	}
	
	/**
	 * This method reads a number of PFMs from a file and return them in an {@link ArrayList} together with some annotation.
	 * 
	 * @param fileName the name of the file containing the PFMs in EMBL format
	 * @param max the maximal number of PFMs that should be read 
	 * 
	 * @return an {@link ArrayList} containing {@link SimpleEntry}s with annotation and PFM
	 * 
	 * @throws IOException if something went wrong while reading the file
	 */
	public static ArrayList<SimpleEntry<String, double[][]>> readPFMsFromEMBL( String fileName, int max ) throws IOException {
		ArrayList<SimpleEntry<String, double[][]>> res = new ArrayList<SimpleEntry<String, double[][]>>();
		ArrayList<double[]> pfm = new ArrayList<double[]>();
		double[][] tooSmall = new double[0][];
		BufferedReader r = new BufferedReader( new FileReader( new File( fileName ) ) );
		String line, annot;
		String[] split;
		double[] positionalFrequencies;
		int a;
		while( res.size() < max && (line = r.readLine()) != null ) {
			if( line.startsWith( "AC" ) ) {
				annot = line.substring( line.indexOf(' ') ).trim();
				while( !(line = r.readLine()).startsWith( "ID" ) ) {}
				annot += " (" + line.substring( line.indexOf(' ') ).trim() + ")";
				while( !(line = r.readLine()).startsWith( "01" ) ) {}
				pfm.clear();
				do {
					split = line.split( "[ ]+" );
					positionalFrequencies = new double[split.length - 2];
					for( a = 0; a < positionalFrequencies.length; a++ ) {
						positionalFrequencies[a] = Double.parseDouble( split[1+a] );
					}
					pfm.add( positionalFrequencies );
				} while( !(line = r.readLine()).startsWith( "XX" ) );
				res.add( new SimpleEntry<String, double[][]>( annot, pfm.toArray( tooSmall ) ) );
				while( !(line = r.readLine()).startsWith( "//" ) ) {}
			}
		}
		r.close();
		return res;
	}
	
	/**
	 * Reads a set of PFM from the Transfac format
	 * @param fileName the Transfac file
	 * @param max the maximum number of PFMs to read
	 * @return the PFMs
	 * @throws IOException if the file could not be read
	 */
	public static ArrayList<SimpleEntry<String, double[][]>> readPFMsFromTransfac( String fileName, int max ) throws IOException {
		ArrayList<SimpleEntry<String, double[][]>> res = new ArrayList<SimpleEntry<String, double[][]>>();
		ArrayList<double[]> pfm = new ArrayList<double[]>();
		double[][] tooSmall = new double[0][];
		BufferedReader r = new BufferedReader( new FileReader( new File( fileName ) ) );
		String line, annot;
		String[] split;
		double[] positionalFrequencies;
		int a;
		while( res.size() < max && (line = r.readLine()) != null ) {
			while( !line.startsWith( "DE" ) ) { line = r.readLine(); }
			annot = line.substring( line.indexOf("DE")+2 ).trim();
			while( !(line = r.readLine()).matches( "^[0-9]+.*" ) ) {}
			pfm.clear();
			do {
				split = line.split( "[ \\t]+" );
				positionalFrequencies = new double[split.length - 2];
				for( a = 0; a < positionalFrequencies.length; a++ ) {
					positionalFrequencies[a] = Double.parseDouble( split[1+a] );
				}
				pfm.add( positionalFrequencies );
			} while( !(line = r.readLine()).startsWith( "XX" ) );
			res.add( new SimpleEntry<String, double[][]>( annot, pfm.toArray( tooSmall ) ) );

		}
		r.close();
		return res;
	}
	
	private static Random r = new Random();
	
	/**
	 * This methods finds for a user specified PFM <code>pfm</code> similar PFMs in a list of known PFMs.
	 * @param abc the alphabet of this PFM
	 * @param pfm the PFM to be compare to the known PFMs
	 * @param knownPFMs a {@link ArrayList} of known PFMs
	 * @param distance a distance measure to assess the similarity of the PFMs
	 * @param minOverlap the minimal number of consecutive positions in an alignment
	 * @param allowMiniShift the minimal number of allowed shifts (even if the minimal overlap is bigger)
	 * @param pValues compute p-values for matrix similarity
	 * @param threshold only hits with a distance below this threshold will be reported
	 * 
	 * @return an array of similar PFMs and their annotations
	 * 
	 * @see PFMComparator#readPFMsFromEMBL(String, int)
	 */
	public static ComparableElement<String, Double>[] find( ComplementableDiscreteAlphabet abc, double[][] pfm, ArrayList<SimpleEntry<String, double[][]>> knownPFMs, PFMDistance distance, int minOverlap, int allowMiniShift, boolean pValues, double threshold ) {
		ArrayList<ComparableElement<String, Double>> list = new ArrayList<ComparableElement<String, Double>>(10);
		double[][] currentCand, pfmRC = getReverseComplement( abc, pfm );
		int myMinimalOverlap, m = Math.min(minOverlap, pfm.length - allowMiniShift ), columns = 0;
		Hashtable<Integer, double[]> scores = new Hashtable<Integer, double[]>();
		double[] values;
		double[][] randomPFM;
		double d;
		if( pValues ) {
			for( int i = 0; i < knownPFMs.size(); i++ ) {
				columns += knownPFMs.get( i ).getValue().length;
			}
		}
		for( int i = 0; i < knownPFMs.size(); i++ ) {
			currentCand = knownPFMs.get(i).getValue();
			myMinimalOverlap = Math.min( m, currentCand.length - allowMiniShift );
			d = Math.min(
					distance.compare( pfm, currentCand, myMinimalOverlap ),
					distance.compare( pfmRC, currentCand, myMinimalOverlap )
				);
			
			if( pValues ) {
				values = scores.get( currentCand.length );
				if( values == null ) {
					values = new double[1000];//TODO
					randomPFM = new double[currentCand.length][];
					for( int n = 0; n < values.length; n++ ) {
						for( int idx, c, len, l = 0; l < randomPFM.length; l++ ) {
							c = r.nextInt( columns );
							idx = 0;
							while( (len = knownPFMs.get( idx ).getValue().length) <= c ) {
								c -= len;
								idx++;
							}
							randomPFM[l] = knownPFMs.get( idx ).getValue()[c];
						}
						values[n] =	Math.min(
								distance.compare( pfm, randomPFM, myMinimalOverlap ),
								distance.compare( pfmRC, randomPFM, myMinimalOverlap )
							);
					}
					Arrays.sort( values );
					scores.put( currentCand.length, values );
				}
				d = PValueComputation.getPValue( values, d );	
			}
			
			if( d <= threshold ) {
				list.add( new ComparableElement<String, Double>( knownPFMs.get(i).getKey(), d ) );
			}
		}
		ComparableElement<String, Double>[] res = list.toArray( new ComparableElement[0] );
		Arrays.sort( res );
		return res;
	}
	
	/**
	 * This interface declares a method for comparing different PFMs.
	 * If the PFMs are identical the value 0 should be returned, otherwise some positive value.
	 * 
	 * Ideally, each distance has to be length-normalized.
	 *   
	 * @author Jens Keilwagen
	 */
	public static abstract class PFMDistance {
		/**
		 * This method computes the distance between two PFMs.
		 * If the PFMs are identical the value 0 should be returned, otherwise some positive value.
		 * 
		 * @param pfm1 the first PFM
		 * @param pfm2 the second PFM
		 * @param offset the offset for the alignment of the PFMs
		 * 
		 * @return the distance between the PFMs
		 */
		public final double getDistance( double[][] pfm1, double[][] pfm2, int offset ) {
			if( offset >= 0 ) {
				return getDistance( pfm1, pfm2, offset, 0 );
			} else {
				return getDistance( pfm1, pfm2, 0, -offset );
			}
		}
		
		/**
		 * Computes the mean distance between the overlapping parts of <code>pfm1</code> and <code>pfm2</code> starting at the offsets
		 * <code>l1</code> and <code>l2</code>, respectively.
		 * @param pfm1 the first PFM
		 * @param pfm2 the second PFM
		 * @param l1 the offset for the first PFM
		 * @param l2 the offset for the second PFM
		 * @return the mean distance of the overlapping part
		 */
		protected abstract double getDistance( double[][] pfm1, double[][] pfm2, int l1, int l2 );
		
		/**
		 * This method compares two PFMs, <code>pfm1</code> and <code>pfm2</code>.
		 * The method tests all alignments of the PFMs with an overlap of <code>mininmalOverlap</code> consecutive positions.
		 * 
		 * @param pfm1 the first PFM
		 * @param pfm2 the second PFM
		 * @param minimalOverlap the minimal number of consecutive positions in an alignment
		 * 
		 * @return the distance
		 */
		public double compare( double[][] pfm1, double[][] pfm2, int minimalOverlap ) {
			int l1 = pfm1.length, l2 = pfm2.length;
			if( l1 < minimalOverlap || l2 < minimalOverlap ) {
				throw new IllegalArgumentException( "The PFM are too small to have a minimal overlap of " + minimalOverlap + "." );
			}
			double current, best = Double.POSITIVE_INFINITY;
			for( int offset = minimalOverlap - l2; offset <= l1 - minimalOverlap ; offset++ ) {
				current = getDistance( pfm1, pfm2, offset );
				if( current < best ) {
					best = current;
				}
			}
			return best;
		}
	}
	
	/**
	 * Wraps a given {@link PFMDistance} and pads the considered PFMs with uniformly distributed positions.
	 * @author Jan Grau
	 *
	 */
	public static class UniformBorderWrapper extends PFMDistance {

		private PFMDistance inner;
		
		/**
		 * Creates a new {@link UniformBorderWrapper} for a given {@link PFMDistance}.
		 * @param inner the internal {@link PFMDistance} to be padded
		 */
		public UniformBorderWrapper(PFMDistance inner){
			this.inner = inner;
		}
		
		@Override
		protected double getDistance(double[][] pfm1, double[][] pfm2, int l1,
				int l2) {
			
			System.out.print(l1+" "+l2+" ");
			
			double tempDist = inner.getDistance(pfm1, pfm2, l1, l2);
			
			System.out.print(tempDist+" ");
			
			double[][] p1 = null;
			double[][] p2 = null;
			
			if(l1 > 0 && l2 > 0){
				int mi = Math.min(l1, l2);
				l1 -= mi;
				l2 -= mi;
			}
			if(l1 > l2){
				double[][] temp = pfm1;
				pfm1 = pfm2;
				pfm2 = temp;
				int tl = l1;
				l1 = l2;
				l2 = tl;
			}
			int maxlen = Math.max(pfm1.length+l2,pfm2.length+l1);//l1 should be zero
			p1 = new double[maxlen][];
			p2 = new double[maxlen][];
			double n1 = ToolBox.sum(pfm1[0]);
			double n2 = ToolBox.sum(pfm2[0]);
			double[] c1 = new double[pfm1[0].length];
			for(int i=0;i<c1.length;i++){
				c1[i] = n1 / c1.length;
			}
			double[] c2 = new double[pfm2[0].length];
			for(int i=0;i<c2.length;i++){
				c2[i] = n2 / c2.length;
			}
			
			for(int i=0;i<maxlen;i++){
				if(i < l2 || i >= pfm1.length){
					p1[i] = c1.clone();
				}else{
					p1[i] = pfm1[i-l2];
				}
				if(i < l1 || i >= pfm2.length){
					p2[i] = c2.clone();
				}else{
					p2[i] = pfm2[i-l1];
				}
			}			
			
			
			
			double dist = inner.getDistance(p1, p2, 0, 0);
			System.out.println(dist);
			return dist;
		}
		
		
		
	}
	
	/**
	 * This class implements the normalized Euclidean distance. 
	 * 
	 * @author Jens Keilwagen
	 */
	public static class NormalizedEuclideanDistance extends PFMDistance {
		protected double getDistance( double[][] pfm1, double[][] pfm2, int l1, int l2 ) {
			double res = 0, sum, diff, sum1, sum2;
			int a, l = 0;
			for( ; l1 < pfm1.length && l2 < pfm2.length; l1++, l2++, l++  ) {
				if( pfm1[l1].length != pfm2[l2].length ) {
					throw new IllegalArgumentException( "The PFMs are not comparable at column " + l1 + " and " + l2 + "." );
				}
				sum1 = sum2 = 0;
				for( a = 0; a < pfm1[l1].length; a++ ) {
					sum1 += pfm1[l1][a];
					sum2 += pfm2[l2][a];
				}
				sum = 0;
				for( a = 0; a < pfm1[l1].length; a++ ) {
					diff = pfm1[l1][a]/sum1 - pfm2[l2][a]/sum2;
					sum += diff*diff;
				}
				res += Math.sqrt( sum );
			}
			return res/(Math.sqrt(2) * l);
		}
	}
	
	/**
	 * This class implements the Pearson correlation coefficient.
	 * For the reason of creating a distance, it returns 1 - CC.
	 *  
	 * @author Jens Keilwagen
	 */
	public static class OneMinusPearsonCorrelationCoefficient extends PFMDistance {
		protected double getDistance( double[][] pfm1, double[][] pfm2, int l1, int l2 ) {
			double res = 0, d1, d2, s, s1, s2, sum1, sum2;
			int a, l = 0;
			for( ; l1 < pfm1.length && l2 < pfm2.length; l1++, l2++, l++  ) {
				if( pfm1[l1].length != pfm2[l2].length ) {
					throw new IllegalArgumentException( "The PFMs are not comparable at column " + l1 + " and " + l2 + "." );
				}
				sum1 = sum2 = 0;
				for( a = 0; a < pfm1[l1].length; a++ ) {
					sum1 += pfm1[l1][a];
					sum2 += pfm2[l2][a];
				}
				s = s1 = s2 = 0;
				for( a = 0; a < pfm1[l1].length; a++ ) {
					d1 = pfm1[l1][a]/sum1 - 0.25;
					d2 = pfm2[l2][a]/sum2 - 0.25;
					s += d1*d2;
					s1 += d1*d1;
					s2 += d2*d2;
				}
				if(s1 > 0 && s2 > 0){//TODO sd=0
					res += s / Math.sqrt( s1 * s2 );
				}
			}
			return 1d - 1d/(double) l * res;
		}
	}
	
	/**
	 * This class implements the symmetric Kullback-Leibler-divergence.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class SymmetricKullbackLeiblerDivergence extends PFMDistance {
		private double ess;
		
		/**
		 * This constructor creates a new instance with a given value for the equivalent sample size.
		 * 
		 * @param ess the equivalent sample size
		 * 
		 * @throws IllegalArgumentException if the <code>ess</code> &lt; 0
		 */
		public SymmetricKullbackLeiblerDivergence( double ess ) throws IllegalArgumentException {
			if( ess < 0 ) {
				throw new IllegalArgumentException( "The ess has to be non-negative." );
			}
			this.ess = ess;
		}
		
		protected double getDistance( double[][] pfm1, double[][] pfm2, int l1, int l2 ) {
			double res = 0, p1_a, p2_a, pc, sum1, sum2;
			int a, l = 0;
			for( ; l1 < pfm1.length && l2 < pfm2.length; l1++, l2++, l++  ) {
				if( pfm1[l1].length != pfm2[l2].length ) {
					throw new IllegalArgumentException( "The PFMs are not comparable at column " + l1 + " and " + l2 + "." );
				}
				sum1 = sum2 = ess;
				pc = ess/(double)pfm1[l1].length;
				for( a = 0; a < pfm1[l1].length; a++ ) {
					sum1 += pfm1[l1][a];
					sum2 += pfm2[l2][a];
				}
				for( a = 0; a < pfm1[l1].length; a++ ) {
					p1_a = (pc+pfm1[l1][a])/(double)sum1;
					p2_a = (pc+pfm2[l2][a])/(double)sum2;
					res += (p1_a - p2_a) * Math.log( p1_a / p2_a );
				}
			}
			return 1d/(2d * l) * res;
		}
	}
	
	
	
	/**
	 * This class implements the BLiC score from the
	 * <a href="http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000010">motif comparison paper</a>. 
	 * 
	 * This implementation uses a Dirichlet prior.
	 * 
	 * @author Jens Keilwagen
	 */
//	public static class BLiC extends PFMDistance {
//
//		private double ess, relaxingFactor;
//		private double[] log_p_bg;
//		
//		/**
//		 * This constructor creates an instance the used a Dirichlet prior (BDeu) with given <code>ess</code>
//		 * and <code>p_bg</code> to compute the BLiC score.
//		 * 
//		 * @param ess the equivalent sample size
//		 * @param p_bg the background distribution
//		 * @param relaxingFactor the factor is used to weight the unaligned positions, in the paper it is set to 0.2
//		 * 
//		 * @throws IllegalArgumentException if the <code>ess</code> &lt; 0
//		 */
//		public BLiC( double ess, double[] p_bg, double relaxingFactor ) throws IllegalArgumentException {
//			if( ess < 0 ) {
//				throw new IllegalArgumentException( "The ess has to be non-negative." );
//			}
//			this.ess = ess;
//			log_p_bg = new double[p_bg.length];
//			double sum = 0;
//			for( int a = 0; a < p_bg.length; a++ ) {
//				sum += p_bg[a];
//				if( p_bg[a] > 0 && p_bg[a] < 1 ) {
//					log_p_bg[a] = Math.log( p_bg[a] );
//				} else {
//					throw new IllegalArgumentException( "The background distribution is not in (0,1) for the " + a + "-th symbol." );
//				}
//			}
//			if( Math.abs( 1d-sum ) > 1E-7 ) {
//				throw new IllegalArgumentException( "The background distribution does not sum to 1." );
//			}
//			if( relaxingFactor < 0 || relaxingFactor > 1 ) {
//				
//			}
//			this.relaxingFactor = relaxingFactor;
//		}
//		
//		public double getDistance( double[][] pfm1, double[][] pfm2, int offset ) {
//			int l1 = 0, l2 = 0;
//			double res = 0;
//			//TODO unaligned position on the right side
//			if( offset > 0 ) {/*
//				while( l1 < offset ) {
//					res += relaxingFactor*getBLiCForPosition( pfm1[l1++], //TODO );
//				}
//				//l1 is now offset;/**/
//				l1 = offset;//TODO remove
//			} else {/*
//				while( l2 < -offset ){
//					res += relaxingFactor*getBLiCForPosition( //TODO, pfm2[l2++] );
//				}
//				//l2 is now -offset;/**/
//				l2 = -offset; //TODO remove
//			}
//			for( ; l1 < pfm1.length && l2 < pfm2.length; l1++, l2++  ) {
//				res += getBLiCForPosition( pfm1[l1], pfm2[l2] );
//			}
//			/*
//			//TODO unaligned position on the right side
//			while( l1 < pfm1.length ) {
//				res += relaxingFactor*getBLiCForPosition( pfm1[l1++], //TODO );
//			}
//			while( l2 < pfm2.length ){
//				res += relaxingFactor*getBLiCForPosition( //TODO, pfm2[l2++] );
//			}
//			/**/
//			return res;
//		}
//		
//		double getBLiCForPosition( double[] distributiuon1, double[] distribution2 ) {
//			if( distributiuon1.length != distribution2.length ) {
//				throw new IllegalArgumentException( "The PFMs are not comparable." );
//			}
//			double sum1 = 0, sum2 = 0, sum3, pc = ess/(double)distributiuon1.length;
//			int a;
//			for( a = 0; a < distributiuon1.length; a++ ) {
//				sum1 += distributiuon1[a];
//				sum2 += distribution2[a];
//			}
//			sum3 = sum1 + sum2 + ess;
//			sum1 += ess;
//			sum2 += ess;
//			
//			double res = 0;
//			for( a = 0; a < distributiuon1.length; a++ ) {
//				res += (distributiuon1[a]+distribution2[a]) * ( 2*Math.log( (distributiuon1[a]+distribution2[a]+pc)/sum3 ) - log_p_bg[a]) 
//					- distributiuon1[a] * Math.log( (distributiuon1[a]+pc)/sum1 )
//					- distribution2[a] * Math.log( (distribution2[a]+pc)/sum2 );
//			}
//			return res;
//		}
//	}
}
