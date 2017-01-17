package projects.dream2016;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.ArbitrarySequence;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Time;

public class GenomeWideScan {
	
	/**
	 * Performs a genomewide scan
	 * 
	 * @param args
	 * 0	= conf file
	 * 1	= out file
	 * 2	= ? or chromsome list
	 * 3...n	= classifier
	 * n+1	= windows
	 * n+2... = feature files
	 * 
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		
		DataParser pars = new DataParser(args[0]);
		AlphabetContainer origCon = null, transCon = null;
		
		ArrayList<AbstractScoreBasedClassifier> cl = new ArrayList<AbstractScoreBasedClassifier>();
		int offset = 3;
		while( new File(args[offset]).exists() ) {
			cl.add( (AbstractScoreBasedClassifier) XMLParser.extractObjectForTags(FileManager.readFile(args[offset]), "classifier") );
			offset++;
		}
		
		HashSet<String> chrs = null;
		if( !args[2].equals("?") ) {
			chrs = new HashSet<String>();
			String[] chr = args[2].split(",");
			for( int i = 0; i < chr.length; i++ ) {
				chrs.add(chr[i]);
			}
			System.out.println( "Try to make predictions only for: " + chrs );
		}
		
		//parse data
		int windows = Integer.parseInt(args[offset++]), c=0, start = 0, last=-1, a=-1, p=-1000;
		String chr = null;
		
		//open all needed files
		int l = args.length-offset;
		BufferedReader[] reader = new BufferedReader[l];
		for(int f=0;f<l;f++){
			GZIPInputStream stream = new GZIPInputStream(new FileInputStream(args[f+offset]));
			reader[f] = new BufferedReader(new InputStreamReader(stream));
		}
		
		String[] line = new String[l], li = new String[l];
		int[] remove = null;
		String[] split = null;
		double[] scores = new double[2];
		double ls = 0;
		boolean any;
		boolean[] read = new boolean[l];
		Arrays.fill( read, true );
		Time t = Time.getTimeInstance(null);
		//BufferedWriter w = new BufferedWriter( new FileWriter( args[1] ) );
		FileOutputStream output = new FileOutputStream(args[1]);
		BufferedWriter w = new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8") );
		
		w.append("#chr\tpos");
		for(int f=3;f<offset;f++){
			w.append("\t"+args[f]);
		}
		w.newLine();
		
		boolean use = false;
		outerloop: while( true ){
			any = false; 
			for(int f=0;f<l;f++){
				if( read[f] ) {
					line[f] = reader[f].readLine();
					read[f] = line[f]!=null && line[f].charAt(0) != '[';
				}
				
				if( !read[f] ) {
					//default values
					if( args[f+offset].endsWith("hg19.fa_other.txt.gz") ) {
						li[f] = "0.5\t0.125\t0.0\t0.0\t0.0";
					} else if( args[f+offset].endsWith("hg19.fa_tracts.txt.gz") ) {
						li[f] = "0\t0\t0\t0";
					} else if( args[f+offset].endsWith("gencode.v19.types.txt.gz") ) {
						li[f] = chr+"\t"+start+"\t----------";
					} else if( args[f+offset].endsWith(".bigwig-interval.txt.gz") ) {
						li[f] = "0\t0\t0\t0";
					} else if( args[f+offset].endsWith(".xml_winscores2.txt.gz") ) {
						li[f] = "Inf\tInf\t0\t0\t..0";
					} else if( args[f+offset].endsWith("_winscores2.txt.gz") ) {
						li[f] = "Inf\tInf\t0\t0\t..0";
					} else if( args[f+offset].endsWith("hg19.genome.fa.seqs.gz") ){
						li[f] = "CTTAGCGGAAATAGGAGAAACTGTACTAGACGTCCTTGATCGTTATTCGG";
					} else if( args[f+offset].endsWith("expression-interval.txt.gz" ) ) {
						li[f] = "0.0\t0.0\t0.0";
					} else if( args[f+offset].endsWith("specific.txt.gz") ) {
						li[f]="";
						int anzahl = args[f+offset].substring(args[f+offset].lastIndexOf("/")+1).startsWith("qN-") ? 2 : 1;
						for( int an = 0; an < anzahl; an++ ) {
							li[f]+=(an==0?"":"\t") + "0.0";
						}
					} else if( args[f+offset].endsWith("conserved.txt.gz") ) {
						li[f]="";
						int anzahl = args[f+offset].substring(args[f+offset].lastIndexOf("/")+1).startsWith("dnase") ? 2 : 1;
						for( int an = 0; an < anzahl; an++ ) {
							li[f]+=(an==0?"":"\t") + "0.0\t0.0";
						}
					} else if( args[f+offset].endsWith("_bw.gz") ) {
						li[f] = "0.0\t0.0\t0.0\t0.0\t0.0\t0.0";
					} else if( args[f+offset].endsWith("DNase-peak2interval.txt.gz") ) {
						li[f] = "0.0\t0.0";
					} else if( args[f+offset].endsWith("fc.signal.bigwig-interval-orange.txt.gz") ) {
						li[f] = "0.0\t0.0\t0.0\t0.0\t0.0";
					} else if( args[f+offset].endsWith("DNase-peakStat2interval.txt.gz") ) {
						li[f] = "0.0\t0.0\t0.0\t0.0";
					} else if( args[f+offset].endsWith(".nearest_tss.txt.gz") ){
						li[f] = chr+"\t"+start+"\t1E6\t-1\t1E6\t-1\t1E6\t-1\t1E6\t-1";
					}
					//PROBLEM
					else {
						System.out.println("WARNING: no default values for " + args[f+offset] );
					}
				} else {
					li[f] = line[f];
				}
				any |= read[f];
			}
			if( !any ){
				if( use ) System.out.println("last window\t" + chr + "\t" + p + "\t" + (p+50) + "\t" + t.getElapsedTime() );
				
				System.out.println(Arrays.toString(line));
				//synchronize files
				boolean changed = false;
				for(int f=0;f<l;f++){
					//System.out.println(f + "\t" + args[f+offset] + "\t" + line[f]);
					int n=0;
					while( line[f] != null && (line[f].charAt(0) != '[' || line[f].indexOf('_')>=0) ) {
						line[f] = reader[f].readLine();
						n++;
					}
					changed |= n>0;
					//if( n > 0 ) System.out.println("->\t" + args[f+offset] + "\t" + line[f] + "\t" + n);
					if( line[f] == null ) {
						break outerloop;
					}
					if( !line[f].equals(line[0]) ) {
						throw new Exception( "Mismatch:"
								+ "\nFile 0: " + args[offset] + "\nLine: " + line[0]
								+ "\nFile "+f+": " + args[f+offset] + "\nLine: " + line[f]
						);
					}
					read[f] = true;
				}
				if( changed )
				{
					System.out.println( "=> " + Arrays.toString(line));
				}
				
				chr = line[0].substring(1, line[0].length()-1);
				use = chrs!=null ? chrs.contains(chr) : true;
				if( !use ) {
					System.out.println("skip " + chr);
				}
				
				c=0;
				start = 0;
				p=-100000;
			} else {
				if( use ) {
					//if( c == 0 ) System.out.println( "values: " + Arrays.toString(line));
					
					//clear uninformative parts at the beginning
					if( remove == null ) {
						remove = new int[l];
						for(int f=0;f<l;f++){
							remove[f] = 0;
							if( li[f].startsWith(chr+"\t") ) {
								int pos = li[f].indexOf('\t')+1;
								li[f]=li[f].substring(pos);
								remove[f]++;
								if( li[f].startsWith(start+"\t") ) {
									pos = li[f].indexOf('\t')+1;
									li[f]=li[f].substring(pos);
									remove[f]++;
								}
							}
						}
					} else {
						for(int f=0;f<l;f++){
							for( int r = 0; r < remove[f]; r++ ) {
								int pos = li[f].indexOf('\t')+1;
								li[f]=li[f].substring(pos);
							}
						}
					}
					
					c++;
					start+=50;
					p=(start-(windows+1)/2*50);
					
					//determine split length
					if( split == null ) {
						//create empty array
						a = 0;
						for( int f = 0; f < li.length; f++ ) {
							a += li[f].split("\t").length;
						}
						split = new String[windows*a];
						last = a*(windows-1);
					} else {
						//move entries left
						System.arraycopy(split, a, split, 0, last);
						//TODO Think about not moving/copying to accelerate the program
					}
					
					//insert
					for(int f=0, b=last;f<l;f++){
						String[] help = li[f].split("\t");
						System.arraycopy(help, 0, split, b, help.length);
						b+=help.length;
					}
					if( c >= windows ) {
						//create initial sequence
						ArbitrarySequence s = pars.parse(origCon, split);
						if( origCon == null ) {
							origCon = s.getAlphabetContainer();
						}
						
						
						// predict
						w.append( chr + "\t" + p );
						for( int d = 0; d < cl.size(); d++ ){
							AbstractScoreBasedClassifier cla =cl.get(d);
							scores[0] = cla.getScore(s, 0);
							scores[1] = cla.getScore(s, 1);
							ls = Normalisation.getLogSum(scores);
							w.append( "\t" + Math.exp(scores[0]-ls) );
						}
						w.newLine();/**/
					}/**/
				}
			}
		}
		System.out.println("elapsed time: " + t.getElapsedTime() );
		
		//close all files
		w.close();
		for(int f=0;f<l;f++){
			reader[f].close();
		}
	}
}
