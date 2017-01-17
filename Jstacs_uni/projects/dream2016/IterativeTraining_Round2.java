package projects.dream2016;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.classifiers.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.classifiers.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.classifiers.performanceMeasures.AucPR;
import de.jstacs.classifiers.performanceMeasures.AucROC;
import de.jstacs.classifiers.performanceMeasures.PRCurve;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.sequences.ArbitrarySequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.FileManager;
import de.jstacs.io.StringExtractor;
import de.jstacs.io.XMLParser;
import de.jstacs.results.ResultSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.IndependentProductDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.continuous.ConstantDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.continuous.GaussianNetwork;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.FixedStructure;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.HomogeneousMMDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.MixtureDiffSM;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.Time;
import de.jstacs.utils.ToolBox;
import projects.dream2016.Aggregation_multi2.Aggregate;
import projects.dream2016.mix.NewMixtureClassifier;
import projects.dream2016.mix.NewMixtureClassifier.OptimizableMSPClassifier;
import projects.dream2016.mix.NewMixtureClassifier.Training;
import projects.dream2016.mix.NewMixtureClassifier.Vote;
import projects.dream2016.mix.SubDiffSM;

/**
 * Performs a iterative training
 * 
 * @author Jens Keilwagen
 */
public class IterativeTraining_Round2 {
	
	//This class reflects a region that is used for training
	static class Region implements Comparable<Region> {

		int pos;
		String chr, line;
		int lab;
		String[] split;
		
		Region( String line, int lab ) {
			this.line = line;
			split = line.split("\t");
			chr = split[0];
			pos = (int) Math.round( 
					//Double.parseDouble(split[1]) //start
					//(Double.parseDouble(split[2])-50+Double.parseDouble(split[1]))/2 //middle
					Double.parseDouble(split[2])//end
					);
			this.lab=lab;
		}
		
		@Override
		public int compareTo(Region o) {
			return Integer.compare( pos, o.pos );
		}		
	}
	
	//load regions from file
	private static void fill( String fName, int lab, HashMap<String, ArrayList<Region>> region ) throws IOException {
		File f = new File( fName);
		if( f.exists() ) {
			BufferedReader r= new BufferedReader( new FileReader(fName) );
			String line;
			while( (line=r.readLine()) != null ) {
				//System.out.println(line);
				Region reg = new Region(line, lab);
				ArrayList<Region> c = region.get(reg.chr);
				if( c == null ) {
					c=new ArrayList<Region>();
					region.put(reg.chr, c);
				}
				c.add( reg );
			}
			r.close();
		}
	}
	
	//sort regions
	private static int sort( HashMap<String, ArrayList<Region>> region ) throws IOException {
		int windows=-1;
		Iterator<String> it = region.keySet().iterator();
		while( it.hasNext() ) {
			String chr = it.next();
			ArrayList<Region> list = region.get(chr);
			if(list.size()>0){
				if( windows<0 ) {
					Region r = list.get(0);
					windows = (int) Math.round( (Double.parseDouble(r.split[2])-Double.parseDouble(r.split[1]))/50 );
				}

				Collections.sort(list);
				sos.writeln(chr + "\t" + list.size() + "\t" + list.get(0).pos + "\t" + list.get(list.size()-1).pos );
			}
		}
		return windows;
	}
	
	//create DataParser
	private static DataParser getDataParser( int offset, String[] args, String confOut ) throws IllegalArgumentException, IOException, DoubleSymbolException {
		int col=0;
		LinkedList<String> conf = new LinkedList<String>();
		for(int f=offset;f<args.length;f++){
			if( args[f].endsWith("hg19.fa_other.txt.gz") ) {
				conf.add( "Percent\t" + col + "," + (col+1) + "\tEach\tMeanCenter3" ); // GC-content, CpG
				col+=2;
				conf.add( "Entropy\t" + (col+2) + "\tEach\tMaxCenter3" );// KL3
				col+=3;
			} else if( args[f].endsWith("hg19.fa_tracts.txt.gz") ) {
				conf.add( "Length\t" + col+","+(col+1)+","+(col+2)+","+(col+3) + "\tEach\tMaxCenter3" );//Poly-A, Poly-C, poly-at, poly-cg
				col+=4;
			} else if( args[f].endsWith("gencode.v19.types.txt.gz") ) {
				conf.add( "Region\t" + col + "\tEach\tMin" );//5 Annotationen auf beiden strands
				col+=1;
			} else if( args[f].endsWith(".bigwig-interval.txt.gz") ) {
				// Min Median Max 25%
				conf.add( "Coverage\t" + col + "\tEach\tMinCenter3" );// Min MinCenter 3
				conf.add( "Coverage\t" + (col+2)	+ "," + (col+3)	+ "\tEach\tMin" ); // Max Min und 25% Min
				conf.add( "Coverage\t" + (col+1)	+ "\tEach\tEach" ); // Median Each
				
				col+=4;
			} else if( args[f].endsWith(".xml_winscores2.txt.gz") ) {
				
				//Sum Max 0.25 0.33 ssp
				conf.add( "Score\t" + (col) + "\tEach\tMaxCenter3" );// Sum MaxCenter3
				conf.add( "Score\t" + (col+1) + "\tEach\tCenter" );// Max Center
				conf.add( "Score\t" + col + "\tEach\tLogSum" );// Sum LogSum
				col+=2;
				//conf.add( "Distance\t" + col + "," + (col+1) + "\tEach\tMean" );//0.25 0.33 Mean
				col+=2;
				
				col++;
				
			} else if( args[f].endsWith("_winscores2.txt.gz") ) {
				conf.add("Score\t" + col + "\tEach\tMax"); // Sum Max
			//	conf.add("Score\t" + col + "\tEach\tMin");
			//	conf.add("Score\t" + col + "\tEach\tMinCenter3");
				col += 2;
				col += 2;
				col++;
			} else if( args[f].endsWith("hg19.genome.fa.seqs.gz") ){
				conf.add("Seq\t" + col + "\tEach\tCenter3");
				col++;
			} else if( args[f].endsWith("expression-interval.txt.gz" ) ) {
				conf.add( "Coverage\t" + col + "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+1) + "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+2) + "\tEach\tCenter" );
				col+=3;
			}
			//NEW
			else if( args[f].endsWith("specific.txt.gz") ) {
				int anz = args[f].substring(args[f].lastIndexOf("/")+1).startsWith("qN-") ? 2 : 1;
				for( int a = 0; a < anz; a++ ) {
					conf.add( "Coverage\t" + col	+ "\tEach\tCenter" );
					col+=1;
				}
			} else if( args[f].endsWith("conserved.txt.gz") ) {
				int anz = args[f].substring(args[f].lastIndexOf("/")+1).startsWith("dnase") ? 2 : 1;
				for( int a = 0; a < anz; a++ ) {
					conf.add( "Coverage\t" + col	+ "\tEach\tCenter" );
					conf.add( "Coverage\t" + (col+1)	+ "\tEach\tCenter" );
					col+=2;
				}
			} else if( args[f].endsWith("_bw.gz") ) {
				conf.add( "Coverage\t" + col	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+1)	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+2)	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+3)	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+4)	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+5)	+ "\tEach\tCenter" );
				col+=6;
			} else if( args[f].endsWith("DNase-peak2interval.txt.gz") ) {
				conf.add( "Coverage\t" + col	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+1)	+ "\tEach\tCenter" );
				col+=2;
			} else if( args[f].endsWith("fc.signal.bigwig-interval-orange.txt.gz") ) {
				conf.add( "Coverage\t" + col	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+1)	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+2)	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+3)	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+4)	+ "\tEach\tCenter" );
				col+=5;
			} else if( args[f].endsWith("DNase-peakStat2interval.txt.gz") ) {
				conf.add( "Coverage\t" + col	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+1)	+ "\tEach\tMaxCenter5" );
				conf.add( "Coverage\t" + (col+2)	+ "\tEach\tCenter" );
				conf.add( "Coverage\t" + (col+3)	+ "\tEach\tMaxCenter5" );
				col += 4;
			} else if( args[f].endsWith(".nearest_tss.txt.gz") ) {
				conf.add( "Distance\t"+(col)+","+(col+2)+","+(col+4)+","+(col+6)+"\tMin\tCenter" );
				col += 8;
			} else if( args[f].endsWith("unicov.gz") ){
				conf.add( "Coverage\t" + (col) + "\tEach\tMinCenter3" );
				conf.add( "Coverage\t" + (col+2) + "\tEach\tMaxCenter3" );
				col += 4;
			}
			//EXCEPTION
			else{
				throw new RuntimeException("Unknown file type: " + args[f]);
				//conf.append( "???\tsubsequent columns might be wrong\n" );
			}
		}
		
		conf.addFirst( 5 + "\t" + col );
		/*
		for( int c = 0; c < conf.size(); c++ ) {
			System.out.println( conf.get(c) );
		}/**/
		
		if( confOut != null ) {
			PrintWriter confw = new PrintWriter(confOut+".conf");
			confw.println(0 + "\t" + col);
			for( int c = 1; c < conf.size(); c++ ) {
				confw.println( conf.get(c) );
			}
			confw.close();
		}
	
		return new DataParser(conf.toArray(new String[0]));
	}

	private static Time time;
	private static AlphabetContainer con;
	private static SafeOutputStream sos = SafeOutputStream.getSafeOutputStream(System.out);
	
	/**
	 * Performs a iterative training
	 * 
	 * @param args
	 * 0	= pos coordinates
	 * 1	= neg coordinates
	 * 2	= ambi (+pos) coordinates
	 * 3	= conf prefix
	 * 4	= blacklistfiltered
	 * 5	= labels
	 * 6	= column in labels
	 * 7	= epigram
	 * 8	= threads
	 * 9... = feature files
	 * 
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		time = Time.getTimeInstance(null);
				
		HashMap<String, ArrayList<Region>> region = new HashMap<String, ArrayList<Region>>();
		fill(args[0], 0, region);
		fill(args[1], 1, region);
		//fill(args[2], 2, region);
		
		int epigram = Integer.parseInt(args[7]);
		int threads = Integer.parseInt(args[8]);
		
		int off=9;
		DataParser pars = getDataParser(off, args, args[3]);
		
		//iterative training:		
		ArrayList<String> chrs = new ArrayList<String>();		
		HashMap<String, ArrayList<int[]>> hash = new HashMap<String, ArrayList<int[]>>();
		BufferedReader r = new BufferedReader( new FileReader(args[4]) );
		String line;
		while( (line = r.readLine()) != null ) {
			String[] split = line.split("\t");
			ArrayList<int[]> current = hash.get(split[0]);
			if( current == null ) {
				current = new ArrayList<int[]>();
				hash.put(split[0], current);
				chrs.add(split[0]);
			}
			current.add(new int[]{Integer.parseInt(split[1]), Integer.parseInt(split[2])} );
		}
		r.close();
		System.out.println("Iterative training using: " + chrs);
		System.out.println();
		
		AbstractScoreBasedClassifier cl=null;
		ArrayList<AbstractScoreBasedClassifier> list = new ArrayList<AbstractScoreBasedClassifier>();
		int stop=10;
		int col = Integer.parseInt(args[6]);
		for( int idx = 0; idx < stop; idx++ ) {
			System.out.println("ITERATION: " + idx);
			//sort
			int windows = sort(region);
			//extract
			extract( off, idx, pars, args, region, windows, null, idx == 0 ? null : chrs );//XXX f�r Jan muss hier, das stehen: idx == 0 ? null : chrs ); //schneller Test: chrs ); //XXX
			if( cl == null ) {
				cl = create("ipsf",false,epigram,threads);
			} else {
				cl = cl.clone();
			}
			//train
			list.add(cl);
			train( args[3]+"_positives"+idx+".txt", args[3]+"_negatives"+idx+".txt", idx, list );
			
			//TODO if we like to modify the positives
			/*
			if( idx == 0 ) {
				boolean f = false;
				Iterator<String> it = chrs.iterator();
				while( it.hasNext() ) {
					String chr = it.next();
					ArrayList<Region> list = region.get(chr);
					for( int i = 0; i < list.size(); ) {
						Region re = list.get(i);
						if( re.fg ) {
							if( f )  {
								System.out.println( re.pos +"\t" +re.line );
							}
							list.remove(i);
						} else {
							i++;
						}
					}
					f=false;
				}
			}/**/
			
			//System.out.println(cl);
			//gws
			extract( off, idx, pars, args, null, windows, cl, chrs );
			
			//delete regions
			Iterator<String> it = region.keySet().iterator();
			while( it.hasNext() ) {
				region.get( it.next() ).clear();
			}
			
			//agg, eval and expand
			int add = evaluate( idx, hash, chrs,  args[3]+".gws.gz", args[5], col, 3, region, windows, false );
			
			if( add == 0 ) {
				break;
			}
		}/**/
	}
	
	static int newNeg=-1;
	
	//extract data or perform gws
	private static void extract( int offset, int idx, DataParser pars, String[] args, HashMap<String, ArrayList<Region>> region, int windows, AbstractScoreBasedClassifier cl, ArrayList<String> chrs ) throws Exception {
		time.reset();

		String outpath = args[3];

		int l = args.length-offset;

		//open all needed files
		BufferedReader[] reader = new BufferedReader[l];
		for(int f=0;f<l;f++){			
			GZIPInputStream stream = new GZIPInputStream(new FileInputStream(args[f+offset]));
			reader[f] = new BufferedReader(new InputStreamReader(stream));
		}
		
		//parse data
		int c=0, start = 0, last=-1, a=-1, p=-1000;
		String chr = null;		
		
		String[] line = new String[l], li = new String[l];
		int[] remove = null;
		String[] split = null;
		double[] scores = new double[2];
		double ls = 0;
		boolean any;
		boolean[] read = new boolean[l];
		Arrays.fill( read, true );
		BufferedWriter[] w;
		int x=0;
		if( cl == null ) {
			w = new BufferedWriter[4+(idx==0?2:0)];
			x=w.length/2;
			w[0] = new BufferedWriter( new FileWriter( outpath+"_positives"+idx+".txt" ) );
			w[1] = new BufferedWriter( new FileWriter( outpath+"_negatives"+idx+".txt" ) );
			w[x] = new BufferedWriter( new FileWriter( outpath+"_positives"+idx+".txt.weights" ) );
			w[x+1] = new BufferedWriter( new FileWriter( outpath+"_negatives"+idx+".txt.weights" ) );
			
			/*if( idx==0 ) {
				w[2] = new BufferedWriter( new FileWriter( outpath+"_ambi.txt" ) );
				w[x+2] = new BufferedWriter( new FileWriter( outpath+"_ambi.txt.weights" ) );
			}*/			
		} else {
			w = null;
		}
		BufferedWriter gws = 
				//without zip cl != null ? new BufferedWriter( new FileWriter(  outpath+"_positives"+idx+".gws" ) ) : null;
				cl != null ? new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream( new FileOutputStream(outpath+".gws.gz"+idx) ), "UTF-8") ) : null;
		ArrayList<Region> reg = null;
		int r = 0;
		boolean use=false;
		int anz = chrs!=null ? chrs.size() : Integer.MAX_VALUE;
		int skip=0, end=Integer.MAX_VALUE;
		outerloop: while( true ){
			any = false; 
			for(int f=0;f<l;f++){
				if( read[f] ) {
					try{
						line[f] = reader[f].readLine();
					}catch(OutOfMemoryError err){
						System.out.println(f+" "+idx);
						System.out.println(args[f+offset]);
						System.out.println(outpath);
						System.out.flush();
						throw err;
					}
					read[f] = line[f]!=null && line[f].charAt(0) != '[';
				}
				
				any |= read[f];
			}
			if( !any ){
				if( use ) {
					System.out.println("last window\t" + chr + "\t" + p + "\t" + (p+50) + "\t" + skip + "\t" + time.getElapsedTime() + "\t" + anz );
				}
				if( anz == 0 ) {
					break outerloop;
				}
				
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
				reg = region == null ? null : region.get(chr);
				r = 0;
				if( reg != null && r < reg.size() ) {
					end = reg.get(r).pos;
				} else {
					end = Integer.MAX_VALUE;
				}
				use = chrs == null || chrs.contains(chr);
				if( !use ) {
					System.out.println("skip " + chr);
				} else {
					anz--;
				}
				c=0;
				start = 0;
				p=-100000;
			} else {
				if( use ) {
					for(int f=0;f<l;f++){
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
							} else if( args[f+offset].endsWith("unicov.gz") ){
								li[f] = "0\t0\t0\t0";
							}
							//PROBLEM
							else {
								System.out.println("WARNING: no default values for " + args[f+offset] );
							}
							
						} else {
							li[f] = line[f];
						}
					}
					
					//if( c == 0 ) System.out.println( "values: " + Arrays.toString(line));

					int st = start;
					start+=50;
					c++;
					p=(start-(windows+1)/2*50);
					
					if( cl != null || start+windows*50 >= end ) 
					{
						//clear uninformative parts at the beginning
						if( remove == null ) {
							remove = new int[l];
							for(int f=0;f<l;f++){
								remove[f] = 0;
								if( li[f].startsWith(chr+"\t") ) {
									int pos = li[f].indexOf('\t')+1;
									li[f]=li[f].substring(pos);
									remove[f]++;
									if( li[f].startsWith(st+"\t") ) {
										pos = li[f].indexOf('\t')+1;
										li[f]=li[f].substring(pos);
										remove[f]++;
									}
								}
							}
						} else {
							for(int f=0;f<l;f++){
								for( int re = 0; re < remove[f]; re++ ) {
									int pos = li[f].indexOf('\t')+1;
									li[f]=li[f].substring(pos);
								}
							}
						}
	
						
						
						boolean out=false;
						Region current=null;
						if( reg != null && r < reg.size() && start == reg.get(r).pos ) {
							//System.out.println("JETZT");
							//System.out.println(chr + "\t" + start);
							current = reg.get(r);
							//System.out.println(current.chr + "\t" + current.start);
							out = true;						
							/*
							r++;
							boolean set=false;
							double wei = Double.parseDouble(current.split[4]);
							while( r < reg.size() && reg.get(r).fg == current.fg && reg.get(r).start==current.start ) {
								current.split[3] = "m";
								wei += Double.parseDouble( reg.get(r).split[4] );
								reg.remove(r);
								set=true;
							}
							if( set ) {
								current.split[4] = ""+wei;
							}
							*/
						}
						
						
						//determine split length
						if( split == null ) {
							//create empty array
							a = 0;
							for( int f = 0; f < li.length; f++ ) {
								a += li[f].split("\t").length;
							}
							split = new String[5+windows*a];
							last = 5+a*(windows-1);
						} else {
							//move entries left
							System.arraycopy(split, 5+a, split, 5, last-5);
							//TODO Think about not moving/copying to accelerate the program
						}
						
						//insert
						for(int f=0, b=last;f<l;f++){
							String[] help = li[f].split("\t");
							System.arraycopy(help, 0, split, b, help.length);
							b+=help.length;
						}
						
						//System.out.println(Arrays.toString(split));
						
						if( (out || cl != null ) && c >= windows ) {
							//create initial sequence
							if( current != null ) {
								System.arraycopy(current.split, 0, split, 0, 5);
								//System.out.println(Arrays.toString(split));
								//System.out.println();
							}
							ArbitrarySequence s = null;
							try{
								s = pars.parse(con, split);
							}catch(WrongAlphabetException exc){
								System.out.println(Arrays.toString(split));
								System.out.println(con);
								System.out.println(pars.getAlphabets(split.length));
								throw exc;
							}
							if( con == null ) {
								con = s.getAlphabetContainer();
								FileManager.writeFile(outpath+"_positives.txt.alpha", con.toXML());
							}
							
							Arrays.fill(split, 0,5,null);
							
							if( cl == null ) {
								//extract
								do {
									current = reg.get(r);
								
									w[current.lab].append(s.toString("\t", 0, s.getLength()));
									w[current.lab].newLine();
									w[current.lab+x].append(current.split[4]);
									w[current.lab+x].newLine();
									
									r++;
								} while( r < reg.size() && /*reg.get(r).lab == current.lab &&*/ reg.get(r).pos==current.pos );
								
								if( reg != null && r < reg.size() ) {
									end = reg.get(r).pos;
								} else {
									end = Integer.MAX_VALUE;
								}
							} else {
								// predict
								gws.append( chr + "\t" + p );
								//for( int d = 0; d < cla.size(); d++ ){
								//	AbstractScoreBasedClassifier cl =cla.get(d);
									scores[0] = cl.getScore(s, 0);
									scores[1] = cl.getScore(s, 1);
									ls = Normalisation.getLogSum(scores);
									gws.append( "\t" + Math.exp(scores[0]-ls) );
								//}
								gws.newLine();/**/
							}/**/
						}
					} else {
						skip++;
					}
				}
			}
		}
		
		System.out.println( (cl==null?"extract":"gws") + " - elapsed time: " + time.getElapsedTime()+ " "+skip );
		System.out.println();
		
		//close all files
		if( w!= null ) {
			for(int i=0;i<w.length;i++){
				if( w[i] != null ) {
					w[i].close();
				}
			}
		}
		if( gws != null ) {
			gws.close();
		}
		for(int f=0;f<l;f++){
			reader[f].close();
		}
	}
	
	//create an classifier
	private static AbstractScoreBasedClassifier create( String type, boolean useDeps, int epigram, int threads ) throws Exception {
		double[] beta = LearningPrinciple.getBeta(LearningPrinciple.MCL);
		double ess = 1;
		
		ArrayList<DifferentiableStatisticalModel> fun = new ArrayList<DifferentiableStatisticalModel>();
		IntList lens = new IntList();
		
		//DNase
		int d = 0;
		while( d < con.getPossibleLength() && !con.isDiscreteAt(d) ) {
			d++;
		}
		int[][] structure = new int[d][];
		if(useDeps){
			structure[1] = new int[0];
			structure[2] = new int[]{1};
			structure[0] = new int[]{2};	
			int c = 3 + (d-3)/2;
			structure[c] = new int[]{0};
			for( int h = 0; c+h+1<d; h++ ) {
				structure[c+h+1] = new int[]{c+h};
				structure[c-h-1] = new int[]{c-h};
			}	
		}else{
			structure = new int[d][0];
		}
		
		fun.add( new GaussianNetwork(con.getSubContainer(0, d), structure) );
		lens.add(d);
		
		//Region
		structure = new int[10][];
		if(useDeps){
			structure[4] = new int[0];
			structure[9] = new int[]{4};
			
			structure[3] = new int[]{4};
			structure[8] = new int[]{9};
			
			structure[2] = new int[]{3};
			structure[7] = new int[]{8};
			
			structure[1] = new int[]{3};
			structure[6] = new int[]{8};
			
			structure[0] = new int[]{2};
			structure[5] = new int[]{7};
		}else{
			structure = new int[10][0];
		}
		fun.add( new BayesianNetworkDiffSM(con.getSubContainer(d, 10), 10, 256.0/*c*/, true, new FixedStructure(structure)) );
		lens.add(10);
		
		//Global
		structure = new int[7][];
		if(useDeps){
			structure[0] = new int[0];
			structure[1] = new int[]{0};
			structure[2] = new int[0];
			structure[5] = new int[]{0};
			structure[3] = new int[]{5};
			structure[6] = new int[]{0};
			structure[4] = new int[]{6};
		}else{
			structure = new int[7][0];
		}
		fun.add( new GaussianNetwork(con.getSubContainer(lens.get(0)+lens.get(1), 7), structure) );
		lens.add(7);
		
		//Epigram
		if( epigram > 0 ) {
			structure = new int[epigram][0]; //independent
			fun.add( new GaussianNetwork(structure) );
			lens.add(epigram);
		}
		
		//Motifs
		int le = 0;
		for( int i = 0; i < lens.length(); i++ ) {
			le += lens.get(i);
		}

		structure = new int[3][];
		if(useDeps){
			structure[2] = new int[0];
			structure[1] = new int[]{2};
			structure[0] = new int[]{1};
		}else{
			structure = new int[3][0];
		}
		int l = con.getPossibleLength();
		while( le < l && !con.isDiscreteAt(le) ) {
			fun.add( new GaussianNetwork(con.getSubContainer(le, 3), structure) );
			lens.add(3);
			le+=3;
		}
		
		//sequence
		if(le < l){
			//fun.add( new MarkovModelDiffSM(con.getSubContainer(le, 50), 50, 4.0, true, new InhomogeneousMarkov(3)) );
			fun.add( new HomogeneousMMDiffSM(con.getSubContainer(le, 150), 3, 256.0, HomogeneousMMDiffSM.getSumOfHyperParameters(3, 150, 256.0), true, true, 1) );
			lens.add(150);
			le += 150;
		}
		
		//expression & new features
		int anz = l-le;
		if( anz > 0 ) {
			fun.add( new GaussianNetwork(con.getSubContainer(le, anz), new int[anz][0]) );
			lens.add(anz);
			le += anz;
		}
		
		//classifier
		//System.out.println(le + "\t" + l );
				
		AbstractScoreBasedClassifier cl;
		GenDisMixClassifierParameterSet ps = new GenDisMixClassifierParameterSet(con, l, Optimizer.QUASI_NEWTON_BFGS, 1E-6, 1E-6, 1, false, KindOfParameter.ZEROS, true, threads); //TODO LAST?
		LogPrior prior = 
				//new SimpleGaussianSumLogPrior(100); 
				DoesNothingLogPrior.defaultInstance;
		
		int start = 0;
		
		GenDisMixClassifier gdsm;
		OptimizableMSPClassifier[] comp;
		DifferentiableStatisticalModel[] mix;
		
		type = type.toLowerCase();
		IndependentProductDiffSM ipsf = new IndependentProductDiffSM(ess, true, fun.toArray(new DifferentiableStatisticalModel[0]),lens.toArray());
		switch( type ) {
			case "ipsf":
				gdsm = new GenDisMixClassifier( ps, prior, beta, ipsf, ipsf );
				//gdsm.setOutputStream(null);
				cl = gdsm;
				break;
			case "bgmix":
				MixtureDiffSM bgMix = new MixtureDiffSM(3, true, ipsf,ipsf);
				gdsm = new GenDisMixClassifier( ps, prior, beta, ipsf, bgMix );
				//gdsm.setOutputStream(null);
				cl = gdsm;
				break;
			case "mix_each":
				comp = new OptimizableMSPClassifier[fun.size()];
				mix = new DifferentiableStatisticalModel[comp.length];
				start=0;
				for( int i = 0; i < comp.length; i++ ) {
					SubDiffSM s = new SubDiffSM( con, l, fun.get(i), start);
					comp[i] = new OptimizableMSPClassifier( ps, prior, s, s );
					comp[i].setOutputStream(null);
					mix[i] = new ConstantDiffSM(con,l);
					
					start+=lens.get(i);
				}							
				cl = new NewMixtureClassifier( threads, Training.COMBINED, 3, mix, comp, Vote.VOC, prior );
				break;
			case "mix_dnase_single":
			case "mix_dnase_block":
				mix = ArrayHandler.createArrayOf(
						type.endsWith("single") 
							? new SubDiffSM( con, l, new GaussianNetwork(new int[1][0] ), 0)
							: new SubDiffSM( con, l, fun.get(0), 0)
						, 2);
				 
				OptimizableMSPClassifier o = new OptimizableMSPClassifier( ps, prior, ipsf, ipsf );
				//o.setOutputStream(null);
				OptimizableMSPClassifier[] cla = ArrayHandler.createArrayOf(o, mix.length);
				cl = new NewMixtureClassifier( threads, Training.COMBINED, 3, mix, cla, Vote.VOC, prior.getNewInstance() );
				break;
			default:
				throw new IllegalArgumentException("Unkown type: " + type );
		}
		return cl;
	}
		
	private static DataSet[] data = new DataSet[2];
	private static double[][] weights;
	
	private static long[][] hist = new long[2][101];
	
	//load data and train classifier
	private static void train( String pFile, String nFile, int idx, ArrayList<AbstractScoreBasedClassifier> list ) throws Exception {
		time.reset();
		
		if( idx ==  0 ) {
			data = new DataSet[]{
					new DataSet(con,new StringExtractor(new File(pFile), 1000, '#'),"\t"),
					new DataSet(con,new StringExtractor(new File(nFile), 1000, '#'),"\t")
				};
			weights = new double[][] {
					DataParser.getWeights(pFile+".weights", DataParser.Weighting.ONE ),
					DataParser.getWeights(nFile+".weights", DataParser.Weighting.DIRECT ),
				};			
			newNeg = (int) Math.round(ToolBox.sum(weights[1])*0.15);//TODO
		} else {
			DataSet part = new DataSet(con,new StringExtractor(new File(nFile), 1000, '#'),"\t");	
			double[] pWeights = DataParser.getWeights(nFile+".weights", DataParser.Weighting.DIRECT );
			//data[1] =part;
			//weights[1] = pWeights;
			
			data[1] = DataSet.union(part,data[1]);
			
			int i = 1;
			double[] help = new double[pWeights.length+weights[i].length];
			System.arraycopy(pWeights, 0, help, 0, pWeights.length);
			System.arraycopy(weights[i], 0, help, pWeights.length, weights[i].length);
			weights[i] = help;
			
			/*
			if( idx == 1 ) {
				data[0] = new DataSet(con,new StringExtractor(new File(pFile), 1000, '#'),"\t");
				weights[0] = new double[data[0].getNumberOfElements()];
			}*/
			double[] scores = new double[2];
			for( int c=0; c<1/*weights.length*/; c++ ) {
				Arrays.fill(hist[0], 0);
				for( i = 0; i < weights[c].length; i++ ) {
					Sequence s = data[c].getElementAt(i);
					double w = 0;
					if( a == Aggregate.MeanProd ) {
						for( int e = 0; e < list.size()-1; e++ ) {
							AbstractScoreBasedClassifier cl = list.get(e);
							scores[c] = cl.getScore(s, 0);
							scores[1] = cl.getScore(s, 1);
							double ls = Normalisation.getLogSum(scores);
							w += Math.exp(scores[0]-ls);
						}
						w /= list.size()-1;
					} else {
						throw new RuntimeException(""+a);
					}
					hist[0][(int)Math.round(w*100)]++;
					weights[c][i]=w;//TODO set this to 1?

					//System.out.println( i + "\t" + weights[c][i] );
				}/**/
				System.out.println(Arrays.toString(hist[0]));
			}
			//}*/
		}
		
		for( int i = 0; i < data.length; i++ ) {
			System.out.println(i + ": #=" + data[i].getNumberOfElements() + ", length=" + data[i].getElementLength() + ", " + data[i].getAnnotation() );
		}
		
		AbstractScoreBasedClassifier cl = list.get(list.size()-1);
		
		cl.train(data,weights);

		System.out.println("train - elapsed time: " + time.getElapsedTime() );
		System.out.println();
		
		if( idx >= 0 ) {
			StringBuffer xml = new StringBuffer();
			XMLParser.appendObjectWithTags(xml, cl, "classifier");
			FileManager.writeFile( pFile + "-classifier.xml", xml);

		}
	}

	//evaluate performance and add data
	private static int evaluate( int idx, HashMap<String, ArrayList<int[]>> hash, ArrayList<String> chrs, String gwsFileName, String truth, int col, int wi, HashMap<String, ArrayList<Region>> region, int windows, boolean curve ) throws IOException {
		//determine threshold
		double th = agg(idx,hash,chrs, gwsFileName, truth,col, wi, Double.NaN, newNeg, null, 0, curve );
		System.out.println("threshold: "+th);
		//add negatives
		return (int) agg(idx, hash,chrs, gwsFileName, truth,col, wi, th, 0, region, windows, curve );
	}
	
	static String type = "EACH";
	static Aggregate a = Aggregate.MeanProd;
	
	//aggregate predictions and add data
	private static double agg( int index, HashMap<String, ArrayList<int[]>> hash, ArrayList<String> chrs, String gwsFileName, String truth, int col, int wi, double th, int anz2, HashMap<String, ArrayList<Region>> region, int bins, boolean curve ) throws NumberFormatException, IOException {		
		//without zip BufferedReader reader = new BufferedReader( new FileReader(gwsFileName) );
		BufferedReader[] reader = new BufferedReader[index+1];
		for( int i = 0; i <= index; i++ ) {
			reader[i] = new BufferedReader( new InputStreamReader(new GZIPInputStream(new FileInputStream(gwsFileName+i))));
		}
		String[] line = new String[reader.length];
		
		//read prediction
		int windows = 2 * wi;	
		int offset = (wi-2)*50;

		String chr=null;
		ArrayList<int[]> intervals = null;
		int[] inter = null;
		HashMap<String, ArrayList<double[][]>> predictions = new HashMap<String, ArrayList<double[][]>>();
		ArrayList<double[][]> currentPred = null;
		double[][] p = null;
		double[] values = new double[reader.length];
		int idx = -1, h = -1;
		int anz = type.equals("EACH") ? reader.length : 1;
		while( (line[0] = reader[0].readLine()) != null ) {
			for( int i = 1; i <= index; i++ ) {
				line[i] = reader[i].readLine();
			}
			if( line[0].charAt(0)=='#' ) {
				continue;
			}
			String[] split = line[0].split("\t");
			
			if( chr == null || !chr.equals(split[0]) ) {
				/*if( chr != null && intervals!= null ) {
					out(idx, chr, inter, h, p, windows);
				}*/
				
				chr = split[0];
				intervals = hash.get(chr);
				if( intervals != null  ) {
					//System.out.println( chr + "\t" + toString( intervals ) );
					currentPred = new ArrayList<double[][]>();
					predictions.put( chr, currentPred );
					idx = 0;
					inter = intervals.get(idx);
					
					if( anz>0 ) {
						p = new double[anz][(inter[1]-inter[0])/50+1+2*(wi-2)];
						currentPred.add( p ) ;
					}
					h=0;
				} else {
					inter = null;
					idx = -1;
					p = null;
				}
			}
			if( intervals != null ) {
				int start = Integer.parseInt( split[1] );
				if( start < inter[0] - offset ) {
					//skip
				} else if( start <= inter[1] + offset ) {
					for( int i = 0; i < values.length; i++ ) {
						values[i] = Double.parseDouble(line[i].split("\t")[2]);
					}
					
					
					switch ( type ) {
						case "MEAN": p[0][h] = ToolBox.mean(0, values.length, values); break;
						case "MAX": p[0][h] = ToolBox.max(0, values.length, values); break;
						case "MIN": p[0][h] = ToolBox.min(0, values.length, values); break;
						case "PROD": 
							p[0][h] = 1;
							for( int k = 0; k <values.length; k++ ) {
								p[0][h] *= (1-values[k]);
							}
							p[0][h] = 1 - p[0][h];
							break;
						case "EACH": 
							for( int k = 0; k <values.length; k++ ) {
								p[k][h] = values[k];
							}
							break;
						default: throw new IllegalArgumentException( "unknown:" + type );
					}
					
					h++;
				} else {
					if( idx+1 < intervals.size() ) {
						//out(idx, chr, inter, h, p, windows);

						idx++;
						if( idx  < intervals.size() ) {
							inter = intervals.get(idx);
							p = new double[anz][(inter[1]-inter[0])/50+1+2*(wi-2)];
							currentPred.add( p ) ;
							h=0;
						} else {
							h=-1;
						}
					}
				}
			}
		}
		for( int i = 0; i <= index; i++ ) {
			reader[i].close();
		}
		//out(idx, chr, inter, h, p, windows);
		/**/
		
		DoubleList pos = new DoubleList();
		DoubleList neg = new DoubleList();

		GZIPInputStream stream = new GZIPInputStream(new FileInputStream(truth));
		BufferedReader r = new BufferedReader(new InputStreamReader(stream));
		String first = r.readLine();
		if( Double.isNaN(th) ) System.out.println( "Compute performance for cell type: " + first.split("\t")[col] );
		int x = 50*(bins-1)/2;
		int[] num=new int[2];
		ArrayList<Region> list=null;
		if( region == null ) {
			Arrays.fill(hist[0], 0);
			Arrays.fill(hist[1], 0);
		}
		double[] weight = new double[p.length];
		Arrays.fill(weight, 1d/p.length);
		for( int c = 0; c < chrs.size(); c++ ) {
			chr = chrs.get(c);
			currentPred = predictions.get(chr);
			if( region != null ) {
				list = region.get(chr);
			}
			if( currentPred != null ) {
				intervals = hash.get(chr);
				for( int i = 0; i < currentPred.size(); i++ ) {
					inter = intervals.get(i);
					p = currentPred.get(i);
					int start = inter[0];
					for( int j = 0; j < p[0].length-windows; j++ ) {
						double ag = -1;
						switch( a ) {
							case Max: ag = ToolBox.max(j, j+windows, p[0]); break;
							case Mean: ag = ToolBox.mean(j, j+windows, p[0]); break;
							case Prod: 
								ag = 1;
								for( int k = j; k <j+windows; k++ ) {
									ag *= (1-p[0][k]);
								}
								ag = 1 - ag;
								break;
							case Median: ag = ToolBox.median(j, j+windows, p[0]); break;
							case MeanProd:
								ag = 0;
								for( int e = 0; e < p.length; e++ ) {
									double ag2 = 1;
									for( int k = j; k <j+windows; k++ ) {
										ag2 *= (1-p[e][k]);
									}
									ag += weight[e]*(1 - ag2);
								}
								break;
							/*case ProdProd:
								ag = 1;
								for( int e = 0; e < p.length; e++ ) {
									double ag2 = 1;
									for( int k = j; k <j+windows; k++ ) {
										ag2 *= (1-p[e][k]);
									}
									ag *= 1 - ag2;
								}
								ag = Math.pow(ag, 1d/p.length);
								break;*/
							default:
								throw new IllegalArgumentException("not implemented: " +a);
						}
						//o.writeln(chr + "\t" + start +"\t" + (start+200) + "\t" + ag );
						
						if( r != null ) {
							String s = r.readLine();
							String[] split = s.split("\t");
							if( j == 0 ) {
								//System.out.println(chr + "\t"+start+"\t"+s);
								
								if( !(chr.equals(split[0]) && split[1].equals(""+start)) ) {
									System.out.println();
									System.out.println( s );
									System.out.println( Arrays.toString(inter) );
									System.out.println( chr + "\t" + start +"\t" + (start+200) );
									System.exit(1);
								}
							}
							if( region == null ) {
								switch( split[col].charAt(0) ) {
									case 'b': case 'B': pos.add(ag); hist[0][(int)Math.round(100*ag)]++; break;
									case 'u': case 'U': neg.add(ag); hist[1][(int)Math.round(100*ag)]++; break;
								}
							} else {								
								if( ag>=th && (split[col].charAt(0)=='U' || split[col].charAt(0)=='u') ) {
									list.add( new Region(chr+"\t"+(start-x) + "\t"+(start+x+50) +"\tn\t1", 1) );
									num[1]++;
								}
							}
						}
						start+=50;
					}
				}
			} else {
				System.out.println("WARNING: Did not find predictions for " + chr);
			}
		}
		r.close();

		if( region == null ) {
			System.out.println( "pos=c"+Arrays.toString(hist[0]).replace('[', '(').replace(']', ')') );
			System.out.println( "neg=c"+Arrays.toString(hist[1]).replace('[', '(').replace(']', ')') );
		}
		if( Double.isNaN(th) ) {
			double[] po = pos.toArray();
			double[] ne = neg.toArray();
			Arrays.sort(po);
			Arrays.sort(ne);
			ResultSet aucpr = new AucPR().compute( po, ne );

			System.out.println( "#positives: " + po.length + "\t" + po[0]+ " .. " + po[po.length-1]);
			System.out.println( "#negatives: " + ne.length + "\t" + ne[0]+ " .. " + ne[ne.length-1] );
				
			System.out.println( "random: " + po.length / (double)(po.length + ne.length) );
			
			System.out.println();
			System.out.println(windows + "\t" + a  + "\t" + Arrays.toString(weight) );
			System.out.println( new AucROC().compute( po, ne ) );
			System.out.println( aucpr );
			
			//recall @ X% fdr
			ResultSet rs = new PRCurve().compute(po, ne);
			DoubleTableResult dtr = (DoubleTableResult)rs.getResultAt(rs.getNumberOfResults()-1);
			double[][] c = dtr.getValue();
			double[] t = {0.95,0.9,0.75,0.5};
			double[] max = new double[t.length];
			Arrays.fill(max, 0);
			BufferedWriter w = null;
			if( curve ) {
				w = new BufferedWriter(new FileWriter(gwsFileName + index +".prcurve"));
			}
			for( int i = 0; i < c.length; i++ ) {
				for( int j = 0; j < t.length; j++ ) {
					if( c[i][1] >= t[j] && c[i][0] > max[j] ) {
						max[j] = c[i][0];
					}
				}
				if( w != null ) {
					w.append( c[i][0] + "\t" + c[i][1] );
					w.newLine();
				}
			}
			if( w != null ) {
				w.close();
			}

			for( int j = 0; j < t.length; j++ ) {
				System.out.println( "Recall at " + Math.round(100*(1-t[j])) + "% FDR: " + max[j] );
			}

			//return Math.max( ne[Math.max(ne.length-anz,0)], po[0] );
			return Math.max( ne[Math.max(ne.length-anz2,0)], po[(int)Math.floor(po.length*0.01)] );//TODO
		} else {
			System.out.println("add " + Arrays.toString(num) );
			System.out.println( list.size() );
			return num[0]+num[1];
		}
	}
}
