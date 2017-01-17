package projects.dream2016;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

import de.jstacs.classifiers.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.classifiers.performanceMeasures.AucPR;
import de.jstacs.classifiers.performanceMeasures.AucROC;
import de.jstacs.classifiers.performanceMeasures.PRCurve;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.ToolBox.TiedRanks;

public class Aggregation_multi {

	static enum Aggregate {
		Mean,
		Median,
		Max,
		Prod,
		MeanProd,
		MeanProd2,
		Test,
		Test2,
		ProdProd;
	}
	
	/**
	 * @param args
	 * 0 file
	 * 1 .. type or type:index1,index2,...
	 * 2 .. window length/2
	 * 3 .. list
	 * 4 .. type
	 * optional:
	 * 5 .. labels
	 * 6 .. column of cell type
	 */
	public static void main(String[] args) throws IOException {
		System.out.println("Bugfix3");
		System.out.println(args[0]);
		//determine files:
		//BufferedReader reader = new BufferedReader( new FileReader(args[0]) );
		GZIPInputStream stream2 = new GZIPInputStream(new FileInputStream(args[0]));
		BufferedReader reader = new BufferedReader(new InputStreamReader(stream2));
		int[] index = null;
		int idx2 = args[1].indexOf(':');
		String type = args[1];
		if( idx2 >= 0 ) {
			String[] split = type.substring(idx2+1).split(",");
			index = new int[split.length];
			for( int i = 0; i<split.length; i++ ) {
				index[i] = Integer.parseInt( split[i] );
			}
			type = type.substring(0, idx2);
		}
		type = type.toUpperCase();
		System.out.println(type + "\t" + Arrays.toString(index));
		
		int m = 2;
		
		//read list of intervals
		
		//order of chr: "Important Note: You have to make sure that the coordinates in the submitted predictions appear in the same order as in the template, otherwise the predictions will not be scored"
		ArrayList<String> chrs = new ArrayList<String>();		
		HashMap<String, ArrayList<int[]>> hash = new HashMap<String, ArrayList<int[]>>();
		BufferedReader r = new BufferedReader( new FileReader(args[m+1]) );
		System.out.println(args[m+1]);
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
		System.out.println(chrs);
		System.out.println();
				
		//read prediction
		int wi = Integer.parseInt(args[m]);
		int windows = 2 * wi;	
		int offset = (wi-2)*50;
		
		String chr=null;
		ArrayList<int[]> intervals = null;
		int[] inter = null;
		HashMap<String, ArrayList<double[][]>> predictions = new HashMap<String, ArrayList<double[][]>>();
		ArrayList<double[][]> currentPred = null;
		double[][] p = null;
		double[] values = null;
		int idx = -1, h = -1, anz = type.equals("EACH") ? -1 : 1;
		String[][] split = new String[Math.max(1, wi-2)][];
		int b = 0;
		boolean firstInter=false;
		while( (line = reader.readLine()) != null ) {
			if( line.charAt(0)=='#' ) {
				continue;
			}
			split[b] = line.split("\t");
			
			if( chr == null || !chr.equals(split[b][0]) ) {
				if( chr != null && intervals!= null ) {
					out(idx, chr, inter, h, p, windows);
				}
				
				chr = split[b][0];
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
				firstInter=true;
			}
			if( intervals != null ) {
				int start = Integer.parseInt( split[b][1] );
				if( start < inter[0] - offset ) {
					//skip
				} else if( start <= inter[1] + offset ) {
					if( values == null ) {
						if( index == null ) {
							index = new int[split[0].length-2];
							for( int i = 2; i < split[0].length; i++ ) {
								index[i-2] = i;
							}
							System.out.println(type + "\t" + Arrays.toString(index));
						}
						if( anz < 0 ) {
							anz=index.length;
							p = new double[anz][(inter[1]-inter[0])/50+1+2*(wi-2)];
							currentPred.add( p ) ;
						}
						values = new double[index.length];
					}
					if( firstInter ) {
						//adjust h for large w
						if( start!= inter[0]-offset ) {
							h = (start - (inter[0]-offset))/50;
							System.out.println(Arrays.toString(inter) + "\t" + start + "\t" + (inter[0]-offset) + "\tset h=" + h );
						}
						firstInter = false;
					}
					
					fill(type, p, index, values, h, split, b);
					h++;
				} else {
					if( idx+1 < intervals.size() ) {
						out(idx, chr, inter, h, p, windows);

						idx++;
						if( idx  < intervals.size() ) {
							inter = intervals.get(idx);
							p = new double[anz][(inter[1]-inter[0])/50+1+2*(wi-2)];
							currentPred.add( p ) ;
							
							if( start < inter[0]-offset ) {
								h=0;
							} else {
								h = (start-(inter[0]-offset))/50+1;							
								//TODO fill p[*][v] for v in [0,h-1] using backup
								for( int c=b, hh = h-1; hh >= 0; hh-- ) {
									fill(type, p, index, values, hh, split, c);
									c--;
									if( c < 0 ) {
										c = split.length-1;
									}
								}
							}
						} else {
							h=-1;
						}
					}
				}
			}
			b++;
			if( b == split.length ) {
				b=0;
			}
		}
		reader.close();
		out(idx, chr, inter, h, p, windows);
		/**/

		Aggregate a = Aggregate.valueOf( args[m+2] );
		//aggregate & write
		SafeOutputStream o = SafeOutputStream.getSafeOutputStream( args.length == m+3 ? new FileOutputStream( args[0]+ "-" + args[m] + "-" + a + ".tab") : null );
		DoubleList pos = new DoubleList();
		DoubleList neg = new DoubleList();
		
		int col;
		if( args.length==m+5 ) {
			col=Integer.parseInt(args[m+4]);
			try {
				GZIPInputStream stream = new GZIPInputStream(new FileInputStream(args[m+3]));
				r = new BufferedReader(new InputStreamReader(stream));
			} catch ( ZipException ze ) {
				r = new BufferedReader(new FileReader(args[m+3]));
			}
			String first = r.readLine();
			System.out.println( "Compute performance for cell type: " + first.split("\t")[col] );
		} else {
			r = null;
			col=-1;
		}
		
		long[][][] hist = new long[1+(a==Aggregate.MeanProd || a==Aggregate.MeanProd2?index.length:0)][2][101];
		for( int c = 0; c < hist.length; c++) {
			Arrays.fill(hist[c][0], 0);
			Arrays.fill(hist[c][1], 0);
		}
		
		for( int c = 0; c < chrs.size(); c++ ) {
			chr = chrs.get(c);
			currentPred = predictions.get(chr);
			if( currentPred != null ) {
				intervals = hash.get(chr);
				for( int i = 0; i < currentPred.size(); i++ ) {
					inter = intervals.get(i);
					p = currentPred.get(i);
					int start = inter[0];
					for( int j = 0; j < p[0].length-windows; j++ ) {
						String s;
						String[] sp;
						DoubleList dl=null;
						int clazz=-1;
						if( r != null ) {
							s = r.readLine();
							sp = s.split("\t");
							if( j == 0 ) {
								//System.out.println(chr + "\t"+start+"\t"+s);
								
								if( !(chr.equals(sp[0]) && sp[1].equals(""+start)) ) {
									System.out.println();
									System.out.println( s );
									System.out.println( Arrays.toString(inter) );
									System.out.println( chr + "\t" + start +"\t" + (start+200) );
									System.exit(1);
								}
							}
							
							switch( sp[col].charAt(0) ) {
								case 'b': case 'B': dl=pos; clazz=0; break;
								case 'u': case 'U': dl=neg; clazz=1; break;
							}
						}
						
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
							case MeanProd:
								ag = 0;
								for( int e = 0; e < p.length; e++ ) {
									double ag2 = 1;
									for( int k = j; k <j+windows; k++ ) {
										ag2 *= (1-p[e][k]);
									}
									ag +=  1 - ag2;
									
									if( clazz >= 0 ) {
										hist[e+1][clazz][(int)Math.round(100*(1-ag2))]++;
									}
								}
								ag /= p.length;
								break;
							case MeanProd2:
								ag = 0;
								for( int e = 0; e < p.length; e++ ) {
									double ag2 = 1;
									for( int k = j; k <j+windows; k++ ) {
										ag2 *= (1-p[e][k]);
									}
									ag +=  Math.pow(1 - ag2,1d/windows);
									
									if( clazz >= 0 ) {
										hist[e+1][clazz][(int)Math.round(100*(1-ag2))]++;
									}
								}
								ag /= p.length;
								break;
							case ProdProd:
								ag = 1;
								for( int e = 0; e < p.length; e++ ) {
									double ag2 = 1;
									for( int k = j; k <j+windows; k++ ) {
										ag2 *= (1-p[e][k]);
									}
									ag *=  1 - ag2;
								}
								break;
							case Test:
								ag = 0;
								for( int e = 0; e < p.length; e++ ) {
									for( int w = 0; w<=wi-2; w++) {
										double ag2 = 1;
										for( int k = j+w; k <j+windows-w; k++ ) {
											ag2 *= (1-p[e][k]);
										}
										ag += Math.pow(1 - ag2,1d/(windows-2d*w));
									}
								}
								ag /= (wi-1)*p.length;
								break;
							case Test2:
								ag = 0;
								for( int e = 0; e < p.length; e++ ) {
									double v = 1;
									for( int w = 0; w<=wi-2; w++) {
										double ag2 = 1;
										for( int k = j+w; k <j+windows-w; k++ ) {
											ag2 *= (1-p[e][k]);
										}
										 v*= Math.pow(1 - ag2,1d/(windows-2d*w));
									}
									
									ag +=Math.pow(v, 1d/(wi-1d));
								}
								ag /= p.length;
								break;
							case Median: ag = ToolBox.median(j, j+windows, p[0]); break;
							default:
								throw new IllegalArgumentException("not implemented: " +a);
						}
						o.writeln(chr + "\t" + start +"\t" + (start+200) + "\t" + ag );
						if( r != null ) {
							if( clazz >= 0 ) {
								dl.add(ag);
								hist[0][clazz][(int)Math.round(100*ag)]++;
							}
						}
						
						start+=50;
					}
				}
			} else {
				System.out.println("WARNING: Did not find predictions for " + chr);
			}
		}
		o.close();
		if( r != null ) {
			r.close();
			
			for( int c = 0; c<hist.length; c++ ) {
				System.out.println( c==0?"combined":"" );
				System.out.println( "pos"+(c==0?"":(c-1))+"=c"+Arrays.toString(hist[c][0]).replace('[', '(').replace(']', ')') );
				System.out.println( "neg"+(c==0?"":(c-1))+"=c"+Arrays.toString(hist[c][1]).replace('[', '(').replace(']', ')') );
			}
		}
		
		if( col>-1 ) {
			double[] po = pos.toArray();
			double[] ne = neg.toArray();
			Arrays.sort(po);
			Arrays.sort(ne);
			
			System.out.println( "#positives: " + po.length + "\t" + po[0]+ " .. " + po[po.length-1]);
			System.out.println( "#negatives: " + ne.length + "\t" + ne[0]+ " .. " + ne[ne.length-1] );
			System.out.println( "random: " + po.length / (double)(po.length + ne.length) );
			
			System.out.println();
			System.out.println(windows + "\t" + args[1] + "\t" + args[m+2] );
			System.out.println( new AucROC().compute( po, ne ) );
			System.out.println( new AucPR().compute( po, ne ) );
			
			//recall @ X% fdr
			ResultSet rs = new PRCurve().compute(po, ne);
			DoubleTableResult dtr = (DoubleTableResult)rs.getResultAt(rs.getNumberOfResults()-1);
			double[][] c = dtr.getValue();
			double[] t = {0.95,0.9,0.75,0.5};
			double[] max = new double[t.length];
			Arrays.fill(max, 0);
			for( int i = 0; i < c.length; i++ ) {
				for( int j = 0; j < t.length; j++ ) {
					if( c[i][1] >= t[j] && c[i][0] > max[j] ) {
						max[j] = c[i][0];
					}
				}
			}
			for( int j = 0; j < t.length; j++ ) {
				System.out.println( "Recall at " + Math.round(100*(1-t[j])) + "% FDR: " + max[j] );
			}
		}
	}
	
	static void fill( String type, double[][] p, int[] index, double[] values, int h, String[][] split, int b ) {
		for( int i = 0; i < index.length; i++ ) {
			values[i] = Double.parseDouble(split[b][index[i]]);
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
	}
	
	static void out( int idx, String chr, int[] interval, int last, double[][] p, int windows ) {
		if( chr != null && interval!= null && p != null && p[0].length != last ) {
			System.out.println(chr + "\t" + idx + "\t" + Arrays.toString(interval) + "\t" + last + "\t" + p[0].length + "\t" + (last!= p[0].length?"WARNING":"") );
		}
	}
	
	static String toString( ArrayList<int[]> list ) {
		StringBuffer sb = new StringBuffer();
		for( int i = 0; i < list.size(); i++ ) {
			sb.append( (sb.length()==0?"":", ") + Arrays.toString( list.get(i) ) );
		}
		return sb.toString();
	}

	static String toString( int start, int end, double[] p ) {
		StringBuffer sb = new StringBuffer();
		while( start < end )  {
			sb.append( (sb.length()==0?"":", ") + p[start++] );
		}
		return sb.toString();
	}
	
	static class ArrayComparator implements Comparator<int[]> {

		@Override
		public int compare(int[] o1, int[] o2) {
			return o1[0] - o2[0];
		}
	}
}