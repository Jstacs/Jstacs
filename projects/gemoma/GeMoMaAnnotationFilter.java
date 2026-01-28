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

package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import javax.script.ScriptException;

import org.graalvm.polyglot.Context;

import de.jstacs.DataType;
import de.jstacs.clustering.kmeans.Kmeans;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.FileExistsValidator;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.parameters.validation.RegExpValidator;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.SafeOutputStream;

/**
 * This class allows to filter GeMoMa annotation files (GFFs) to obtain the most relevant predictions.
 * 
 * @author Jens Keilwagen
 * 
 * @see GeMoMa
 */
public class GeMoMaAnnotationFilter extends GeMoMaModule {
	
	private static String replace( String cds, String oldString, String newString, int idx ) {
		if( (idx=cds.indexOf(oldString,idx))>0 ) {
			char c = cds.charAt(idx-1);
			if( c=='\t' || c==';' ) {
				cds = cds.substring(0, idx)
						+ newString
						+ cds.substring(idx+oldString.length());
			} else {
				cds = replace( cds, oldString, newString, idx+1 );
			}
		}
		return cds;
	}
	
	private static int MAX;
	
	private static HashMap<String,int[]> del = new HashMap<String,int[]>();
	static {
		del.put("Parent",new int[1]);
		del.put("alternative",new int[1]);
		del.put("evidence",new int[1]);
		del.put("sumWeight",new int[1]);
	}
	
	static class UserSpecifiedComparator implements Comparator<Prediction>{
		static HashSet<String> asDouble = new HashSet<String>();
		static {
			asDouble.add("aa");
			asDouble.add("raa");
			asDouble.add("score");
			asDouble.add("bestScore");
			asDouble.add("maxScore");
			
			asDouble.add("tae");
			asDouble.add("tde");
			asDouble.add("tie");
			asDouble.add("minSplitReads");
			asDouble.add("tpc");
			asDouble.add("minCov");
			asDouble.add("avgCov");
			
			asDouble.add("iAA");
			asDouble.add("pAA");
			asDouble.add("lpm");
			asDouble.add("maxGap");
			
			asDouble.add("evidence");
			asDouble.add("sumWeight");
		}
		
		String[] keys;
		HashSet<String> errors;
		Protocol protocol;
		Context context;
		
		public UserSpecifiedComparator( String specified, Context context, Protocol protocol ) {
			keys=specified.split(",");
			for( int i = 0; i < keys.length; i++ ) {
				keys[i] = keys[i].trim();
			}
			errors = new HashSet<String>();
			this.context = context;
			this.protocol = protocol;
		}
		
		public static Double parse( String d ) {
			if( d.equals( "NA" ) ) {
				return Double.NaN;
			} else {
				return -Double.parseDouble(d);
			}
		}
		
		public void prepare( Prediction p ) {
			p.comp = new Comparable[keys.length+1];
			for( int i = 0; i < keys.length; i++ ) {
				String v = p.hash.get(keys[i]);
				boolean asDoubleNow = asDouble.contains(keys[i]);
				if( v == null ) {
					try {
						v=Tools.eval(context, keys[i], p.hash);
						asDoubleNow=true;
					} catch( ScriptException se ) {
						if( !errors.contains(keys[i])) {
							errors.add(keys[i]);
							protocol.appendWarning("Could not find key (" + keys[i] + ") that is used for sorting\n");
						}
					}
				}
				if( v!=null && asDoubleNow ) {
					p.comp[i] = parse(v);
				} else {
					p.comp[i] = v;
				}
			}
			p.comp[keys.length]=p.id;
		}

		@SuppressWarnings("unchecked")
		public int compare(Prediction o1, Prediction o2) {
			int d=0;
			for( int i=0; i < o1.comp.length && d==0; i++ ) {
				if( o1.comp[i]!=null && o2.comp[i]!=null ) {
					d = o1.comp[i].compareTo( o2.comp[i] );
				}
			}
			return d;
		}
	}
	
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String tempD) throws Exception {
		gene=h=0;
		String tag = parameters.getParameterForName("tag").getValue().toString();
		int maxTranscripts = (Integer) parameters.getParameterForName("maximal number of transcripts per gene").getValue();
		double cbTh = (Double) parameters.getParameterForName("common border filter").getValue();
		String[] defAttributes = parameters.getParameterForName("default attributes").getValue().toString().split(",");
		for( int i = 0; i < defAttributes.length; i++ ) {
			defAttributes[i] = defAttributes[i].trim();
			if( defAttributes[i].length() == 0 ) {
				defAttributes[i] = null;
			}
		}
		
		Context context = Context.newBuilder("js").build();
		

		String filter = Tools.prepareFilter( (String) parameters.getParameterForName("filter").getValue() );
		UserSpecifiedComparator comp = new UserSpecifiedComparator( (String) parameters.getParameterForName("sorting").getValue(), context, protocol );
		Double d = null;
		Parameter pa = parameters.getParameterForName("length difference");
		if( pa != null ) {
			d = (Double) pa.getValue();
		}
		double lenPerc = d==null ? Double.POSITIVE_INFINITY : d;

		String altFilter = Tools.prepareFilter( (String) parameters.getParameterForName("alternative transcript filter").getValue() );
		boolean addAltTransIDs = (Boolean) getParameter( parameters, "add alternative transcripts" ).getValue();
		boolean addAdd = (Boolean) getParameter( parameters, "transfer features" ).getValue();
		
		pa = getParameter( parameters, "kmeans" );
		int minNum=0, cluster=0, good = 0;
		String[] ex = null;
		int margin=-1;
		double quantile=-1;
		
		ParameterSet ps = (ParameterSet) pa.getValue();
		if( ps.getNumberOfParameters()>0 ) {
			minNum = (Integer) ps.getParameterForName("minimal number of predictions").getValue();
			pa = ps.getParameterForName("cluster");
			if( pa != null ) {
				cluster = (Integer) pa.getValue();
			}
			pa = ps.getParameterForName("good cluster");
			if( pa != null ) {
				good = (Integer) pa.getValue();
			}
			if( good >= cluster ) {
				throw new IllegalArgumentException("The number of good clusters must be smaller than the total number of clusters");
			}
			
			/*
			ExpandableParameterSet att = (ExpandableParameterSet) ps.getParameterForName("cluster attribute").getValue();
			ex = new String[att.getNumberOfParameters()];
			for( int j = 0; j < ex.length; j++ ) {
				ex[j] = ((SimpleParameterSet) att.getParameterAt(j).getValue()).getParameterAt(0).getValue().toString();
			}
			*/
			ex = new String[]{
					//"Math.abs(maxScore-Math.max(0,score))/maxScore"
					//to be tested: 
					"Math.max(0,1 - Math.max(0,score)/maxScore)"
			};
			
			ps = (ParameterSet) ps.getParameterForName("trend").getValue();
			if( ps.getNumberOfParameters()>0 ) {
				//local
				margin = (Integer) ps.getParameterForName("margin").getValue();
				quantile = (Double) ps.getParameterForName("quantile").getValue();
			}
		}
		
		boolean kmeans = cluster>0;
			    
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterForName("predicted annotation").getValue();
		
		String line, t;
		String[] split;
				
		HashMap<String,ArrayList<Prediction>> pHash = new HashMap<String,ArrayList<Prediction>>();
		ArrayList<Prediction> pred;
		Prediction current = null;

		//read genes in GeMoMa gff format!
		BufferedReader r;
		MAX = eps.getNumberOfParameters();
		stats = new int[MAX][5];
		evidenceStat = new int[MAX][3];
		combinedEvidence = new boolean[MAX];
		String[] prefix = new String[MAX];
		double[] weight = new double[MAX];
		ArrayList<String> allInfos = new ArrayList<String>();
		HashSet<String> hash = new HashSet<String>();
		HashSet<String> ids = new HashSet<String>();
		HashMap<String,String[]> annotInfo = new HashMap<String,String[]>();
		boolean rnaSeq=false;
		StringBuffer sb = new StringBuffer();
		String[] att = del.keySet().toArray(new String[0]);
		for( int k = 0; k < MAX; k++ ) {
			SimpleParameterSet sps = ((SimpleParameterSet)eps.getParameterAt(k).getValue());
			String h = (String) sps.getParameterAt(0).getValue();
			prefix[k] = h==null?"":h;
			String sampleInfo = k + (prefix[k]==null||prefix[k].length()==0?"":(" (" + prefix[k] + ")"));
			weight[k] = (Double) sps.getParameterAt(1).getValue();
			protocol.append("species " + sampleInfo + "\n");
					
			annotInfo.clear();
			int transcript = -1, go = -1, defline=-1;
			if( sps.getParameterAt(3).isSet() ) {
				//read annotation information
				r = new BufferedReader(new FileReader(sps.getParameterAt(3).getValue().toString()));
				String[] header = r.readLine().split("\t");
				for( int i = 0; i < header.length; i++ ) {
					switch( header[i] ) {
						case "transcriptName": transcript=i; break;
						case "GO": go=i; break;
						default: 
							if( header[i].endsWith("-defline") ) {
								defline=i;
							}
					}
				}
				if( transcript<0 || go<0 || defline<0 ) {
					r.close();
					throw new IllegalArgumentException("annotation info " + sampleInfo + " must be a tab-delimited file with at least the following columns: transcriptName, GO, and .*defline");
				}
				while( (line = r.readLine()) != null ) {
					split = line.split("\t");
					annotInfo.put(split[transcript].toUpperCase(), split);
				}
				r.close();
			}
			
			//gff needs to be "clustered" (all features of one transcript should be in one block)
			String fName = sps.getParameterAt(2).getValue().toString();
			//System.out.println(fName);
			hash.clear();
			ArrayList<Prediction> list = new ArrayList<Prediction>();
			r = new BufferedReader(new FileReader(fName));
			while( (line=r.readLine()) != null ) {
				if( line.length() == 0 ) continue;
				if( line.charAt(0)=='#' ) {
					if( line.startsWith(INFO) ) {
						hash.add(line);
					}
					continue;
				}
				
				split = line.split("\t");
				
				t = split[2];
				if( t.equals("gene") ) {
					//TODO WARNING
				} else if( t.equals(tag) ) { //tag: prediction/mRNA/transcript
					current = new Prediction(split, MAX, k, prefix[k], del, annotInfo, go, defline, defAttributes);
					if( ids.contains( current.id ) ) {
						//System.out.println(line);
						//System.out.println(current);
						r.close();
						throw new IllegalArgumentException("The id (" + current.id + ") has been used before. You can try to use a prefix (cf. parameters)." );
					} else {
						ids.add(current.id);
					}
					list.add(current);
					rnaSeq |= split[8].indexOf(";tie=")>0;
				} else {
					if( current==null ) {
						r.close();
						throw new NullPointerException("There is no gene model. Please check parameter \"tag\" and the order within your annotation file "+sampleInfo+": " + fName );
					}
					if( split[8].contains("Parent="+current.oldId) ) {
						if( t.equals( "CDS" ) ) current.addCDS( line );
						else if( addAdd ) current.addAdd( line );
					} else {
						//TODO WARNING
						r.close();
						throw new IllegalArgumentException("The GFF has to be clustered, i.e., all features of a transcript must be adjacent lines:\n" + line);
					}
				}
			}
			r.close();
			if( hash.size() != 1 ) {
				if( hash.size() == 0  ) {
					protocol.appendWarning("Evidence file " + sampleInfo + " contains no parameter description\n" );
				} else {
					protocol.appendWarning("Evidence file " + sampleInfo + " contains a different parameter description\n" );
				}
			}
			sb.delete(0,sb.length());
			for( int i =0; i < att.length; i++ ) {
				int[] v = del.get(att[i]);
				if( v[0] > 0 ) {
					sb.append( (sb.length()==0?"":", ") + v[0] + " \"" + att[i] + "\"" );
					//delete values => prepare next round
					v[0]=0;
				}
			}
			if( sb.length()>0 ) {				
				protocol.appendWarning("Delete attributes in evidence file " + sampleInfo + ": " + sb  + "\n");
			}
			Iterator<String> it = hash.iterator();
			while( it.hasNext() ) {
				allInfos.add(it.next());
			}
			
			//ML-filter of predictions
			if( kmeans && list.size() >= minNum ) {
				Collections.sort(list);
				double[][] matrix = new double[list.size()][ex.length];
				boolean anyNAN = false;
				for( int i = 0; i < list.size(); i++ ) {
					Prediction cur = list.get(i);
					for( int j = 0; j< ex.length; j++ ) {
						matrix[i][j] = Double.parseDouble(Tools.eval(context, ex[j], cur.hash));
						anyNAN |= Double.isNaN( matrix[i][j] ) ;
					}
				}
				
				if( matrix.length>1 && !anyNAN ) {
					if( margin>0 ) {	
						double[] quant = new double[list.size()];
						
						//compute quantile using margin and quantile
						int first=0, last=0, old=-1;
						String oldChr="";
						for( int i = 0; i < list.size(); i++ ) {
							Prediction cur = list.get(i), help;
							boolean changed=false;
							if( !cur.split[0].equals(oldChr) || old != cur.start ) {
								//find start
								while( !(help=list.get(first)).split[0].equals(cur.split[0]) 
										|| help.start+margin<cur.start ) {
									first++;
									changed=true;
								}
								//find end
								while( last < list.size() && (help=list.get(last)).split[0].equals(cur.split[0]) 
										&& help.start<cur.start+margin ) {
									last++;
									changed=true;
								}
							}
							
							if( !changed ) {
								//same predictions lead to same value 
								quant[i] = quant[i-1];
							} else {
								if( last-first==1 ) {
									//only one prediction detrend does not make sense
									quant[i]=0;
								} else {
									//find quantile between first and last
									
									//naive implementation: TODO improve
									double[] forSorting = new double[last-first];
									for( int j = first; j < last; j++ ) {
										forSorting[j-first] = matrix[j][0];
									}
									Arrays.sort(forSorting);
									
									quant[i] = forSorting[ (int) Math.ceil((forSorting.length-1)*quantile) ];
								}
							}
//System.out.println(first + ".." + i + ".."+last + "\t" + quant[i] );
							
							oldChr=cur.split[0];
							old=cur.start;
						}
						
						//use local trend
						for( int i = 0; i < list.size(); i++ ) {
							matrix[i][0] -= quant[i];
						}
					}
					
					double[][] data = ex.length>1 ? Kmeans.rescale(matrix) : matrix;
					int[] assignment = new Kmeans().cluster(data, cluster, 10);
					double[][] mean = new double[cluster][ex.length];
					int[] stat = new int[cluster];
					for( int i = 0; i < matrix.length; i++ ) {
						for( int j = 0; j< ex.length; j++ ) {
							mean[assignment[i]][j] += matrix[i][j];
						}
						stat[assignment[i]]++;
					}
					protocol.append("cluster assignment: " + Arrays.toString(stat) + "\n" );
					double[] dist = new double[mean.length];
					for( int i = 0; i < mean.length; i++ ) {
						mean[i][0]/=stat[i];
						dist[i]=mean[i][0]; /*0
						for( int j = 0; j< ex.length; j++ ) {
							mean[i][j] /=stat[i];
							dist[i]+=mean[i][j]*mean[i][j];
						}*/
					}
					double[] clone = dist.clone();
					Arrays.sort(clone);
					double threshold = clone[good];
					boolean[] use = new boolean[dist.length];
					for( int i = 0; i < dist.length; i++ ) {
						use[i] = dist[i]<threshold;
						protocol.append( i + /*"\t" + Arrays.toString(mean[i]) +*/ "\t" + dist[i] + "\t" + use[i] + "\t" + stat[i] + "\n" );
					}
					
					ArrayList<Prediction> mlFiltered = new ArrayList<Prediction>();
					for( int i = 0; i < list.size(); i++ ) {
						Prediction cur = list.get(i);
						if( use[assignment[i]] ) {
							mlFiltered.add(cur);
						}
					}
					protocol.append( "ml-filtered: " + mlFiltered.size() + "/" + list.size() + "\n" );
					
					list.clear();
					list=mlFiltered;
				} else {
					protocol.append( "ml-filter not possible\n" );
				}
			}
			
			for( Prediction p: list ) {
				pred = pHash.get( p.contigStrand );
				if( pred == null ) {
					pred = new ArrayList<Prediction>();
					pHash.put(p.contigStrand,pred);
				}
				pred.add( p );
			}
		}
		
		Set<String> s = pHash.keySet();
		int anz = 0;
		for( String chrSt: s ) {
			anz += pHash.get(chrSt).size();
		}
		protocol.append( "all: " + anz + "\n" );

		int filtered = 0;
		if( anz>0 ) {
			//per contig/chromosome and strand combination => could be parallelized
			for( String chrSt: s ) {
				pred = pHash.get(chrSt);
				
				//combine redundant
				Collections.sort(pred);
				Prediction last = pred.get(pred.size()-1);
				for( int i = pred.size()-2; i >= 0; i-- ){
					Prediction cur = pred.get(i);
					if( cur.compareTo(last) == 0 ) {
						//identical
						if( cur.score > last.score || (cur.score== last.score && cur.id.compareTo(last.id)<0) ) {
							cur.combine(last, addAltTransIDs);
							pred.remove(i+1);
							last=cur;
						} else {
							last.combine(cur, addAltTransIDs);
							pred.remove(i);
						}
					} else {
						last=cur;
					}
				}
				
				//filter: user-specified
				for( int i = pred.size()-1; i >= 0; i-- ){
					Prediction p = pred.get(i);
					p.setEvidenceAndWeight( weight );
					if( !Tools.filter(context, filter, p.hash) ) {
						//System.out.println(p.hash.toString());
						pred.remove(i);
					} else {
						fillEvidence(p.evidence, 0);
/*
String[] array = {"aa","raa","score","tie","tpc","pAA","iAA","lpm","maxGap"};
System.out.print(p.id);
for( int z=0; z < array.length; z++ ) {
	System.out.print("\t" + p.hash.get(array[z]));
}
System.out.println();
*/
					}
				}
				filtered += pred.size();
			}
		}
		protocol.append( "filtered: " + filtered + "\n" );

		File out = Tools.createTempFile("GAF-filtered", tempD);
		pa = parameters.getParameterForName("intermediate result");
		File out2 = (pa==null || !((Boolean)pa.getValue()) ) ? null : Tools.createTempFile("GAF-combined", tempD);
		if( filtered >=0 ) {
			//per contig/chromosome and strand combination => could be parallelized
			BufferedWriter w = new BufferedWriter( new FileWriter(out) );
			SafeOutputStream sos = SafeOutputStream.getSafeOutputStream(out2==null ? null : new FileOutputStream(out2) );
			w.append("##gff-version 3"); w.newLine();
			sos.writeln("##gff-version 3");
			for( int i = 0; i < allInfos.size(); i++ ) {
				String st = allInfos.get(i);
				w.append( st );
				w.newLine();
				sos.writeln(st);
			}
			w.append(INFO + getShortName() + " " + getToolVersion() + "; ");
			String info = JstacsTool.getSimpleParameterInfo(parameters);
			if( info != null ) {
				w.append("SIMPLE PARAMETERS: " + info );
			}
			w.newLine();
			sos.writeln(INFO + getShortName() + " " + getToolVersion() + "; " + (info != null ? ("SIMPLE PARAMETERS: " + info) : "") );
			StringBuffer anno = new StringBuffer();
			for( String chrSt: s ) {
				pred = pHash.get(chrSt);
				if( pred.size() > 0 ) {
					if( !sos.doesNothing() ) {
						for( Prediction p : pred ) {
							anno.delete(0, anno.length());
							p.write(anno, null, gene);
							sos.write(anno);
						}
					}
					//cluster predictions
					split(addAdd, clustered, pred, comp, context, lenPerc, altFilter, rnaSeq, maxTranscripts, w, cbTh, protocol);
				}
			}
			w.close();
			sos.close();
		}
	
		protocol.append( "clustered: " + clustered[0] + "\n\n" );
		protocol.append( "genes: " + gene + "\n" );
		protocol.append( "transcripts: " + transcripts + "\n\n" );
		
		protocol.append( "\tgenes\tgenes with maxTie=1\ttranscripts\ttranscripts with tie=1\ttranscripts with tie=NA, tpc=1\n");
		for( int i  = 0; i < stats.length; i++ ) {
			protocol.append("(max)evidence=" + (i+1));
			for( int j  = 0; j < stats[i].length; j++ ) {
				protocol.append("\t" + stats[i][j]);
			}
			protocol.append("\n");
		}
		
		if( MAX > 1 ) {
			protocol.append("\ninput\tfiltered transcripts\tfinal transcripts supported\tfinal genes supported\n" );
			for( int i = 0; i < MAX; i++ ) {
				protocol.append( "input " + i + (prefix[i].length()>0?" ("+prefix[i]+")":"") + "\t" + evidenceStat[i][0] + "\t" + evidenceStat[i][1]  + "\t" + evidenceStat[i][2] + "\n" );
			}
		}
		ArrayList<Result> res = new ArrayList<Result>();
		res.add( new TextResult("filtered predictions", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "gff", getToolName(), null, true) );
		if( out2!= null ) res.add( new TextResult("intermediate predictions", "Result", new FileParameter.FileRepresentation(out2.getAbsolutePath()), "gff", getToolName(), null, true) );
		return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
	}	
	
	int h;
	int transcripts=0;
	int[] clustered = new int[1];
	
	static int[][] stats;
	int gene;
	
	void split( boolean addAdd, int[] clustered, ArrayList<Prediction> pred, UserSpecifiedComparator comp, Context context, double lenPerc, String altFilter, boolean rnaSeq, int maxTranscripts, BufferedWriter w, double cbTh, Protocol protocol ) throws Exception {
		int i = 0;
		Prediction next = pred.get(i), current;
		ArrayList<Prediction> currentCluster = new ArrayList<Prediction>();
		while( i < pred.size() ) {
			//determine next cluster
			currentCluster.clear();
			current = next;
			currentCluster.add(current);			
			int end = current.end;
			i++;	
			while( i < pred.size() && end>(next = pred.get(i)).start ) {
				end = Math.max(end,next.end);
				currentCluster.add(next);
				i++;
			}
			if( clustered!= null ) clustered[0]++;
			
			//sort and alternative filter
			for( Prediction p: currentCluster ) {
				comp.prepare(p);
				p.altCand = Tools.filter(context, altFilter, p.hash);
			}
			Collections.sort( currentCluster, comp );
			for( Prediction p: currentCluster ) {
				p.comp=null;
			}
			
			//select predictions
			create( addAdd, rnaSeq, maxTranscripts, currentCluster, w, lenPerc, cbTh, protocol );
		}
	}
	
	void create( boolean addAdd, boolean rnaSeq, int maxTranscripts, ArrayList<Prediction> list, BufferedWriter w, double lenPerc, double cbTh, Protocol protocol ) throws Exception {
		ArrayList<ArrayList<Prediction>> sel = new ArrayList<ArrayList<Prediction>>();
		ArrayList<BitSet> usedBp = new ArrayList<BitSet>();
		int clusterStart = Integer.MAX_VALUE;
		int clusterEnd = 0;
		for( int i = 0; i < list.size(); i++ ) {
			Prediction p = list.get(i);
			if( p.start < clusterStart ) {
				clusterStart = p.start;
			}
			if( p.end > clusterEnd ) {
				clusterEnd = p.end;
			}
		}
		
		//XXX change
		int minOverlapThreshold = 0; // can later be extended to small overlap
		for( int i = 0; i < list.size(); i++ ) {
			Prediction n = list.get(i);
			int overlapIdx = -1, enoughOverlap = 0, maxOverlap=0;
			for( int j = 0; j < sel.size(); j++ ) {
				//check overlap
				int o = n.overlap(usedBp.get(j),clusterStart,clusterEnd);
				if( o > minOverlapThreshold ) {
					enoughOverlap++;
					if( o>maxOverlap ) {
						maxOverlap=o;
						overlapIdx=j;
					}
				}
			}
			if( enoughOverlap==0 ) {
				//new gene
				ArrayList<Prediction> selection = new ArrayList<Prediction>();
				selection.add(n);
				n.used=true;
				sel.add( selection );
				
				BitSet usedCov = new BitSet(clusterEnd-clusterStart+1);
				n.fillNuc(usedCov,clusterStart,clusterEnd);
				usedBp.add(usedCov);
			} else {
				//new alternative transcript?
				if( n.altCand && enoughOverlap==1 ) {
					ArrayList<Prediction> selection = sel.get(overlapIdx);
					
					//TODO
					//for( int k = 0; k < selection.size(); k++ ) { //test all members
					int k=0; //test cluster representative
						Prediction x = selection.get(k);
						double diff = Math.abs( x.length - n.length ) / ((double)x.length)*100;
						if( diff <= lenPerc ) {
							double min=2d*Math.min(x.cds.size(),n.cds.size());
							double overlap = x.commonBorders(n)/min;
							//double overlap = n.phaseCompatible(x)/(double) n.length;
							if( overlap>=cbTh ) {
								selection.add(n);
								n.used=true;
								n.fillNuc(usedBp.get(overlapIdx),clusterStart,clusterEnd);
								//break;
							}
						} else {
//System.out.println(x.id + "\t" + x.length + "\t" + n.id + "\t" + n.length + "\t" + diff);
						}
					//}
				}
			}
		}			
		
		//write
		Prediction n=null;
		StringBuffer anno = new StringBuffer();
		for( int j = 0; j < sel.size(); j++ ) {
			ArrayList<Prediction> selection = sel.get(j);
			int start=Integer.MAX_VALUE, end=Integer.MIN_VALUE;
			maxEvidence=0;
			Arrays.fill( combinedEvidence, false );
			maxTie=-1;
			complete=0;
			int i = 0;
			//write transcripts and determine start and end
			anno.delete(0, anno.length());
			for( ; i < selection.size() && i < maxTranscripts; i++ ) {
				n = selection.get(i);
				n.write(anno, protocol, gene);
				
				start = Math.min(start, addAdd ? Integer.parseInt(n.split[3]) : n.start);
				end = Math.max(end, addAdd ? Integer.parseInt(n.split[4]) : n.end );
			}
			transcripts+=i;
			if( i>0 ) {
				//write gene
				fillEvidence(combinedEvidence, 2);
				w.append(n.split[0] + "\tGAF\tgene\t" + start + "\t" + end  + "\t.\t" + n.split[6] + "\t.\tID=gene_"+gene+";transcripts=" + i + ";complete="+complete+ (rnaSeq?";maxTie=" + (maxTie<0?"NA":maxTie):"") +(MAX>1?";maxEvidence="+maxEvidence+";combinedEvidence=" + count(combinedEvidence):"") );
				w.newLine();
				w.append(anno);
				
				stats[maxEvidence-1][0]++;
				if( maxTie == 1 ) {
					stats[maxEvidence-1][1]++;
				}
				gene++;
			}
		}
	}
	
	static int maxEvidence, complete;
	static boolean[] combinedEvidence;
	static double maxTie;
	static int[][] evidenceStat;
	
	static int count( boolean[] evidence ) {
		int count=0;
		if( evidence != null ) {
			for( int i = 0; i < evidence.length; i++ ) {
				if( evidence[i] ) {
					count++;
				}
			}
		} else {
			count=1;
		}
		return count;
	}
	
	static double sum( boolean[] evidence, double[] weight ) {
		double w=0;
		for( int i = 0; i < evidence.length; i++ ) {
			if( evidence[i] ) {
				w += weight[i];
			}
		}
		return w;
	}
	
	static void fillEvidence( boolean[] evidence, int idx ) {
		for( int i = 0; i < MAX; i++ ) {
			if( evidence[i] ) {
				evidenceStat[i][idx]++;
			}
		}
	}
	
	static class Prediction implements Comparable<Prediction>{
		static boolean first=true;
		
		static class CDS {
			String line;
			int start, end;
			CDS( String line ) {
				this.line=line;
				String[] split = line.split("\t");
				start=Integer.parseInt(split[3]);
				end=Integer.parseInt(split[4]);
			}
		}
		
		boolean altCand=false, used=false;
		String[] split;
		ArrayList<CDS> cds;
		ArrayList<String> add;
		HashMap<String,String> hash;
		int length, start, end, score;
		String contigStrand;
		HashSet<String> alternative, alternativeTranscript;
		boolean[] evidence;
		String prefix, oldId, id;
		HashSet<String> gos, defl;
		Comparable[] comp;
		
		public Prediction( String[] split, int max, int index, String prefix, HashMap<String, int[]> del ) {
			this(split, max, index, prefix, del, null, -1, -1, null);
		}
		
		public Prediction( String[] split, int max, int index, String prefix, HashMap<String,int[]> del, HashMap<String,String[]> annotInfo, int go, int defline, String[] defAttributes ) {
			evidence=new boolean[max];
			Arrays.fill(evidence, false);
			evidence[index]=true;
			
			this.split = split;
			contigStrand = split[0]+split[6];
			
			String s;
			split = split[8].split(";");
			boolean lastNotSemi = this.split[8].charAt(this.split[8].length()-1)!=';';
			hash = new HashMap<String, String>();
			String par=null;
			for( int i = 0; i < split.length; i++ ) {
				int idx = split[i].indexOf('=');
				s = split[i].substring(idx+1);
				
				String key = split[i].substring(0,idx);
				int[] stat = del.get(key);
				if( stat!= null ) {
					stat[0]++;
					if( key.equals("Parent") ) {
						par=s;
					}
					String delete = split[i] +((i+1==split.length && lastNotSemi)?"":";");
					this.split[8]=this.split[8].replace( delete, "" );
				} else {
					if( s.charAt(0)=='?' ) {
						s = "NA";
					}
					hash.put( key, s );
				}
			}			
			String sc = hash.get("score");
			if( sc == null ) {
				score = Integer.MIN_VALUE;
			} else {
				try {
					score = Integer.parseInt(sc);
				}catch(NumberFormatException nfe) {
					score = (int) Double.parseDouble(sc);
				}
			}
			oldId = hash.get("ID");
			String rg = hash.get("ref-gene");
			if( rg == null ) {
				rg = par!= null ? par : oldId;
				hash.put("ref-gene", rg);
				this.split[8] = this.split[8].replace("ID="+oldId,"ID=" + oldId + ";ref-gene=" + hash.get("ref-gene"));
			}
			if( prefix != null && prefix.length()>0 ) {
				this.prefix = prefix;
				if( prefix.charAt(prefix.length()-1) != '_' ) {
					this.prefix += "_";
				}
				id = this.prefix + oldId;
				
				this.split[8] = this.split[8].replace("ID="+oldId,"ID=" + id);
				this.split[8] = this.split[8].replace("ref-gene="+rg,"ref-gene=" + this.prefix + rg);
				hash.put("ref-gene",this.prefix + rg);
			} else {
				id=oldId;
			}
			
			//avoid Exceptions with filter conditions that access attributes that are not defined
			if( defAttributes!= null ) {
			String nan= "" + Double.NaN;
			for( int i = 0; i <defAttributes.length; i++ ) {
				if( defAttributes[i] != null ) {
					if( !hash.containsKey(defAttributes[i]) ) hash.put(defAttributes[i], nan);
				}
			}
			}
			
			if( annotInfo!= null && annotInfo.size()>0 ) {
				gos = new HashSet<String>();
				defl = new HashSet<String>();

				int i = -1;
				String[] tab = null;
				while( (i=oldId.indexOf(".",i+1))>=0 && (tab=annotInfo.get(oldId.substring(0, i)))==null );
				if( tab!= null ) {
					String goValue = tab.length>go ? tab[go] : null;
					if( goValue != null ) {
						for( String e: goValue.split(",") ) {
							if( e.length()> 0 ) {
								gos.add(e);
							}
						}
					}
					if( tab.length>defline ) {
						defl.add( tab[12].replaceAll(";", ",") );
					}
				}  else {
					System.out.println("Could not identify annotation info for: " + oldId); //TODO
				}
			}
			
			cds = new ArrayList<CDS>();
			length = 0;
			start = Integer.MAX_VALUE;//Integer.parseInt(split[3]);
			end = Integer.MIN_VALUE; //Integer.parseInt(split[4]);
			alternative = new HashSet<String>();
			alternativeTranscript = new HashSet<String>();
			add=null;
		}
		
		public int getIndex() {
			int i = 0; 
			while( i < evidence.length && !evidence[i] ) {
				i++;
			}
			return i;
		}
		
		void addAdd( String line ) {
			if( add==null ) {
				add=new ArrayList<String>();
			}
			add.add(line);
		}
				
		public void setEvidenceAndWeight(double[] weight) {
			hash.put("evidence", ""+count(evidence));
			hash.put("sumWeight", ""+sum( evidence, weight ));
		}

		void addCDS( String cds ) {
			if( prefix != null && prefix.length()>0 ) {
				//bugfix: if the id could be interpreted as regex, we had some problems
				
				//old:
				//cds = cds.replaceAll("([;\t])(ID|Parent)="+oldId, "$1$2=" + id);
				
				//middle
				//cds = cds.replace("\tParent="+oldId, "\tParent=" + id);
				//cds = cds.replace(";Parent="+oldId, ";Parent=" + id);
				//cds = cds.replace("\tID="+oldId, "\tID=" + id);
				//cds = cds.replace(";ID="+oldId, ";ID=" + id);
				
				//new:
				cds=replace( cds, "Parent="+oldId, "Parent="+id, 0 );
				cds=replace( cds, "ID="+oldId, "ID="+id, 0 );
			}
			CDS current = new CDS( cds );
			this.cds.add( current );
			length+= (current.end-current.start+1);
			start = Math.min(start, current.start);
			end = Math.max(end, current.end);
		}

		@Override
		public int compareTo(Prediction o) {
			int res = split[0].compareTo(o.split[0]);
			if( res == 0 ) {
				res = Integer.compare( start, o.start );
				if( res == 0 ) {
					res = Integer.compare( end, o.end );
					if( res == 0 ) {
						res = Integer.compare( cds.size(), o.cds.size() );
						if( res == 0 ) {
							CDS part1 = null, part2 = null;
							for( int i = 0; i < cds.size(); i++ ) {
								part1 = cds.get(i);
								part2 = o.cds.get(i);
								if( part1.start == part2.start ) {
									if( part1.end != part2.end ) {
										res = Integer.compare(part1.end,part2.end);
										break;
									}
								} else {
									res = Integer.compare(part1.start,part2.start);
									break;
								}
							}
						}
					}
				}
			}
			return res;
		}
		
		public void combine( Prediction n, boolean addAltTransIDs ) {
			String rg = n.hash.get("ref-gene");
			if( !rg.equals( hash.get("ref-gene") ) ) {
				alternative.add(rg);
			}
			for( int i = 0; i < evidence.length; i++ ) {
				evidence[i] |= n.evidence[i];
			}
			
			alternative.addAll( n.alternative );
			if( addAltTransIDs ) {
				alternativeTranscript.add( n.id );
				alternativeTranscript.addAll( n.alternativeTranscript );
			}
			
			if( gos != null ) {
				if( n.gos!=null ) {
					gos.addAll(n.gos);
					defl.addAll(n.defl);
				}
			} else {
				if( n.gos!=null ) {
					gos=n.gos;
					defl=n.defl;
				}
			}
		}
		
		public void write( StringBuffer w, Protocol p, int gene ) throws IOException, NumberFormatException, ScriptException {			
			int count = Integer.parseInt(hash.get("evidence"));
			String we = hash.get("sumWeight");
			String t = hash.get("tie");
			String tpc = hash.get("tpc");

			if( p != null ) {
				for( int i = 0; i < evidence.length; i++ ) {
					combinedEvidence[i] |= evidence[i];
				}
				stats[count-1][2]++;
				if( t == null || t.equals("NA") ) {
					if( tpc!= null && !tpc.equals("NA") && Double.parseDouble(tpc)==1 ) {
						stats[count-1][4]++;
					}
				} else {
					if( t != null ) {
						double tie = Double.parseDouble(t);
						if( tie == 1d ) {
							stats[count-1][3]++;
						}
						maxTie = Math.max(tie, maxTie);
					}
				}
			}
			for( int i = 0; i < split.length; i++ ) {
				if( add==null && (i==3 || i==4) ) {
					w.append( "\t"+(i==3?start:end) );
				} else {
					w.append( (i==0?"":"\t") + split[i] );
				}
			}
			String h = split[split.length-1];
			if( h.charAt(h.length()-1)!=';') {
				w.append(";");
			}
			w.append( "evidence=" + count + ";" );
			if( p!=null ) w.append( "Parent=gene_"+gene + ";");
			w.append( "sumWeight=" + we + ";" );
			
			if( p!= null ) {
				maxEvidence = Math.max(count, maxEvidence);
				if( hash.containsKey("start") && hash.containsKey("stop") ) {
					complete += (hash.get("start").charAt(0)=='M' && hash.get("stop").charAt(0)=='*') ? 1 : 0;
				} else {
					if( first ) {
						p.appendWarning("Could not find attributes 'start' or 'stop' for at least one prediction. You can use AnnotationEvidence to add these to the input file." );
						first = false;
					}
				}
			}
			
			if( alternative.size() > 0 ) {
				String[] it = alternative.toArray(new String[0]);
				Arrays.sort(it);
				w.append( "alternative=\"" + it[0] );
				for( int i = 1; i < it.length; i++ ) {
					w.append("," + it[i] );
				}
				w.append( "\";" );
			}
			if( alternativeTranscript.size() > 0 ) {
				String[] it = alternativeTranscript.toArray(new String[0]);
				Arrays.sort(it);
				w.append( "alternativeTranscript=\"" + it[0] );
				for( int i = 1; i < it.length; i++ ) {
					w.append("," + it[i] );
				}
				w.append( "\";" );
			}
			if( gos!= null ) {
				String[] goTerms = gos==null ? new String[0] : gos.toArray(new String[0]);
				if( goTerms != null && goTerms.length > 0 )
				{
					h = Arrays.toString(goTerms).replace("[","\"").replace("]", "\"");
					w.append("potential_GO=" + h+ ";");
				}
				String[] defTerms = defl==null ? new String[0] : defl.toArray(new String[0]);
				if( defTerms != null && defTerms.length > 0 )
				{
					h = Arrays.toString(defTerms).replace("[","\"").replace("]", "\"");
					w.append("potential_defline=" + h);
				}
			}
			w.append("\n");//w.newLine();
			if( add!= null ) {
				for( int i = 0; i < add.size(); i++ ) {
					String line = add.get(i);
					//adjust Parent
					if( oldId!=id ) {
						line = line.replace("Parent="+oldId,"Parent=" + id);
					}
					w.append( line );
					w.append("\n");//w.newLine();
				}
			}
			for( int i = 0; i < cds.size(); i++ ) {
				w.append( cds.get(i).line );
				w.append("\n");//w.newLine();
			}
			
			if( p!=null ) fillEvidence(evidence, 1);
		}
		
		public String toString() {
			String r = "";
			for( int i = 0; i < split.length; i++ ) {
				r += (i==0?"":"\t") + split[i];
			}
			return r;
		}
		
		HashSet<Integer> s, e;
		
		public double commonBorders(Prediction o) {
			int anz = 0;
			CDS c;
			if( s== null ) {
				s = new HashSet<Integer>();
				e = new HashSet<Integer>();
				for( int i = 0; i < cds.size(); i++ ) {
					c = cds.get(i);
					s.add(c.start);
					e.add(c.end);
				}
			}
			for( int i = 0; i < o.cds.size(); i++ ) {
				c = o.cds.get(i);
				
				if( s.contains(c.start) ) {
					anz++;
				}
				if( e.contains(c.end) ) {
					anz++;
				}
			}
			return anz;
		}

		void fillNuc( BitSet nuc, int s, int e ) {
			for( int i = 0; i < cds.size(); i++ ) {
				CDS c = cds.get(i);
				try {
					nuc.set( c.start-s, c.end-s+1 );
				} catch( IndexOutOfBoundsException ioobe ) {
					System.out.println(Arrays.toString(split));
					System.out.println(c.line);
					throw ioobe;
				}
			}
		}

		int overlap( BitSet b, int clusterStart, int clusterEnd ) {
			int dir = split[6].charAt(0)=='+'?+1:-1;

			int overlap = 0;
			for( int i = dir==1?0:cds.size()-1; dir==1 ? (i < cds.size()) : (i>=0); i+=dir ) {
				CDS c = cds.get(i);
				int s = c.start-clusterStart;
				int e = c.end-clusterStart;
				
				int set = b.nextSetBit( s );
				//overlap
				while( set >=0 && set < e ) {
					int unset = b.nextClearBit(set);
					if( unset > e ) {
						overlap += (e-set+1);
						break;
					} else {
						overlap += (unset-set);
						set = b.nextSetBit( unset );
					}
				}
				if( set < 0 ) break;
			}
			return overlap;
		}
		
		/*int phaseCompatible( Prediction p ) {
			//TODO read this from gff?
			int phase = 0;
			int pPhase = 0;

			int overlap=0, common=0;
			int i = 0, j = 0;
			while( i<cds.size() && j<p.cds.size() ) {
				CDS c = null;
				CDS d = p.cds.get(j);
				while( i<cds.size() && (c=cds.get(i)).end < d.start ) {
					phase = (phase+ (c.end-c.start+1) ) % 3;
					i++;
					common=0;
				}
				if( i==cds.size() ) break;
				
				while( j<p.cds.size() && (d=p.cds.get(j)).end < c.start ) {
					pPhase = (pPhase+ (d.end-d.start+1) ) % 3;
					j++;
					common=0;
				}
				if( j==p.cds.size() ) break;

				//TODO test
				//if( d.end < c.start ) //does not happen
				if( c.end < d.start ) {
					//skip
					common=0;
				} else {
					//overlap
					//c.start,d.start < c.end,d.end
					int start = Math.max(d.start,c.start);
					int end = Math.min(d.end,c.end);
					int len = end-start;
					
					//if( c.start==d.start && common>0 ) common=;
					
					int compPhase = (phase+start-c.start) % 3;
					int compPPhase = (pPhase+start-d.start) % 3;
					
					if( compPhase == compPPhase ) {
						overlap += ((len-compPhase)/3)*3;
					}
					
					//if( c.end==d.end ) common=;
					
					if( c.end==end ) {
						phase = (phase+ (c.end-c.start+1) ) % 3;
						i++;
					}
					if( d.end==end ) {
						pPhase = (phase+ (d.end-d.start+1) ) % 3;
						j++;
					}
				}
			}
			return overlap;
		}*/
	}


	@Override
	public ToolParameterSet getToolParameters() {
		try{
			ExpandableParameterSet attribute = new ExpandableParameterSet( new SimpleParameterSet(		
					new SimpleParameter(DataType.STRING,"attribute","An attribute used for the kmeans clustering", true, new RegExpValidator("[a-zA-Z 0-9\\.,()><=!'\\-\\+\\*\\/]*") )
				), "attribute", "", 1 );
			((SimpleParameterSet)attribute.getParameterAt(0).getValue()).getParameterAt(0).setValue("Math.abs(maxScore-Math.max(0,score))/maxScore");
			//attribute.addParameterToSet();
			//((SimpleParameterSet)attribute.getParameterAt(1).getValue()).getParameterAt(0).setValue("Math.abs(raa-aa)/raa");			
			
			return
				new ToolParameterSet( getShortName(),
					new SimpleParameter(DataType.STRING,"tag","the tag used to read the GeMoMa annotations",true,GeMoMa.TAG),
					new ParameterSetContainer( "predicted annotation", "", new ExpandableParameterSet( new SimpleParameterSet(		
							new SimpleParameter(DataType.STRING,"prefix","the prefix can be used to distinguish predictions from different input files", false, new RegExpValidator("\\w*") ),
							new SimpleParameter(DataType.DOUBLE,"weight","the weight can be used to prioritize predictions from different input files; each prediction will get an additional attribute sumWeight that can be used in the filter", false, new NumberValidator<Double>(0d, 1000d), 1d),
							new FileParameter( "gene annotation file", "GFF file containing the gene annotations (predicted by GeMoMa)", "gff,gff3",  true, new FileExistsValidator(), true ),
							new FileParameter( "annotation info", "annotation information of the reference, tab-delimited file containing at least the columns transcriptName, GO and .*defline", "tabular", false, new FileExistsValidator() )
						), "gene annotations", "", 1 ) ),
					new SimpleParameter(DataType.STRING,"default attributes","Comma-separated list of attributes that is set to NaN if they are not given in the annotation file. "
							+ "This allows to use these attributes for sorting or filter criteria. "
							+ "It is especially meaningful if the gene annotation files were received fom different software packages (e.g., GeMoMa, Braker, ...) having different attributes.", true, "tie,tde,tae,iAA,pAA,score,lpm,maxGap,bestScore,maxScore,raa,rce" ),				
					new SelectionParameter(DataType.PARAMETERSET, new String[] {"NO","YES"},
							new ParameterSet[] {
									new SimpleParameterSet(),
									new SimpleParameterSet( 
											new SimpleParameter( DataType.INT, "minimal number of predictions", "only gene sets with at least this number of predictions will be used for clustering", true, new NumberValidator<Integer>(0, 100000000), 1000),
											new SimpleParameter( DataType.INT, "cluster", "the number of clusters to be used for kmeans", true, new NumberValidator<Integer>(2, 100), 2),
											new SimpleParameter( DataType.INT, "good cluster", "the number of good clusters, good clusters are those with small mean, all members of a good cluster are further used", true, new NumberValidator<Integer>(1, 99), 1),
											//new ParameterSetContainer( "cluster attribute", "", attribute )
											
											new SelectionParameter(DataType.PARAMETERSET, new String[] {"GLOBAL","LOCAL"},
													new ParameterSet[] {
															new SimpleParameterSet(),
															new SimpleParameterSet( 
																	new SimpleParameter( DataType.INT, "margin", "the number of bp upstream and downstream of a predictions used to identify neighboring predictions for the statistics", true, new NumberValidator<Integer>(0, 100000000), 1000000),
																	new SimpleParameter( DataType.DOUBLE, "quantile", "the quantile used for the local trend", true, new NumberValidator<Double>(0d,1d), 0.2)
															)
													},
													"trend",
													"whether a local component should be used for the cluster attribute (might be helpful for regions with different conservation (e.g. introgressions in chromosomes))",
													true )
									)
							},
							"kmeans",
							"whether kmeans should be performed for each input file and clusters with large mean distance to the origin will be discarded",
							true ),
					new SimpleParameter(DataType.STRING,"filter","A filter can be applied to transcript predictions to receive only reasonable predictions. The filter is applied to the GFF attributes. "
							+ "The default filter decides based on the completeness of the prediction (start=='M' and stop=='*') and the relative score (score/aa>=0.75) whether a prediction is used or not. In addition, predictions without score (isNaN(score)) will be used as external annotations do not have the attribute score. "
							+ "Different criteria can be combined using 'and' and 'or'. "
							+ "A more sophisticated filter could be applied for instance using the completeness, the relative score, the evidence as well as statistics based on RNA-seq: start=='M' and stop=='*' and score/aa>=0.75 and (evidence>1 or tpc==1.0)", false, new RegExpValidator("[a-zA-Z 0-9\\.()><=!'\\-\\+\\*\\/]*"), "start=='M' and stop=='*' and (isNaN(score) or score/aa>=0.75)" ),
					new SimpleParameter(DataType.STRING,"sorting","comma-separated list that defines the way predictions are sorted within a cluster", true, "sumWeight,score,aa"),
					new SimpleParameter(DataType.BOOLEAN, "intermediate result", "a switch to decide whether an intermediate result of filtered predictions that are not combined to genes should be returned", true, false),
					new SimpleParameter(DataType.DOUBLE,"length difference","maximal percentage of length difference between the representative transcript and an alternative transcript, alternative transcripts with a higher percentage are discarded", false, new NumberValidator<Double>(0d, 10000d) ),
					new SimpleParameter(DataType.STRING,"alternative transcript filter","If a prediction is suggested as an alternative transcript, this additional filter will be applied. The filter works syntactically like the 'filter' parameter and allows you to keep the number of alternative transcripts small based on meaningful criteria. "
							+ "Reasonable filter could be based on RNA-seq data (tie==1) or on sumWeight (sumWeight>1). "
							+ "A more sophisticated filter could be applied combining several criteria: tie==1 or sumWeight>1", false, new RegExpValidator("[a-zA-Z 0-9\\.()><=!?:'\\-\\+\\*\\/_]*"), "tie==1 or sumWeight>1" ),
					new SimpleParameter(DataType.DOUBLE,"common border filter","the filter on the common borders of transcripts, the lower the more transcripts will be checked as alternative splice isoforms", true, new NumberValidator<Double>(0d, 1d), 0.75 ),
					new SimpleParameter(DataType.INT,"maximal number of transcripts per gene","the maximal number of allowed transcript predictions per gene", true, new NumberValidator<Comparable<Integer>>(1, Integer.MAX_VALUE), Integer.MAX_VALUE ),
					new SimpleParameter(DataType.BOOLEAN, "add alternative transcripts", "a switch to decide whether the IDs of alternative transcripts that have the same CDS should be listed for each prediction", true, false),
					new SimpleParameter( DataType.BOOLEAN, "transfer features", "if true, additional features like UTRs will be transfered from the input to the output. Only features of the representatives will be transfered. The UTRs of identical CDS predictions listed in \"alternative\" will not be transfered or combined", true, false )
				);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public String getToolName() {
		return "GeMoMa Annotation Filter";
	}

	@Override
	public String getShortName() {
		return "GAF";
	}

	@Override
	public String getDescription() {
		return "filters most reliable gene predictions";
	}

	@Override
	public String getHelpText() {
		return
				"This tool combines and filters gene predictions from different sources yielding a common gene prediction."
				+ " The tool does not modify the predictions, but filters redundant and low-quality predictions and selects relevant predictions."
				+ " In addition, it adds attributes to the annotation of transcript predictions.\n\n"
				
				+ "The algorithm does the following:\n"
				+ "First, redundant predictions are identified (and additional attributes (evidence, sumWeight) are introduced).\n"
				+ "Second, the predictions are filtered using the user-specified criterium based on the attributes from the annotation.\n"
				+ "Third, clusters of overlapping predictions are determined, the predictions are sorted within the cluster and relevant predictions are extracted.\n\n"
				
				+ "Optionally, annotation info can be added for each reference organism enabling a functional prediction for predicted transcripts based on the function of the reference transcript.\n"
				+ "Phytozome provides annotation info tables for plants, but annotation info can be used from any source as long as they are tab-delimited files with at least the following columns: transcriptName, GO and .*defline.\n\n"
				
				+ "Initially, GAF was build to combine gene predictions obtained from **GeMoMa**. It allows to combine the predictions from multiple reference organisms, but works also using only one reference organism. "
				+ "However, GAF also allows to integrate predictions from ab-initio or purely RNA-seq-based gene predictors as well as manually curated annotation. "
				+ "If you like to do so, we recommend to run **AnnotationEvidence** for each of these input files to add additional attributes that can be used for sorting and filtering within **GAF**. "
				+ "The sort and filter criteria need to be carefully revised in this case. Default values can be set for missing attributes."
				
				+ MORE;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", "filtered predictions"),
		};
	}
	
	@Override
	public ToolResult[] getTestCases( String path ) {
		try {
			return new ToolResult[]{new ToolResult(FileManager.readFile(path+File.separator+"tests/gemoma/xml/gaf-test.xml"))};
		} catch( Exception e ) {
			e.printStackTrace();
			return null;
		}
	}
}