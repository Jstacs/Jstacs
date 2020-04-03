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
import java.util.Map.Entry;

import javax.script.Bindings;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import de.jstacs.DataType;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.parameters.validation.RegExpValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.IntList;

/**
 * This class allows to filter GeMoMa annotation files (GFFs) to obtain the most relevant predictions.
 * 
 * @author Jens Keilwagen
 * 
 * @see GeMoMa
 */
public class GeMoMaAnnotationFilter extends GeMoMaModule {
	
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
			asDouble.add("score");
			
			asDouble.add("tae");
			asDouble.add("tde");
			asDouble.add("tie");
			asDouble.add("minSplitReads");
			asDouble.add("tpc");
			asDouble.add("minCov");
			asDouble.add("avgCov");
			
			asDouble.add("iAA");
			asDouble.add("pAA");
			
			asDouble.add("evidence");
		}
		
		String[] keys;
		HashSet<String> errors;
		Protocol protocol;
		
		public UserSpecifiedComparator( String specified, Protocol protocol ) {
			keys=specified.split(",");
			for( int i = 0; i < keys.length; i++ ) {
				keys[i] = keys[i].trim();
			}
			errors = new HashSet<String>();
			this.protocol = protocol;
		}
		
		public static double parse( String d ) {
			if( d.equals( "NA" ) ) {
				return Double.NaN;
			} else {
				return Double.parseDouble(d);
			}
		}
		
		public static int compareAsDoubles( String d1, String d2 ) {
			return -Double.compare(parse(d1), parse(d2));
		}
		
		public int compare(Prediction o1, Prediction o2) {
			int d=0, i = 0;
			while( i < keys.length && d==0 ) {
				String v1 = o1.hash.get(keys[i]);
				String v2 = o2.hash.get(keys[i]);
				if( v1 != null && v2 !=null ) {
					if( asDouble.contains(keys[i]) ) {
						d = compareAsDoubles(v1, v2);
					} else {
						d = v1.compareTo( v2 );
					}
				} else {
					if( !errors.contains(keys[i])) {
						errors.add(keys[i]);
						protocol.appendWarning("Could not find key (" + keys[i] + ") that is used for sorting\n");
					}
				}
				i++;
			}
			
			if( d == 0 ) {
				d = o1.id.compareTo(o2.id);
			}
			return d;
		}
	}
	
	private static String prepareFilter( String filter ) {
		if( filter == null ) {
			filter = "";
		} else {
			filter = filter.trim();
		}
		filter = filter.replaceAll( " or ", " || " );
		filter = filter.replaceAll( " and ", " && " );
		return filter;
	}
	
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		
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
		String filter = prepareFilter( (String) parameters.getParameterForName("filter").getValue() );
		UserSpecifiedComparator comp = new UserSpecifiedComparator( (String) parameters.getParameterForName("sorting").getValue(), protocol );
		String altFilter = prepareFilter( (String) parameters.getParameterForName("alternative transcript filter").getValue() );
				
		ScriptEngineManager mgr = new ScriptEngineManager();
		ScriptEngine engine = mgr.getEngineByName("nashorn");
	    
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterForName("predicted annotation").getValue();
		
		String line, t;
		String[] split;
				
		ArrayList<Prediction> pred = new ArrayList<Prediction>();
		ArrayList<Prediction> currentCluster = new ArrayList<Prediction>();
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
			weight[k] = (Double) sps.getParameterAt(1).getValue();
			
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
					throw new IllegalArgumentException("annotation info ("+ k+ ") must be a tab-delimited file with at least the following columns: transcriptName, GO, and .*defline");
				}
				while( (line = r.readLine()) != null ) {
					split = line.split("\t");
					annotInfo.put(split[transcript].toUpperCase(), split);
				}
				r.close();
			}
			
			String fName = sps.getParameterAt(2).getValue().toString();
			//System.out.println(fName);
			hash.clear();
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
				if( t.equalsIgnoreCase( tag ) ) { //tag: prediction/mRNA/transcript
					current = new Prediction(split, k, prefix[k], annotInfo, go, defline, defAttributes);
					if( ids.contains( current.id ) ) {
						//System.out.println(line);
						//System.out.println(current);
						r.close();
						throw new IllegalArgumentException("The id (" + current.id + ") has been used before. You can try to use a prefix (cf. parameters)." );
					} else {
						ids.add(current.id);
					}
					pred.add( current );
					rnaSeq |= current.hash.containsKey("tie");
				} else if( t.equals( "CDS" ) ) {
					current.addCDS( line );
				}
				
			}
			r.close();
			if( hash.size() != 1 ) {
				if( hash.size() == 0  ) {
					protocol.appendWarning("Evidence file " + k + " contains no parameter description\n" );
				} else {
					protocol.appendWarning("Evidence file " + k + " contains a different parameter description\n" );
				}
			}
			int delta = 0;
			sb.delete(0,sb.length());
			for( int i =0; i < att.length; i++ ) {
				int[] v = del.get(att[i]);
				delta += v[0];
				sb.append( (i==0?"":", ") + v[0] + " \"" + att[i] + "\"" );
				//delete values => prepare next round
				v[0]=0;
			}
			if( delta>0 ) {				
				protocol.appendWarning("Delete attributes in evidence file " + k + ": " + sb  + "\n");
			}
			Iterator<String> it = hash.iterator();
			while( it.hasNext() ) {
				allInfos.add(it.next());
			}
			
			//System.out.println( pred.size() + "\t" + k + "\t" +current.index);
		}
		protocol.append( "all: " + pred.size() + "\n" );

		//combine redundant
		Collections.sort(pred);
		Prediction last = pred.get(pred.size()-1);
		for( int i = pred.size()-2; i >= 0; i-- ){
			Prediction cur = pred.get(i);
			if( cur.compareTo(last) == 0 ) {
				//identical
				if( cur.score > last.score || (cur.score== last.score && cur.id.compareTo(last.id)<0) ) {
					cur.combine(last);
					pred.remove(i+1);
					last=cur;
				} else {
					last.combine(cur);
					pred.remove(i);
				}
			} else {
				last=cur;
			}
		}
		
		//filter: user-specified
		for( int i = pred.size()-1; i >= 0; i-- ){
			Prediction p = pred.get(i);
			if( !filter(filter, engine, p, weight) ) {
				//System.out.println(p.hash.toString());
				pred.remove(i);
			} else {
				fillEvidence(p.evidence, p.index, 0);
			}
		}
		protocol.append( "filtered: " + pred.size() + "\n" );
		
		File out = Tools.createTempFile("GAF-filtered");
		int clustered=0, transcripts = 0;
		if( pred.size() > 0 ) {
			//cluster predictions			
			BufferedWriter w = new BufferedWriter( new FileWriter(out) );
			w.append("##gff-version 3");
			w.newLine();
			for( int i = 0; i < allInfos.size(); i++ ) {
				w.append( allInfos.get(i) );
				w.newLine();
			}
			w.append(INFO + getShortName() + " " + getToolVersion() + "; ");
			String info = JstacsTool.getSimpleParameterInfo(parameters);
			if( info != null ) {
				w.append("SIMPLE PARAMETERS: " + info );
			}
			w.newLine();
			
			int i = 0;
			Prediction next = pred.get(i);
			while( i < pred.size() ) {
				current = next;
				int end = current.end, alt = i;
				i++;	
	
				currentCluster.clear();
				currentCluster.add(current);
				while( i < pred.size() && (next = pred.get(i)).contigStrand.equals(current.contigStrand)
						//&& (end-Integer.parseInt(next.split[3])+1d)/(end-start+1d) > 0.1
						&& end>Integer.parseInt(next.split[3])
				) {
					end = Math.max(end,Integer.parseInt(next.split[4]));
					currentCluster.add(next);
					i++;
				}
				
				Collections.sort( currentCluster, comp );
				int p=create( rnaSeq, maxTranscripts, currentCluster, w, altFilter, engine, weight, cbTh, protocol );
				clustered++;
				transcripts+=p;
				
				transcripts += checkNestedSameStrand( rnaSeq, maxTranscripts, pred, alt, i, w, altFilter, engine, weight, cbTh, comp, protocol );
			}
			w.close();
		}
				
		protocol.append( "clustered: " + clustered + "\n\n" );
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
		
		
		return new ToolResult("", "", null, new ResultSet(new TextResult("filtered predictions", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "gff", getToolName(), null, true)), parameters, getToolName(), new Date());
	}
	
	static int checkNestedSameStrand( boolean rnaSeq, int maxTranscripts, ArrayList<Prediction> list, int start, int end, BufferedWriter w, String altFilter, ScriptEngine engine, double[] weight, double cbTh, Comparator<Prediction> comp, Protocol protocol ) throws Exception {
		IntList used = new IntList(), notUsed = new IntList();
		int s = Integer.MAX_VALUE, e = 0;
		for( int i = start; i < end; i++ ) {
			Prediction p = list.get(i);
			if( p.start < s ) {
				s = p.start;
			}
			if( p.end > e ) {
				e = p.end;
			}
			if( p.used ) {
				used.add(i);
			} else {
				notUsed.add(i);
			}
		}
		
		if( used.length() > 0 && notUsed.length()>0 ) {
			//fill BitSet
			for( int i = start; i < end; i++ ) {
				list.get(i).fillNuc(s,e);
			}
			
			ArrayList<Prediction> cand = new ArrayList<Prediction>();
			for( int n = 0; n < notUsed.length(); n++ ) {
				Prediction not = list.get(notUsed.get(n));
				boolean overlap = false;
				for( int u = 0; u < used.length(); u++ ) {
					Prediction use = list.get(used.get(u));
					overlap |= use.overlap(not);
				}
				if( !overlap ) {
					not.discard=false;
					cand.add(not);
				}
			}
			
			//delete BitSet
			for( int i = start; i < end; i++ ) {
				list.get(i).clear();
			}		
			
			if( cand.size() > 0 ) {
				Collections.sort( cand, comp );
				return create(rnaSeq, maxTranscripts, cand, w, altFilter, engine, weight, cbTh, protocol );
			}
		}
		return 0;
	}
	
	
	static int h = 0;
	
	static int[][] stats;
	static int gene=0;
	
	static int create( boolean rnaSeq, int maxTranscripts, ArrayList<Prediction> list, BufferedWriter w, String filter, ScriptEngine engine, double[] weight, double cbTh, Protocol protocol ) throws Exception {
		if( list.size() == 0 ) {
			return 0;
		}
		Prediction current = list.remove(0);
		ArrayList<Prediction> used = new ArrayList<Prediction>();
		used.add( current );
		current.discard = true;
		current.used=true;

		int count = 0, i = 0;
		while( i < list.size() ) {
			Prediction n = list.get(i);
			if( n.end < current.start || current.end < n.start ) {//no overlap
				count++;
				i++;
			} else { //overlapping
				double min=2d*Math.min(n.cds.size(),current.cds.size());
				double cb = current.commonBorders(n);
				
				int j = 0;
				if( cb/min>=cbTh ) {//hinreichende Überlappung
					while( j < used.size() ) {
						Prediction x = used.get(j); 
						cb = x.commonBorders(n);
//						double max=2d*Math.max(n.cds.size(),x.cds.size());								
						j++;
					}
				}
				
				if( j == used.size() ) {
					if( filter(filter, engine, n, weight) ) {
						used.add(n);
						n.used=true;
					}
				}
				list.remove(i);
			}
		}
		
		maxEvidence=0;
		Arrays.fill( combinedEvidence, false );
		maxTie=-1;
		complete=0;
		st=Integer.MAX_VALUE;
		en=Integer.MIN_VALUE;
		Prediction n=null;
		
		for( i = 0; i < used.size() && i < maxTranscripts; i++ ) {
			n = used.get(i);
			n.write(w, protocol);
		}
		if( i>0 ) {
			//write
			fillEvidence(combinedEvidence, -999, 2);
			
			w.append(n.split[0] + "\tGAF\tgene\t" + st + "\t" + en  + "\t.\t" + n.split[6] + "\t.\tID=gene_"+gene+";transcripts=" + i + ";complete="+complete+ (rnaSeq?";maxTie=" + (maxTie<0?"NA":maxTie):"") +(MAX>1?";maxEvidence="+maxEvidence+";combinedEvidence=" + count(combinedEvidence):"") );
			w.newLine();
			
			stats[maxEvidence-1][0]++;
			if( maxTie == 1 ) {
				stats[maxEvidence-1][1]++;
			}
			gene++;
		}
		
		if( count>0  ) {
			i += create(rnaSeq, maxTranscripts, list, w, filter, engine, weight, cbTh, protocol);
		}

		return i;
	}
	
	static int maxEvidence, st, en, complete;
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
	
	static double sum( boolean[] evidence, int index, double[] weight ) {
		double w=0;
		if( evidence != null ) {
			for( int i = 0; i < evidence.length; i++ ) {
				if( evidence[i] ) {
					w += weight[i];
				}
			}
		} else {
			w=weight[index];
		}
		return w;
	}
	
	static boolean filter( String filter, ScriptEngine engine, Prediction toBeChecked, double[] weight ) throws ScriptException {
		int count=count(toBeChecked.evidence);
		toBeChecked.hash.put("evidence", ""+count);
		double w = sum( toBeChecked.evidence, toBeChecked.index, weight );
		toBeChecked.hash.put("sumWeight", ""+w);
		
		if( filter.length()== 0 ) {
			return true;
		} else {
			Bindings b = engine.createBindings();
			Iterator<Entry<String,String>> it = toBeChecked.hash.entrySet().iterator();
			while( it.hasNext() ) {
				Entry<String,String> e = it.next();
				String key = e.getKey();
				if( filter.indexOf(key)>=0 ) {
					String val = e.getValue();
					b.put(key, val.equals("NA")?""+Double.NaN:val);
				}
			}			
			
			String s = engine.eval(filter, b).toString();
			Boolean bool = Boolean.parseBoolean( s );
//System.out.println(toBeChecked.id + "\t" + f + "\t" + s + "\t" + bool);
			return bool;
		}
	}
	
	static void fillEvidence( boolean[] evidence, int index, int idx ) {
		if( evidence == null ) {
			evidenceStat[index][idx]++;
		} else {
			for( int i = 0; i < MAX; i++ ) {
				if( evidence[i] ) {
					evidenceStat[i][idx]++;
				}
			}
		}
	}
	
	static class Prediction implements Comparable<Prediction>{
		static boolean first=true;
		
		boolean discard = false, used=false;
		String[] split;
		ArrayList<String[]> cds;
		HashMap<String,String> hash;
		int length, start, end, score;
		String contigStrand;
		HashSet<String> alternative;
		int index;
		boolean[] evidence;
		String prefix, oldId, id;
		HashSet<String> gos, defl;
		
		public Prediction( String[] split, int index, String prefix, HashMap<String,String[]> annotInfo, int go, int defline, String[] defAttributes ) {
			this.index = index;
			
			this.split = split;
			contigStrand = split[0]+split[6];
			start = Integer.parseInt(split[3]);
			end = Integer.parseInt(split[4]);
			
			String s;
			split = split[8].split(";");
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
					String delete = split[i] +(i+1==split.length?"":";");
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
				score = Integer.parseInt(sc);
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
			String nan= "" + Double.NaN;
			for( int i = 0; i <defAttributes.length; i++ ) {
				if( defAttributes[i] != null ) {
					if( !hash.containsKey(defAttributes[i]) ) hash.put(defAttributes[i], nan);
				}
			}
			
			if( annotInfo.size()>0 ) {
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
						defl.add( tab[12] );
					}
				}  else {
					System.out.println(oldId); //TODO
				}
			}
			
			cds = new ArrayList<String[]>();
			length = 0;
			alternative = new HashSet<String>();
		}
				
		void addCDS( String cds ) {
			if( prefix != null && prefix.length()>0 ) {
				//bugfix: if the id could be interpreted as regex, we had some problems
				
				//old:
				//cds = cds.replaceAll("([;\t])(ID|Parent)="+oldId, "$1$2=" + id);
				
				//new:
				cds = cds.replace("\tParent="+oldId, "\tParent=" + id);
				cds = cds.replace(";Parent="+oldId, ";Parent=" + id);
				cds = cds.replace("\tID="+oldId, "\tID=" + id);
				cds = cds.replace(";ID="+oldId, ";ID=" + id);
			}
			String[] split = cds.split("\t");
			this.cds.add( split );
			length+= (Integer.parseInt(split[4])-Integer.parseInt(split[3])+1);
		}

		@Override
		public int compareTo(Prediction o) {
			int res = contigStrand.compareTo(o.contigStrand);
			if( res == 0 ) {
				res = Integer.compare( start, o.start );
				if( res == 0 ) {
					res = Integer.compare( end, o.end );
					if( res == 0 ) {
						res = Integer.compare( cds.size(), o.cds.size() );
						if( res == 0 ) {
							int i = 0;
							String[] part1 = null, part2 = null;
							while( i < cds.size() ) {
								part1 = cds.get(i);
								part2 = o.cds.get(i);
								if( part1[3].equals(part2[3]) && part1[4].equals(part2[4]) ) {
									i++;
								} else {
									break;
								}
							}
							if( i < cds.size() ) {
								res = (part1[3]+"\t"+part1[4]).compareTo( part2[3]+"\t"+part2[4] );
							}
						}
					}
				}
			}
			return res;
		}
		
		public void combine( Prediction n ) {
			String rg = n.hash.get("ref-gene");
			if( !rg.equals( hash.get("ref-gene") ) ) {
				alternative.add(rg);
			}
			if( evidence == null ) {
				evidence = new boolean[MAX];
				Arrays.fill(evidence, false);
				if( n.evidence == null ) {
					evidence[n.index]=true;
				} else {
					System.arraycopy(n.evidence, 0, evidence, 0, evidence.length );
				}
				evidence[index]=true;
			} else {
				if( n.evidence == null ) {
					evidence[n.index]=true;
				} else {
					for( int i = 0; i < evidence.length; i++ ) {
						evidence[i] |= n.evidence[i];
					}
				}
			}
			
			alternative.addAll( n.alternative );
			
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
		
		public void write( BufferedWriter w, Protocol p ) throws IOException, NumberFormatException, ScriptException {
			if( evidence != null ) {
				for( int i = 0; i < evidence.length; i++ ) {
					combinedEvidence[i] |= evidence[i];
				}
			} else {
				combinedEvidence[index]=true;
			}
		
			int count = Integer.parseInt(hash.get("evidence"));
			String we = hash.get("sumWeight");
			String t = hash.get("tie");
			String tpc = hash.get("tpc");

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
			for( int i = 0; i < split.length; i++ ) {
				w.append( (i==0?"":"\t") + split[i] );
			}
			String h = split[split.length-1];
			if( h.charAt(h.length()-1)!=';') {
				w.append(";");
			}
			w.write( "evidence=" + count + ";Parent=gene_"+gene + ";sumWeight=" + we + ";" );
			
			maxEvidence = Math.max(count, maxEvidence);
			st=Math.min(st, start);
			en=Math.max(en, end);
			if( hash.containsKey("start") && hash.containsKey("stop") ) {
				complete += (hash.get("start").charAt(0)=='M' && hash.get("stop").charAt(0)=='*') ? 1 : 0;
			} else {
				if( first ) {
					p.appendWarning("Could not find attributes 'start' or 'stop' for at least one prediction. You can use AnnotationEvidnece to add these to the input file." );
					first = false;
				}
			}
			
			if( alternative.size() > 0 ) {
				String[] it = alternative.toArray(new String[0]);
				Arrays.sort(it);
				w.write( "alternative=\"" + it[0] );
				for( int i = 1; i < it.length; i++ ) {
					w.write("," + it[i] );
				}
				w.write( "\";" );
			}
			if( gos!= null ) {
				String[] goTerms = gos==null ? new String[0] : gos.toArray(new String[0]);
				if( goTerms != null && goTerms.length > 0 )
				{
					h = Arrays.toString(goTerms).replace("[","\"").replace("]", "\"");
					w.write("potential_GO=" + h+ ";");
				}
				String[] defTerms = defl==null ? new String[0] : defl.toArray(new String[0]);
				if( defTerms != null && defTerms.length > 0 )
				{
					h = Arrays.toString(defTerms).replace("[","\"").replace("]", "\"");
					w.write("potential_defline=" + h);
				}
			}
			w.newLine();
			for( int i = 0; i < cds.size(); i++ ) {
				String[] split = cds.get(i);
				for( int j = 0; j < split.length; j++ ) {
					w.append( (j==0?"":"\t") + split[j] );
				}
				w.newLine();
			}
			
			fillEvidence(evidence, index, 1);
		}
		
		public String toString() {
			String r = "";
			for( int i = 0; i < split.length; i++ ) {
				r += (i==0?"":"\t") + split[i];
			}
			return r;
		}
		
		HashSet<String> s, e;
		
		public double commonBorders(Prediction o) {
			int anz = 0;
			String[] split;
			if( s== null ) {
				s = new HashSet<String>();
				e = new HashSet<String>();
				for( int i = 0; i < cds.size(); i++ ) {
					split = cds.get(i);
					s.add(split[3]);
					e.add(split[4]);
				}
			}
			for( int i = 0; i < o.cds.size(); i++ ) {
				split = o.cds.get(i);
				
				if( s.contains(split[3]) ) {
					anz++;
				}
				if( e.contains(split[4]) ) {
					anz++;
				}
			}
			return anz;
		}
		
		BitSet nuc;
		
		void fillNuc(int s, int e) {
			nuc = new BitSet(e-s+1);
			for( int i = 0; i < cds.size(); i++ ) {
				String[] split = cds.get(i);
				try {
					nuc.set( Integer.parseInt(split[3])-s, Integer.parseInt(split[4])-s );
				} catch( IndexOutOfBoundsException ioobe ) {
					System.out.println(Arrays.toString(this.split));
					System.out.println(Arrays.toString(split));
					throw ioobe;
				}
			}
		}
		
		public boolean overlap(Prediction o) {
			return nuc.intersects(o.nuc);
		}
		
		public void clear() {
			nuc=null;
		}
	}


	@Override
	public ToolParameterSet getToolParameters() {
		try{
			return
				new ToolParameterSet( getShortName(),
					new SimpleParameter(DataType.STRING,"tag","the tag used to read the GeMoMa annotations",true,GeMoMa.TAG),
					new SimpleParameter(DataType.DOUBLE,"common border filter","the filter on the common borders of transcripts, the lower the more transcripts will be checked as alternative splice isoforms", true, new NumberValidator<Double>(0d, 1d), 0.75 ),
					new SimpleParameter(DataType.INT,"maximal number of transcripts per gene","the maximal number of allowed transcript predictions per gene", true, new NumberValidator<Comparable<Integer>>(1, Integer.MAX_VALUE), Integer.MAX_VALUE ),
					new ParameterSetContainer( "predicted annotation", "", new ExpandableParameterSet( new SimpleParameterSet(		
							new SimpleParameter(DataType.STRING,"prefix","the prefix can be used to distinguish predictions from different input files", false),
							new SimpleParameter(DataType.DOUBLE,"weight","the weight can be used to prioritize predictions from different input files; each prediction will get an additional attribute sumWeight that can be used in the filter", false, new NumberValidator<Double>(0d, 1000d), 1d),
							new FileParameter( "gene annotation file", "GFF files containing the gene annotations (predicted by GeMoMa)", "gff",  true ),
							new FileParameter( "annotation info", "annotation information of the reference, tab-delimted file containing at least the columns transcriptName, GO and .*defline", "tabular",  false )
						), "gene annotations", "", 1 ) ),
					new SimpleParameter(DataType.STRING,"default attributes","Comma-separated list of attributes that is set to NaN if they are not given in the annotation file. "
							+ "This allows to use these attributes for sorting or filter criteria. "
							+ "It is especially meaningful if the gene annotation files were received fom different software packages (e.g., GeMoMa, Braker, ...) having different attributes.", true, "tie,tde,tae,iAA,pAA,score" ),				
					new SimpleParameter(DataType.STRING,"filter","A filter can be applied to transcript predictions to receive only reasonable predictions. The filter is applied to the GFF attributes. "
							+ "The deault filter decides based on the completeness of the prediction (start=='M' and stop=='*') and the relative score (score/aa>=0.75) whether a prediction is used or not. Different criteria can be combined using 'and' and 'or'. In addition, you can check for NaN, e.g., 'isNaN(score)'. "
							+ "A more sophisticated filter could be applied for instance using the completeness, the relative score, the evidence as well as statistics based on RNA-seq: start=='M' and stop =='*' and score/aa>=0.75 and (evidence>1 or tpc==1.0)", false, new RegExpValidator("[a-zA-Z 0-9\\.()><=!-\\+\\*/]*"), "start=='M' and stop =='*' and score/aa>=0.75" ),
					new SimpleParameter(DataType.STRING,"sorting","comma-separated list that defines the way predictions are sorted within a cluster", true, "evidence,score"),
					new SimpleParameter(DataType.STRING,"alternative transcript filter","If a prediction is suggested as an alternative transcript, this additional filter will be applied. The filter works syntactically like the 'filter' parameter and allows you to keep the number of alternative transcripts small based on meaningful criteria. "
							+ "Reasonable filter could be based on RNA-seq data (tie==1) or on evidence (evidence>1). "
							+ "A more sophisticated filter could be applied combining several criteria: tie==1 or evidence>1", false, new RegExpValidator("[a-zA-Z 0-9\\.()><=!-\\+\\*/]*"), "tie==1 or evidence>1" )

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
				"**What it does**\n\nThis tool combines and filters gene predictions from different sources yielding a common gene prediction."
				+ "The tool does not modify the predictions, but filters redundant and low-quality predictions and selects relevant predictions.\n\n"
				
				+ "The algorithm does the following:\n"
				+ "First, redundant predictions are identified (and additional attributes (evidence, sumWeight) are introduced).\n"
				+ "Second, the predictions are filtered using the user-specified criterium based on the attributes from the annotation.\n"
				+ "Third, clusters of overlapping predictions are determined, the predictions are sorted within the cluster and relevant predictions are extracted.\n\n"
				
				+ "Optionally, annotation info can be added for each reference organism enabling a functional prediction for predicted transcripts based on the function of the reference transcript.\n"
				+ "Phytozome provides annotation info tables for plants, but annotation info can be used from any source as long as they are tab-delimited files with at least the following columns: transcriptName, GO and .*defline.\n\n"
				
				+ "Initially, GAF was build to combine gene predictions obtained from GeMoMa. It allows to combine the predictions from multiple reference organisms, but works also using only one reference organism. "
				+ "However, GAF also allows to integrate predictions from ab-initio or purely RNA-seq-based gene predictors as well as manually curated annotation. "
				+ "If you like to do so, we recommend to run AnnotationEvidence for each of these input files to add additional attributes that can be used for sorting and filtering within GAF. "
				+ "The sort and filter criteria need to be carefully revised in this case.\n\n"
				
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