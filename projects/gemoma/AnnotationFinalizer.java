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
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.FileExistsValidator;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.parameters.validation.RegExpValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;

/**
 * This class finalizes an annotation. It predicts UTRs for a given set of CDS annotations using RNA-seq evidence and allows to generate generic gene and transcript names.
 * 
 * @author Jens Keilwagen
 */
public class AnnotationFinalizer extends GeMoMaModule {
	
	static boolean rename;
	
	static int set( Feature[] array, ArrayList<? extends Feature> list, int start ) {
		if( list == null ) return start;
		for( int i = 0; i < list.size(); i++, start++ ) {
			array[start]=list.get(i);
		}
		return start;
	}
	
	static class FeatureComparator implements Comparator<Feature> {

		static final FeatureComparator DUMMY = new FeatureComparator();
		
		@Override
		public int compare(Feature first, Feature other) {
			int diff = first.start-other.start;
			if( diff==0 ) {
				diff=(first.end-first.start) - (other.end-other.start);
			}
			if( first.split[6].charAt(0)=='-' ) {
				diff*=-1;
			}
			return diff;
		}
	}
	
	static abstract class Feature{
		int start, end;
		String[] split;
		ArrayList<Feature> add;
		
		Feature( String[] split ) {
			this.split=split;
			this.start=Integer.parseInt(split[3]);
			this.end=Integer.parseInt(split[4]);
		}
		
		public abstract Feature addFeature( String[] split );
		
		public void add( String[] split ) {
			if( add == null ) {
				add = new ArrayList<Feature>();
			}
			add.add(new AddFeature(split));
		}
		
		public void rename( String tag, String oldN, String newN ) {
			String sep = split[8].indexOf(';')>=0 ? ";": "";
			split[8]=split[8].replace(tag+"="+oldN+sep, tag+"="+newN+sep);
		}

		static<T extends Feature> void rename( String tag, String oldN, String newN, ArrayList<T> l ) {
			if( l==null ) return;
			for( int i = 0; i < l.size(); i++ ) {
				l.get(i).rename(tag, oldN, newN);
			}
		}
		
		public String toString() {
			return toString(add!=null);
		}
		
		public String toString( boolean addAdd ) {
			StringBuffer all = new StringBuffer();
			for( int i = 0; i < split.length; i++ ) {
				all.append( split[i] + (i+1==split.length?"\n":"\t") );
			}
			if( addAdd ) {
				for( int i = 0; i < add.size(); i++ ) {
					all.append(add.get(i).toString());
				}
			}
			return all.toString();
		}

	}
	
	//additional feature
	static class AddFeature extends Feature {

		AddFeature(String[] split) {
			super(split);
		}

		@Override
		public Feature addFeature(String[] split) {
			throw new RuntimeException("IllegalOperation");
		}
	}
	
	
	static class Gene extends Feature implements Comparable<Gene> {
		ArrayList<Transcript> list;
		int min=-100, max=-100;
		
		Gene( String[] split ) {
			super( split );
			list = new ArrayList<Transcript>();
		}
		
		public int strand() {
			return split[6].charAt(0)=='+' ? 1 : -1;
		}
		
		public Transcript addFeature( String[] split ){
			Transcript t = new Transcript(split);
			list.add(t);
			return t;
		}

		void determine() {
			min = Integer.MAX_VALUE;
			for( int i = 0; i < list.size(); i++ ) {
				min = Math.min( list.get(i).determine(), min);
			}
			max = Integer.MIN_VALUE;
			for( int i = 0; i < list.size(); i++ ) {
				max = Math.max( Integer.parseInt(list.get(i).split[4]), max);
			}			
		}
		
		void set() {
			determine();
			split[3] = "" + min;
			split[4] = "" + max;
			//not needed: Collections.sort(list);
		}
		
		@Override
		public int compareTo(Gene g) {
			if( min == -100 ) {
				determine();
				//Collections.sort(list);
			}
			if( g.min == -100 ) {
				g.determine();
				//Collections.sort(g.list);
			}
			return Integer.compare( min, g.min );
		}
		
		public String toString() {
			StringBuffer all = new StringBuffer();
			all.append(super.toString());
			for( int i = 0; i < list.size(); i++ ) {
				all.append(list.get(i).toString());
			}
			return all.toString();
		}

		public void genericName(String name, boolean asName) {
			if( asName ) {
				//simple: new attribute "Name"
				split[8]="Name=" + name + ";" + split[8];
				for( int i = 0; i < list.size(); i++ ) {
					Transcript t = list.get(i);
					t.split[8]="Name=" + name + "." + (i+1)+";" + t.split[8];
				}
			} else {
				//complex: ID, Parent
				int i1 = split[8].indexOf("ID=");
				int i2 = split[8].indexOf(";",i1);
				String old = split[8].substring(i1+3,i2);
				split[8]=split[8].replace("ID="+old+";", "ID="+name+";");
				Feature.rename("Parent", old, name, add);
				for( int i = 0; i < list.size(); i++ ) {
					Transcript t = list.get(i);
					
					i1 = t.split[8].indexOf("ID=");
					i2 = t.split[8].indexOf(";",i1);
					String oldT = t.split[8].substring(i1+3,i2);
					String newT = name + "." + (i+1);
					
					t.rename("Parent",old,name);
					t.rename("ID",oldT,newT);
					
					Feature.rename("Parent", oldT, newT, t.upUTR);
					Feature.rename("Parent", oldT, newT, t.list);
					Feature.rename("Parent", oldT, newT, t.downUTR);
					Feature.rename("Parent", oldT, newT, t.add);
				}
			}
		}
	}
	
	static class Transcript extends Feature implements Comparable<Transcript> {
		ArrayList<CDS> list;
		ArrayList<UTR> upUTR;
		ArrayList<UTR> downUTR;
		String id;
		double score;
		
		Transcript( String[] split ) {
			super( split );
			String[] attr = split[8].split(";");
			int l = 0;
			while( l < attr.length && !attr[l].startsWith("ID=") ) {
				l++;
			}
			id=attr[l].substring(3);
			
			l = 0;
			while( l < attr.length && !attr[l].startsWith("score=") ) {
				l++;
			}
			if( l<attr.length ) {
				try {
					score = Double.parseDouble(attr[l].substring(6));
				} catch( Exception e ) {
					score = Double.NaN;
				}
			} else {
				score = 0;
			}
			
			list = new ArrayList<CDS>();
			upUTR = new ArrayList<UTR>();
			downUTR = new ArrayList<UTR>();
		}
		
		public void add( String[] split ) {
			switch( split[2] ) {
				case "three_prime_UTR": downUTR.add(new UTR(split) ); break;
				case "five_prime_UTR": upUTR.add(new UTR(split) ); break;
				default: super.add(split);
			}
		}
		
		public int strand() {
			return split[6].charAt(0)=='+' ? 1 : -1;
		}

		public int determine() {
			int min = getMin();
			int max = getMax();
			split[3] = "" + min;
			split[4] = "" + max;
			return min;
		}
		
		public int getMin() {
			int min = Integer.MAX_VALUE;
			if( add!= null ) for( int i = 0; i < add.size(); i++ ) {
				min = Math.min(min, add.get(i).start );
			}

			for( int i = 0; i < upUTR.size(); i++ ) {
				min = Math.min(min, upUTR.get(i).getMin() );
			}
			for( int i = 0; i < list.size(); i++ ) {
				min = Math.min(min, list.get(i).getMin() );
			}
			for( int i = 0; i < downUTR.size(); i++ ) {
				min = Math.min(min, downUTR.get(i).getMin() );
			}
			return min;
		}
		
		public int getMax() {
			int max = Integer.MIN_VALUE;
			if( add!=null ) for( int i = 0; i < add.size(); i++ ) {
				max = Math.max(max, add.get(i).end );
			}
			
			for( int i = 0; i < upUTR.size(); i++ ) {
				max = Math.max(max, upUTR.get(i).getMax() );
			}
			for( int i = 0; i < list.size(); i++ ) {
				max = Math.max(max, list.get(i).getMax() );
			}
			for( int i = 0; i < downUTR.size(); i++ ) {
				max = Math.max(max, downUTR.get(i).getMax() );
			}
			return max;
		}

		@Override
		public CDS addFeature(String[] split) {
			CDS c = new CDS(split);
			list.add(c);
			return c;
		}
		
		public void addUTR( boolean up, String[] split ) {
			UTR utr = new UTR(split);
			if( up ) {
				upUTR.add( utr );
			} else {
				downUTR.add( utr );
			}
		}
		
		public String toString() {
			StringBuffer all = new StringBuffer();
			all.append(super.toString(false));

			int num = (add==null?0:add.size()) + upUTR.size() + list.size() + downUTR.size();
			Feature[] array = new Feature[num];
			int start = 0;
			start = set(array, add, start);
			start = set(array, upUTR, start);
			start = set(array, list, start);
			start = set(array, downUTR, start);
			Arrays.sort(array,FeatureComparator.DUMMY);
			
			for( int j = 0; j < array.length; j++ ) {
				all.append(array[j].toString());
			}
			return all.toString();
		}

		@Override
		public int compareTo(Transcript t) {
			return -Double.compare(score, t.score);
		}
	}
	
	static class CDS extends Feature {
		
		CDS( String[] split ) {
			super( split );
		}

		@Override
		public Feature addFeature(String[] split) {
			throw new RuntimeException("not supported");
		}
		
		public int getMin() {
			return Integer.parseInt(split[3]);
		}
		
		public int getMax() {
			return Integer.parseInt(split[4]);
		}
	}
	
	static class UTR extends Feature {
		UTR( String[] split ) {
			super( split );
		}
		
		@Override
		public Feature addFeature(String[] split) {
			throw new RuntimeException("not supported");
		}
		
		public int getMin() {
			return Integer.parseInt(split[3]);
		}
		
		public int getMax() {
			return Integer.parseInt(split[4]);
		}
	}
	
	static HashMap<String, HashMap<String,Gene>> read( String fName, Protocol protocol, String tag, StringBuffer begin, boolean add ) throws IOException {
		String[] tags = {"gene", tag, "CDS"};//GeMoMa specific
		//read annotation
		BufferedReader r = new BufferedReader( new FileReader( fName ) );
		ArrayList<String[]>[] list = new ArrayList[tags.length + (add?1:0)];
		HashMap<String,Integer> tagIndex = new HashMap<String, Integer>();
		for( int i = 0; i < tags.length; i++ ) {
			list[i] = new ArrayList<String[]>();
			tagIndex.put( tags[i],  i );
		}
		if( add ) list[tags.length] = new ArrayList<String[]>();
		
		String line;
		boolean gff=true;//TODO currently only for gff
		String[] split;
		boolean start = true;
		while( (line=r.readLine()) != null ) {
			if( gff && line.equalsIgnoreCase("##FASTA") ) {//http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
				//protocol.append("Stop reading the annotation file because of '##FASTA'\n"); 
				break;  
			}
			if( line.length() == 0 ) continue;
			if( line.startsWith("#") ) {
				if( start && !line.startsWith("##sequence-region")) {
					begin.append(line+"\n");
				}
				continue;
			}
			start=false;
			if( !gff ) {
				int index = line.indexOf('#');
				if( index > 0 ) {
					line = line.substring(0, index);
				}
			}
			
			split = line.split("\t");
			Integer i = tagIndex.get( split[2] );
			if( i==null && add ) {
				i=list.length-1;
			}
			if( i!= null ) {
				list[i].add(split);
			}
		}		
		r.close();
		
		//parse annotation
		HashMap<String,HashMap<String,Gene>> res = new HashMap<String, HashMap<String,Gene>>();
		HashMap<String, Feature> all = new HashMap<String,Feature>();
		//int old=0;
		int[] w = new int[2];
		for( int i = 0; i < tags.length; i++ ) {
			int anz = 0;
			for( int j=0; j < list[i].size(); j++ ) {
				split = list[i].get(j);
				
				if( i == 0 ) {
					String chr = split[0];
					HashMap<String, Gene> h = res.get(chr);
					if( h == null ) {
						h = new HashMap<String, Gene>();
						res.put(chr, h);
					}
					String[] attr = split[8].split(";");
					int k = 0;
					while( k < attr.length && !attr[k].startsWith("ID=") ) {
						k++;
					}
					if( k < attr.length ) {
						String id = attr[k].substring(3);
						Gene g = new Gene(split);
						anz++;
						h.put( id, g );
						all.put( id, g );
					} else {
						System.out.println("WARNING1: " + Arrays.toString(split) );
						w[0]++;
					}
				} else {
					String[] attr = split[8].split(";");
					int k = 0;
					while( k < attr.length && !attr[k].startsWith("Parent=") ) {
						k++;
					}
					int l = 0;
					while( l < attr.length && !attr[l].startsWith("ID=") ) {
						l++;
					}
					if( k < attr.length ) {
						String parent = attr[k].substring(7);
						Feature f = all.get(parent);
						Feature n = f.addFeature(split);
						anz++;
						if( l < attr.length ) {
							String id = attr[l].substring(3);
							
							//System.out.println(id + "\t" + parent );
							
							all.put( id, n );
						}
					} 
					
					if( k==attr.length || (i+1!=tags.length && l==attr.length) ) {
						System.out.println("WARNING2: " + Arrays.toString(split) );
						w[1]++;
					}
				}
			}
			protocol.append( "#" + tags[i] + "s: " + anz + "\n" );
			protocol.append( "#warnings: " + Arrays.toString(w) + "\n" );
			//old=all.size();
		}
		
		if( add ) {
			for( int j=0; j < list[tags.length].size(); j++ ) {
				split = list[tags.length].get(j);
				String parent = get(PARENT,split);
				if( parent!=null ) {
					//has Parent feature
					Feature f = all.get(parent);
					if( f==null ) {
						throw new NullPointerException("Could not find Parent: " + Arrays.toString(split) );
					} else {
						f.add(split);
					}
				} else {
					throw new NullPointerException("Could not find Parent attribute: " + Arrays.toString(split) );
				}
			}
		}
		
		return res;
	}
	
	
	static final String PARENT = "Parent=";
	
	static String get( String tag, String[] split ) {
		String[] attr = split[8].split(";");
		int k = 0;
		while( k < attr.length && !attr[k].startsWith(tag) ) {
			k++;
		}
		if( k < attr.length ) {
			//has tag
			return attr[k].substring(tag.length());
		} else {
			return null;
		}
	}
	
	static class SequenceIDComparator implements Comparator<String> {
	    public int compare(String o1, String o2) {
	        int diff = Long.compare(extractLong(o1),extractLong(o2));
	        if( diff == 0 ) {
	        	diff = o1.compareTo(o2);
	        }
	        return diff;
	    }

	    long extractLong(String s) {
	        String num = s.replaceAll("\\D", "");
	        // return 0 if no digits found
	        if( num.isEmpty() ) {
	        	return 0;
	        } else {
	        	try {
	        		long l = Long.parseLong(num);
	        		return l;
	        	} catch (NumberFormatException nfe) {
	        		return 0;
	        	}
	        }
	    }
	}
	
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp ) throws Exception {
		return run(parameters, protocol, progress, threads, null, null, temp);
	}
	
	ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, ToolParameterSet description, String version, String temp ) throws Exception {
		progress.setIndeterminate();
		SimpleParameterSet utrPs = (SimpleParameterSet) parameters.getParameterForName("UTR").getValue();
		if( GeMoMa.seqs == null  ) {
			String genome = (String) parameters.getParameterForName("genome").getValue();

			int reads = 1;
			ExpandableParameterSet introns=null, coverage=null;
			
			if( utrPs.getNumberOfParameters()>0 ) {
				reads = (Integer) utrPs.getParameterForName("reads").getValue();
				introns = (ExpandableParameterSet)((ParameterSetContainer)utrPs.getParameterAt(0)).getValue();
				coverage = (ExpandableParameterSet)((ParameterSetContainer)utrPs.getParameterAt(2)).getValue();
			}
			GeMoMa.fill(protocol, false, 0, genome, null, reads, introns, coverage );				
			protocol.append("\n");
		}
		
		String source=null;
		if( utrPs.getNumberOfParameters()>0 ) {
			Parameter p = utrPs.getParameterForName("additional source suffix");
			source = p==null ? null : (String) p.getValue();
			
			if( source!=null && source.length()==1) {
				source=null;
			}
		}
		source= toolName+(source==null?"":("_"+source));

		SimpleParameterSet renamePS = (SimpleParameterSet) parameters.getParameterForName("rename").getValue();
		int n = renamePS.getNumberOfParameters();
		int digits=-100;
		String prefix, infix=null, suffix = "";
		String search = null, replace = null;
		rename = n>0;
		if( rename ) {
			prefix = renamePS.getParameterForName("prefix").getValue().toString();
			digits = (Integer) renamePS.getParameterForName("digits").getValue();
			if( n>2 ) {
				infix = renamePS.getParameterForName("infix").getValue().toString();
				suffix = renamePS.getParameterForName("suffix").getValue().toString();
				search = renamePS.getParameterForName("contig search pattern").getValue().toString();
				replace = renamePS.getParameterForName("contig replace pattern").getValue().toString();
			}
		} else {
			prefix=null;
		}
		boolean asName = (Boolean) parameters.getParameterForName("name attribute").getValue();
		
		Parameter p = parameters.getParameterForName("transfer features");
		boolean addAdditional=p==null?false:(Boolean)p.getValue();
		
		//annotation
		String tag = parameters.getParameterForName("tag").getValue().toString();
		StringBuffer begin = new StringBuffer();
		HashMap<String, HashMap<String,Gene>> annotation = 
				//Extractor.read( parameters.getParameterForName("annotation").getValue().toString(), null, protocol);
				read( parameters.getParameterForName("annotation").getValue().toString(), protocol, tag, begin, addAdditional );

		
		//iterate over all chromosomes/contigs
		File gffFile = Tools.createTempFile(getShortName(),temp);
		BufferedWriter w = new BufferedWriter(new FileWriter(gffFile));
		String info;
		if( description == null ) {
			w.append(begin);
			w.append(INFO + getShortName() + " " + getToolVersion() + "; ");
			info = JstacsTool.getSimpleParameterInfo(parameters);
		} else {
			w.append("##gff-version 3");
			w.newLine();
			w.append(INFO + description.getToolName() + " " + version + "; ");
			info = JstacsTool.getSimpleParameterInfo(description);
		}
		if( info != null ) {
			w.append("SIMPLE PARAMETERS: " + info );
		}
		w.newLine();
		
		int num = 1, utrBoth=0, utr=0, utr5=0, utr3=0;
		String[] keys = new String[annotation.size()];
		annotation.keySet().toArray(keys);
		Arrays.sort( keys, new SequenceIDComparator() );
		for( String c : keys ) {
			w.append("##sequence-region " + c + " 1 " + GeMoMa.seqs.get(c).length() );
			w.newLine();
			
			String newC = c;
			String pref = null;
			if( rename ) {
				if( n==2 ) {
					pref = prefix;
				} else {
					newC=newC.replaceAll(search, replace);
					pref = prefix + newC + infix;
					num=1;
				}
			}
			
			Collection<Gene> gg = annotation.get(c).values();
			ArrayList<Gene> genes = new ArrayList<Gene>(gg);
			Collections.sort(genes);
					
			int[][] fwdCov, revCov;
			if( GeMoMa.coverage !=  null ) { 
				fwdCov = GeMoMa.coverage[0].get(c);
				revCov = GeMoMa.coverage[1].get(c);
			} else {
				fwdCov=revCov=null;
			}
			
			int[][][] acc, don;
			if( GeMoMa.acceptorSites!= null ) {
				acc = GeMoMa.acceptorSites.get(c);
				don = GeMoMa.donorSites.get(c);
			} else {
				acc=don=null;
			}
			
			//System.out.println(c + "\t" + newC + "\t" + fwdCov + "\t" + revCov + "\t" + acc + "\t" + don );
			
			if( fwdCov!= null || revCov != null ) {
				for( Gene g : genes ) {
					//get coverage for current strand
					int strand = g.strand();
					int[][] a = acc==null?null:acc[strand==1?0:1];
					int[][] d = don==null?null:don[strand==1?0:1];
					int[][] cov = strand == 1 ? fwdCov : revCov;
					if( cov != null ) {
						for( Transcript t : g.list ) {
//if( c.equals("1") ) System.out.println(t);
							boolean u5 = extendUTR( t, -1, cov, a, 1, source ); //upstream
							boolean u3 = extendUTR( t, 1, cov, d, 0, source ); //downstream
//if( c.equals("1") ) System.exit(1);
							if( u5 ) {
								utr5++;
							}
							if( u3 ) {
								utr3++;
							}
							if( u5 || u3 ) {
								utr++;
								//System.exit(0);
							}
							if( u5 && u3 ) {
								utrBoth++;
								//System.exit(0);
							}
						}
					}
					g.set();
				}
			}
			
			Collections.sort(genes);
			for( Gene g : genes ) {
				if( rename ) {
					g.genericName( pref + String.format("%0"+digits+"d", num++) + suffix, asName );
				}
				w.append(g.toString());
			}
		}
		w.close();
		
		if( GeMoMa.seqs != null ) {
			protocol.append("\n#transcripts with 5'-UTR annotation: " + utr5 + "\n");
			protocol.append("#transcripts with 3'-UTR annotation: " + utr3 + "\n");
			protocol.append("#transcripts with some UTR annotation: " + utr + "\n");
			protocol.append("#transcripts with 5'- and 3'-UTR annotation: " + utrBoth + "\n");
		}
		ArrayList<TextResult> res = new ArrayList<TextResult>();
		res.add( new TextResult(defResult, "Result", new FileParameter.FileRepresentation(gffFile.getAbsolutePath()), "gff", getToolName(), null, true) );
		
		return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());

	}
	
	private static boolean extendUTR( Transcript t, int d, int[][] cov, int[][] splice, int sIdx, String source ) {
		int strand = t.strand();
		int min = t.getMin();
		int max = t.getMax();
		int dir = d*strand;
				
		String utr=(d==1?"three_prime_UTR":"five_prime_UTR");
		LinkedList<String[]> utrs = new LinkedList<String[]>();
		
		int start = dir==-1 ? min : max;
		boolean extend = d==1 ? t.downUTR.size()==0 : t.upUTR.size()==0; //only extend if there is no annotation available yet (cf. transfer feature) 
//System.out.println(t.id);
		while( extend ) {
			int firstBase = start+dir;
			extend=false;
			
			//find stretch that is has some positive coverage
			int h=firstBase;
			int idx = Arrays.binarySearch(cov, new int[]{firstBase}, GeMoMa.IntArrayComparator.comparator[2] );
			if( idx < 0 ) {
				idx = -(idx+1);
				idx = Math.max(0, idx-1);//has to be reduced since we did not find a match
			}
			if( cov[idx][0] <= firstBase && firstBase <= cov[idx][1] ) {
				//System.out.println(dir);
				do {
					h = cov[idx][dir==-1?0:1];
					idx+=dir;
				} while( idx >= 0 && idx < cov.length && cov[idx][dir==-1?1:0] == h+dir );
			} else {
				break;
			}
			int lastBase = h;
			if( dir==-1 ) {
				int xxx = lastBase;
				lastBase = firstBase;
				firstBase = xxx;
			}
//System.out.println(utr + "\t" + firstBase + "\t" + lastBase + "\t" + dir);
			
			
			//find splice site in interval [firstBase-1, lastBase+1]
			if( splice != null && splice[sIdx]!= null) {
				int st = firstBase+dir;
				idx = Arrays.binarySearch(splice[sIdx], st);
				if( idx > 0 ) {
					while( idx>0 &&  splice[sIdx][idx-1] == st ) {
						idx--;
					}
				}
				if( idx < 0 ) {
					idx = -(idx+1);
					//we don't need this here: idx = Math.max(0, idx-1);//has to be reduced since we did not find a match
				}
				int maxSplitReads = -1, maxIdx=-1;
				int end = lastBase + dir;
//if( idx >= 0 && idx < splice[sIdx][0] ) System.out.println("splice: ["+st+","+end+"]\t"+splice[sIdx][idx] + " <= " + end);
				while( idx < splice[sIdx].length && splice[sIdx][idx] <= end ) {
					int i=Arrays.binarySearch(cov, new int[]{splice[sIdx][idx]+dir}, GeMoMa.IntArrayComparator.comparator[2] );
					int c;
					if( i < 0 ) {				
						i = -(i+1);
						i = Math.max(0, idx-1);//has to be reduced since we did not find a match
						if( cov[i][1] > splice[sIdx][idx]+dir ) {
							c=cov[i][2];
						} else {
							c=0;
						}
					} else {
						c=cov[i][2];
					}

//System.out.println("looking for a splice site (" + sIdx + ") " + firstBase + "\t" + idx + "\t" + splice[0][idx] + ".." + splice[1][idx] + " ("+splice[2][idx]+")\t" + c + "\t" + (splice[2][idx] >= c) + "\t" + (splice[2][idx] > maxSplitReads) );
					
					if( splice[2][idx] >= c //check whether its better to use the split read (=intron) than following the coverage
						&& splice[2][idx] > maxSplitReads //use the intron with the most split reads
					) {
						maxSplitReads = splice[2][idx];
						maxIdx = idx;
					}
					idx++;
				}
				
				if( maxIdx>= 0 ) {
					extend = true;
					int add = (dir==1) ? -1 : 0;//(strand==1?0:-1);
					if( dir == -1 ) {
						//set new start base
						firstBase = splice[sIdx][maxIdx]+add;
					} else {
						//set new last base
						lastBase = splice[sIdx][maxIdx]+add;
					}
					start = splice[1-sIdx][maxIdx]+add;
				}
			}
						
			//add utr;
			if( firstBase <= lastBase ) {
				String[] split = {
					t.split[0],
					source,
					utr,
					""+ firstBase,
					""+ lastBase,
					".",
					t.split[6],
					".",
					null
				};
				
				if( d == -1 ) {
					utrs.addFirst( split );
				} else {
					utrs.addLast( split );
				}
			}
		}
		boolean hasUTR = utrs.size()>0;
		
		//int u = 1;
		while( utrs.size() > 0 ) {
			String[] split;
			split = utrs.removeFirst();
			split[8] = /*"ID=" + t.id + "." + utr + "."+ u++ +";"+*/ "Parent="+t.id; //+";";
			t.addUTR(d==-1, split);
		}
		
		return hasUTR;
	}
	
	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getShortName(),
					new FileParameter( "genome", "The genome file (FASTA), i.e., the target sequences in the blast run. Should be in IUPAC code", "fasta,fa,fas,fna,fasta.gz,fa.gz,fas.gz,fna.gz", true, new FileExistsValidator(), true ),
					new FileParameter( "annotation", "The predicted genome annotation file (GFF)", "gff,gff3", true, new FileExistsValidator(), true ),
					
					new SimpleParameter( DataType.STRING, "tag", "A user-specified tag for transcript predictions in the third column of the returned gff. It might be beneficial to set this to a specific value for some genome browsers.", true, GeMoMa.TAG ),					
					new SimpleParameter( DataType.BOOLEAN, "transfer features", "if true other features than gene, &lt;tag&gt; (default: mRNA), and CDS of the input will be written in the output", true, false ),

					new SelectionParameter( DataType.PARAMETERSET, 
							new String[]{"NO", "YES"},
							new Object[]{
								//no UTRs
								new SimpleParameterSet(),
								//UTR prediction
								new SimpleParameterSet(
									new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
											new FileParameter( "introns file", "Introns (GFF), which might be obtained from RNA-seq", "gff,gff3", true, new FileExistsValidator(), true )
										), "introns", "", 1 ) ),
									new SimpleParameter( DataType.INT, "reads", "if introns are given by a GFF, only use those which have at least this number of supporting split reads", true, new NumberValidator<Integer>(1, Integer.MAX_VALUE), 1 ),
				
									new ParameterSetContainer( "coverage", "", new ExpandableParameterSet( new SimpleParameterSet(	
										new SelectionParameter( DataType.PARAMETERSET, 
												new String[]{"NO", "UNSTRANDED", "STRANDED"},
												new Object[]{
													//no coverage
													new SimpleParameterSet(),
													//unstranded coverage
													new SimpleParameterSet(
															new FileParameter( "coverage_unstranded", "The coverage file contains the unstranded coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true, new FileExistsValidator() )
													),
													//stranded coverage
													new SimpleParameterSet(
															new FileParameter( "coverage_forward", "The coverage file contains the forward coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true, new FileExistsValidator() ),
															new FileParameter( "coverage_reverse", "The coverage file contains the reverse coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true, new FileExistsValidator() )
													)
												},  "coverage file", "experimental coverage (RNA-seq)", true
										)
									), "coverage", "", 1 ) ),
									
									new SimpleParameter( DataType.STRING, "additional source suffix", "a suffix for source values of UTR features", true, new RegExpValidator("\\w*"), "" )
								),
							},
							new String[] {
								"no UTR prediction",
								"UTR prediction"
							},
							"UTR", "allows to predict UTRs using RNA-seq data", true
					),
					
					new SelectionParameter( DataType.PARAMETERSET, 
							new String[]{"COMPOSED", "SIMPLE", "NO"},
							new Object[]{
								//generic renaming
								new SimpleParameterSet(
										new SimpleParameter( DataType.STRING, "prefix", "the prefix of the generic name", true, new RegExpValidator("\\w+") ),
										new SimpleParameter( DataType.STRING, "infix", "the infix of the generic name", true, "G" ),
										new SimpleParameter( DataType.STRING, "suffix", "the suffix of the generic name", true, "0" ),
										new SimpleParameter( DataType.INT, "digits", "the number of informative digits", true, new NumberValidator<Integer>(4, 10), 5 ),
										new SimpleParameter( DataType.STRING, "contig search pattern", "search string, i.e., a regular expression for search-and-replace parts of the contig/scaffold/chromosome names, the modified string is used as infix for the gene name", true, "" ),
										new SimpleParameter( DataType.STRING, "contig replace pattern", "replace string, i.e., a regular expression for search-and-replace parts of the contig/scaffold/chromosome names, the modified string is used as infix for the gene name", true, new RegExpValidator("[^;^=^\"^\\s]*"), "" )
								),
								//simple renaming
								new SimpleParameterSet(
										new SimpleParameter( DataType.STRING, "prefix", "the prefix of the generic name", true, new RegExpValidator("\\w+") ),
										new SimpleParameter( DataType.INT, "digits", "the number of informative digits", true, new NumberValidator<Integer>(4, 10), 5 )
								),
								//no renaming
								new SimpleParameterSet()
							},
							new String[] {
								"composed renaming: <prefix><contig/chr*><infix><number><suffix>; recommended for whole genome (re-)annotations",
								"simple renaming: <prefix><number>; recommended for gene family annotations",
								"no renaming"
							},
							"rename", "allows to generate generic gene and transcripts names (cf. parameter &quot;name attribute&quot;)", true
					),
					new SimpleParameter( DataType.BOOLEAN, "name attribute", "if true the new name is added as new attribute &quot;Name&quot;, otherwise &quot;Parent&quot; and &quot;ID&quot; values are modified accordingly", true, true )
					
			);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	static String toolName = "AnnotationFinalizer";
	
	public String getToolName() {
		return toolName;
	}
		
	public String getShortName() {
		return toolName;
	}

	public String getDescription() {
		return "predicts UTRs for a given set of CDS annotations";
	}

	public String getHelpText() {
		return 
			"This tool finalizes an annotation."
			+ " It allows to predict for UTRs for annotated coding sequences and to generate generic gene and transcript names. UTR prediction might be negatively influenced (i.e. too long predictions) by genomic contamination of RNA-seq libraries, overlapping genes or genes in close proximity as well as unstranded RNA-seq libraries."
			+ " Please use **ERE** to preprocess the mapped reads."
			+ MORE;
	}

	static final String defResult = "final_annotation";
	
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", defResult),
		};
	}
	
	@Override
	public ToolResult[] getTestCases( String path ) {
		try {
			return new ToolResult[]{
				new ToolResult(FileManager.readFile(path+File.separator+"tests/gemoma/xml/af-test-nothing.xml")),
				new ToolResult(FileManager.readFile(path+File.separator+"tests/gemoma/xml/af-test-utr.xml"))
			};
		} catch( Exception e ) {
			e.printStackTrace();
			return null;
		}
	}
}