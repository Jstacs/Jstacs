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

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import de.jstacs.DataType;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.FileExistsValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.IntList;
import de.jstacs.utils.SafeOutputStream;
import projects.gemoma.Tools.Ambiguity;

/**
 * This class extracts the information need to run GeMoMa.
 * 
 * @author Jens Keilwagen
 * 
 * @see GeMoMa
 */
public class Extractor extends GeMoMaModule {
	
	private int maxSize;
	
	public Extractor( int maxSize ) {
		this.maxSize = maxSize;
	}
	
	private static void getOut( String prefix, List<File> file, List<SafeOutputStream> out, String temp ) throws IOException {
		File f = prefix == null ? null : Tools.createTempFile("Extractor-" + prefix, temp);
		BufferedOutputStream b = (f == null) ? null : new BufferedOutputStream( new FileOutputStream( f ) );
		file.add(f);
		out.add(SafeOutputStream.getSafeOutputStream(b));
	}

	static String[] name = {"cds-parts", "assignment", "proteins", "cds", "genomic", "introns","identical"};
	static String[] type;
	static {
		type = new String[name.length];
		for( int i = 0; i < name.length; i++ ) {
			type[i] = (i!=1?"fasta":"tabular");
		};
		type[5]="gff";
		type[6]="tabular";
	}
	
	private int[] problem = new int[10];
	private int repair = 0;
	private boolean rep;
	StringBuffer shortInfo = new StringBuffer(), discarded = new StringBuffer(); 
	
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp) throws Exception {
		shortInfo.delete(0, shortInfo.length());
		discarded.delete(0, discarded.length());
		progress.setIndeterminate();
		
		HashMap<String, String> selected = null;
		BufferedReader r;
		String line, comment=null;
		Parameter p = parameters.getParameterForName("selected"); 
		if( p.isSet() ) {
			selected = Tools.getSelection( p.getValue().toString(), maxSize, protocol );
			protocol.append("selected: " + selected.size() + "\t"+ selected+"\n");
		}
		
		HashMap<String, HashMap<String,Gene>> annot = read( (boolean) parameters.getParameterForName("upcase IDs").getValue(), parameters.getParameterForName("annotation").getValue().toString(), selected, protocol );
		
		InputStream in = Tools.getInputStream( parameters.getParameterForName("genetic code"), "projects/gemoma/test_data/genetic_code.txt" );
		HashMap<String,Character> code = Tools.getCode( in );
		
		Ambiguity ambi = (Ambiguity) parameters.getParameterForName("Ambiguity").getValue();
		boolean discardPreMatureStop = (Boolean) parameters.getParameterForName("discard pre-mature stop").getValue();
		boolean stopCodonEx = (Boolean) parameters.getParameterForName("stop-codon excluded from CDS").getValue();
		boolean fullLength = (Boolean) parameters.getParameterForName("full-length").getValue();
		boolean verbose = (Boolean) parameters.getParameterForName("verbose").getValue();
		boolean longComment = (Boolean) getParameter(parameters, "long fasta comment").getValue();
		rep = (Boolean) parameters.getParameterForName("repair").getValue();
		
		p = parameters.getParameterForName("replace");
		String replacement = (p!=null && p.isSet()) ? p.getValue().toString() : null;

		HashSet<String> geneIDs = new HashSet<String>(10000);
		
		ArrayList<File> file = new ArrayList<File>();
		out = new ArrayList<SafeOutputStream>();
		getOut( name[0], file, out, temp );
		getOut( name[1], file, out, temp );
		for( int i = 2; i < name.length; i++ ) {
			p = parameters.getParameterForName(name[i]);
			getOut( (p!=null && ((Boolean)p.getValue())) ? name[i] : null, file, out, temp );
		}
	
		out.get(1).writeln("#geneID\ttranscript\tcds-parts\tphases\tchr\tstrand\tstart\tend\tfull-length\tlongest intron\tsmallest exon\tsplit AA" );
		out.get(6).writeln("#used transcript\tdiscarded transcript" );
		
		//read genome contig by contig
		r = Tools.openGzOrPlain( parameters.getParameterForName("genome").getValue().toString() );
		
		StringBuffer seq = new StringBuffer();
		donor = new HashMap<String, int[]>();
		acceptor = new HashMap<String, int[]>();
		count = new HashMap<Integer, int[]>();
		Arrays.fill(problem, 0);
		int[] info = new int[3];
		Arrays.fill(info, 0);
		
		HashSet<String> unUsedChr = new HashSet<String>( annot.keySet() );
		String ref = "";
		while( (line=r.readLine()) != null ) {
			if( line.startsWith(">") ) {
				//do
				extract( stopCodonEx, fullLength, ambi, protocol, verbose, comment, info, seq, annot, code, discardPreMatureStop, longComment, replacement, geneIDs );
				unUsedChr.remove(comment);
				//clear
				comment = line.substring(1).replaceAll("\\s.*", "");
				ref += (ref.length()>0 ? ", " : "") + comment;
				seq.delete(0, seq.length());
			} else {
				//add
				seq.append( line.trim().toUpperCase() );
			}
		}
		//do
		extract( stopCodonEx, fullLength, ambi, protocol, verbose, comment, info, seq, annot, code, discardPreMatureStop, longComment, replacement, geneIDs );
		unUsedChr.remove(comment);
		r.close();

		ArrayList<TextResult> res = new ArrayList<TextResult>();
		for( int i = 0; i < file.size(); i++ ) {
			File current = file.get(i);
			if( current != null ) {
				out.get(i).close();
				res.add( new TextResult(name[i], "Result", new FileParameter.FileRepresentation(current.getAbsolutePath()), type[i], getToolName(), null, true) );
				//current.deleteOnExit();
			}
		}
		
		shortInfo.append( "\ngenes\t" + info[0] +"\n");
		shortInfo.append( "identical CDS of same gene\t" + info[1] +"\n");
		shortInfo.append( "transcripts\t" + info[2]+"\n\n");
		
		shortInfo.append( "reasons for discarding transcripts:\n");
		shortInfo.append( "ambiguous nucleotide\t" + problem[0] +"\n");
		shortInfo.append( "start phase not zero\t" + problem[1]+"\n");
		shortInfo.append( "missing start\t" + problem[2]+"\n");
		shortInfo.append( "missing stop\t" + problem[3] +"\n");
		shortInfo.append( "premature stop\t" + problem[4]+"\n");
		shortInfo.append( "no DNA\t" + problem[5]+"\n");
		shortInfo.append( "wrong phase\t" + problem[6]+"\n");
		shortInfo.append( "conflicting phase\t" + problem[7]+"\n");
		shortInfo.append( "not linear\t" + problem[8]+"\n\n");
		shortInfo.append( "unexpected error\t" + problem[9]+"\n\n");
		
		shortInfo.append( "repaired\t" + repair+"\n\n");

		if( discarded.length()>0 ) {
			shortInfo.append( "discarded transcript IDs: " + discarded + "\n\n");		
		}
		if( unUsedChr.size() > 0 ) {
			shortInfo.append( "WARNING: There are gene annotations on chromosomes/contigs with missing reference sequence: " + unUsedChr + "\n");
			shortInfo.append( "reference sequence IDS: " + ref + "\n");
		}
		shortInfo.append("\n");
		
		protocol.append( shortInfo.toString() );
		
		//log
		Integer[] array = new Integer[count.size()];
		if( array.length > 0 ) {
			protocol.append("coding exons\t#\n");
			count.keySet().toArray(array);
			Arrays.sort(array);
			for( int i = 0; i < array.length; i++ ) {
				protocol.append( array[i] + "\t" + count.get(array[i])[0] + "\n");
			}
		}
		
		Iterator<Entry<String,int[]>> it;
		Entry<String,int[]> e;
		if( acceptor.size() > 0 ) {
			protocol.append("\nacceptor\t#\n");
			it = acceptor.entrySet().iterator();
			while( it.hasNext() ) {
				e = it.next();
				protocol.append( "acceptor\t" + e.getKey() + "\t" + e.getValue()[0] + "\t" + (e.getValue()[0]/(double)a) + "\n");
			}
		}
		if( donor.size() > 0 ) {
			protocol.append("\ndonor\t#\n");
			it = donor.entrySet().iterator();
			while( it.hasNext() ) {
				e = it.next();
				protocol.append( "donor\t" + e.getKey() + "\t" + e.getValue()[0] + "\t" + (e.getValue()[0]/(double)d) + "\n");
			}
		}
		
		if( intronL != null && intronL.size() > 0 ) {
			protocol.append("\nintron length\t#\tcumulative\n");
			Integer[] il = new Integer[intronL.size()];
			intronL.keySet().toArray(il);
			Arrays.sort(il);
			double all = 0;
			for( int j = 0; j < il.length; j++ ) {
				int[] stat = intronL.get(il[j]);
				all += stat[0];
				protocol.append(il[j] + "\t" + stat[0] + "\t" + (all/anz) + "\n");
			}
		}
		return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
	}
	
	private HashMap<Integer,int[]> intronL = new HashMap<Integer, int[]>();
	private long anz = 0;

	private static String par = "Parent=";
	private static String gID = "gene_id \"";
	private static String tID = "transcript_id \"";
	
	static void add( HashSet<String> errors, HashMap<String, Gene> trans, HashMap<String, HashMap<String,Gene>> annot, String evidence, String c, String strand, String geneID, String transcriptID, String attributes ) {
		HashMap<String,Gene> chr = annot.get(c);
		if( chr == null ) {
			chr = new HashMap<String,Gene>();
			annot.put(c, chr);
		}
		
		Gene gene = chr.get(geneID);			
		if( gene == null ) {
			gene = new Gene(evidence,geneID,strand);
			chr.put(geneID, gene);
			
		} else {
			if( gene.strand != (strand.charAt(0)=='+' ? 1: -1) ) {
				errors.add(c+";"+geneID);
			}
		}
		gene.add( transcriptID, attributes );
		trans.put(transcriptID, gene);
	}
	
	//gff has to be sorted
	/**
	 * This method reads an annotation file and create a data structure to be used in a GeMoMa module.
	 * 
	 * @param upcaseIDs	whether IDs in the input should be upcased
	 * @param input the annotation file either GFF or GTF
	 * @param selected the user-specified selected transcript IDS, if <code>null</code> all transcripts are used  
	 * @param protocol the protocol for reporting
	 * 
	 * @return a {@link HashMap} with contigs/chromosome IDs as keys and another {@link HashMap} as values.
	 * 		The values {@link HashMap} used gene ID as keys and {@link Gene} as values.
	 * 
	 * @throws IOException if the input can not be read
	 */
	public static HashMap<String, HashMap<String,Gene>> read( boolean upcaseIDs, String input, HashMap<String,String> selected, Protocol protocol ) throws IOException {
		HashMap<String, HashMap<String,Gene>> annot = new HashMap<String, HashMap<String,Gene>>();
		Gene gene = null;
		BufferedReader r;
		String line;
		String[] split;
		int idx, h;
		
		//read transcripts
		r = Tools.openGzOrPlain(input);
		HashMap<String, Gene> trans = new HashMap<String, Gene>();
		HashMap<String, ArrayList<String[]>> additional = new HashMap<String, ArrayList<String[]>>();
		ArrayList<String[]> cds = new ArrayList<String[]>();
		boolean first = true, gff = true;
		HashSet<String> errors = new HashSet<String>();
		while( (line=r.readLine()) != null ) {
			if( gff && line.equalsIgnoreCase("##FASTA") ) {//http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
				protocol.append("Stop reading the annotation file because of '##FASTA'\n"); 
				break;  
			}
			if( line.length() == 0 || line.startsWith("#") ) continue;
			if( !gff ) {
				int index = line.indexOf('#');
				if( index > 0 ) {
					line = line.substring(0, index);
				}
			}
			
			split = line.split("\t");
			if( split.length!=9 ) {
				throw new IllegalArgumentException("This line does not seem to be a valid (tab-delimited) line of the annotation with 9 columns: " + line);
			}
			
			if( first ) {
				gff = split[8].indexOf('=')>0;//!( split[8].indexOf(tID)>=0 && split[8].indexOf(gID)>=0 );
				protocol.append("detected annotation format: " + (gff?"GFF":"GTF") + "\n");
				first = false;
			}
			/*
			idx = line.indexOf('\t')+1;
			idx = line.indexOf('\t',idx)+1;
			end = line.indexOf('\t',idx); 
			
			t = line.substring(idx,end);
			switch( t ) {
			*/
			switch( split[2] ) {
				case "CDS":
					if( gff ) {
						boolean add=true;
						if( split[8].indexOf(par)<0 ) {//new version 1.3.3
							idx = split[8].indexOf("ID=");
							if( idx>= 0 ) {
								split[8] = split[8].replace("ID=", par);
								
								idx = split[8].indexOf(par) + par.length();
								h = split[8].indexOf(';',idx);
								
								/*
								String transcriptID = split[8].substring(idx, h>0?h:split[8].length() );
								String geneID = transcriptID+".gene";
								
								//System.out.println(Arrays.toString(split));
								
								add(trans, annot, split[0], split[6], geneID, transcriptID);*/
							} else {
								add=false;
							}
						}
						if( add ) {
							cds.add(split);
						}
					} else {
						idx = split[8].indexOf(tID) + tID.length();
						String transcriptID = split[8].substring(idx, split[8].indexOf('\"',idx));
						String tr = transcriptID;
						if( upcaseIDs ) tr = tr.toUpperCase();
						split[8] = split[8].replace(tID+transcriptID+"\"", par+tr);
						transcriptID=tr;
						
						idx = split[8].indexOf(gID) + gID.length();
						String geneID = split[8].substring(idx, split[8].indexOf('\"',idx));
						if( geneID.length() == 0 ) {
							geneID = transcriptID+".gene";
						}
						add(errors,trans, annot, split[1], split[0], split[6], geneID, transcriptID, null);
						cds.add(split);
					}
					break;
				case "mRNA": case "transcript": case "prediction":
					if( gff ) {
						idx = split[8].indexOf("ID=")+3;
						h = split[8].indexOf(';',idx);
						String transcriptID = split[8].substring(idx, h>0?h:split[8].length() );
						if( upcaseIDs ) transcriptID=transcriptID.toUpperCase();
						if( selected == null || selected.containsKey(transcriptID) ) {
							idx = split[8].indexOf(par);
							if( idx>=0 ) {
								idx+=par.length();
								h = split[8].indexOf(';',idx);
							}
							String geneID = idx<0 ? transcriptID+".gene" : split[8].substring(idx, h>0?h:split[8].length() );
							if( geneID.indexOf(',')>= 0 ) {
								protocol.appendWarning("Could not parse line (multiple parents): " + line + "\n" );
							}
							
							add(errors,trans, annot, split[1], split[0], split[6], geneID, transcriptID, split[8]);
						}
					}
					break;
				default:
					if( gff ) {
						idx = split[8].indexOf(par);
						if( idx>=0 ) {
							idx+=par.length();
							h = split[8].indexOf(';',idx);
							String parentID = split[8].substring( idx, h>0?h:split[8].length() );
							if( upcaseIDs ) parentID = parentID.toUpperCase();
							ArrayList<String[]> l = additional.get(parentID);
							if( l==null ) {
								l = new ArrayList<String[]>();
								additional.put( parentID, l );
							}
							l.add(split);
						}
					}
			}
		}
		r.close();
		
		//read cds
		protocol.append("number of detected CDS lines: " + cds.size() + "\n");
		HashSet<Gene> usedG = new HashSet<Gene>();
		HashSet<String> usedT = new HashSet<String>();
		for( int i = 0 ; i < cds.size(); i++ ) {
			split = cds.get(i);
			//split = line.split("\t");
				
			idx = split[8].indexOf(par)+par.length();
			h = split[8].indexOf(';',idx);
			String[] parent = split[8].substring(idx, h>0?h:split[8].length() ).split(",");
			for( int j = 0; j < parent.length; j++ ) {
				if( upcaseIDs ) parent[j] = parent[j].toUpperCase();
				gene = trans.get(parent[j]);
				if( selected==null || selected.containsKey(parent[j]) ) {
					if( gene == null  ) {
						add(errors, trans, annot, split[1], split[0], split[6], parent[j]+".gene", parent[j], null);
						gene = trans.get(parent[j]);
					}
					usedG.add(gene);
					usedT.add(parent[j]);
					int s = split[6].charAt(0)=='+'?1:-1;
					if( gene.strand != s ) errors.add(split[0]+";"+gene.id);
					gene.add( parent[j], new int[]{
							s, //strand
							Integer.parseInt( split[3] ), //start
							Integer.parseInt( split[4] ), //end
							split[7].charAt(0)=='.' ? Part.NO_PHASE : Integer.parseInt(split[7]) //phase
					} ); 
				}
			}			
		}
		
		//add additional
		if( additional.size()>0 ) {
			Iterator<Gene> it = usedG.iterator();
			while( it.hasNext() ) {
				Gene g = it.next();
				Iterator<String> itT = g.transcript.keySet().iterator();
				while( itT.hasNext() ) {
					String tid = itT.next(); 
					ArrayList<String[]> add = additional.get(tid);
					if( add != null ) {
						Transcript t = g.transcript.get(tid);
						t.add=add;
					}
				}
			}
		}
		
		if( errors.size()>0 ) {//remove errors from add-method
			protocol.append("number of genes with errors: " + errors.size() + "\n");
			Iterator<String> e = errors.iterator();
			while( e.hasNext() ) {
				String er = e.next();
				split = er.split(";");
				HashMap<String,Gene> x = annot.get(split[0]);
				Gene g = x.remove(split[1]);
				
				usedG.remove(g);
				usedT.removeAll(g.transcript.keySet());
			}
		}
		
		protocol.append("number of detected genes: " + usedG.size() + "\n");
		protocol.append("number of detected transcripts: " + usedT.size() + "\n");
		return annot;
	}
	
	public static class Transcript {
		IntList b;
		String attributes;
		ArrayList<String[]> add;
		
		public Transcript( String attributes ) {
			this.attributes = attributes;
			b = new IntList();
			add=null;
		}
	}
	
	/**
	 * This class represents a gene.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class Gene implements Comparable<Gene>{
		/**
		 * The {@link HashMap} contains all transcripts for this gene. 
		 * The keys are the transcript IDs and the values are {@link Transcript}. Each {@link Transcript} contains an IntList of the indices of the (coding) exons this specific transcript.  
		 */
		HashMap<String,Transcript> transcript;
		ArrayList<String> attributes;
		/**
		 * The list of (coding) exons of this gene. The indices of the exon included in this list are used in {@link #transcript}.
		 */
		ArrayList<int[]> exon;
		int start, end;
		/**
		 * The strand of the gene
		 */
		int strand;
		String id, evidence;
		
		public String toString() {
			return id + ": " + transcript.size() + " transcripts";
		}
		
		Gene(String evidence, String id, String strand) {
			this.evidence=evidence;
			transcript = new HashMap<String, Transcript>();
			exon = new ArrayList<int[]>();
			this.id = id;
			this.strand = strand.charAt(0)=='+' ? 1: -1;
			start = end = -1;
			attributes = new ArrayList<String>();
		}
		
/*		Gene(String id, String start, String end, String strand) {
			this( id, strand );
			this.start = Integer.parseInt(start);
			this.end = Integer.parseInt(end);
		}*/
		
		void add( String t, String attributes ) {
			Transcript tr = transcript.get(t);
			if( tr == null ) {
				transcript.put(t,new Transcript(attributes));
			}
		}
	
		void add( String t, int[] border ) {
			IntList x = transcript.get(t).b;
			int i = 0, j;
			while( i < exon.size() ) {
				int[] c = exon.get(i); 
				j = 0;
				while( j < c.length && c[j] == border[j] ) {
					j++;
				}
				if( j < c.length ) {
					i++;
				} else {
					break;
				}
			}
			if( i == exon.size() ) {
				//current exon has not been seen before
				exon.add(border);
			}
			x.add(i);
		}
		
		void reduce( int[] info, SafeOutputStream out ) throws IOException {
			info[0]++;
			String[] s = new String[transcript.size()];
			transcript.keySet().toArray(s);
			Arrays.sort(s);
			boolean[] in = new boolean[s.length];
			Arrays.fill( in, true );
			for( int i = 0; i < s.length; i++ ) {
				Transcript t = transcript.get( s[i] );
				if( t!= null ) {
					IntList il = t.b;
					for( int j = i+1; in[i] && j < s.length; j++ ) {
						if( in[j] && identical(transcript.get(s[j]).b,il) ) {
							in[j] = false;
							out.writeln(s[i] + "\t" + s[j] );
							transcript.remove(s[j]);
							info[1]++;
						}
					}
				}
			}
			sortExons();
		}
		
		void sortExons() {
			Iterator<Entry<String,Transcript>> it = transcript.entrySet().iterator();
			Entry<String,Transcript> e;
			IntList il;
			int[] ids, current, compare;
			boolean swapped = false, swap;
			while( it.hasNext() ) {
				e=it.next();
				il = e.getValue().b;
				ids = il.toArray();
				for( int i = 0; i < ids.length; i++ ) {
					current = exon.get(ids[i]);
					int j=i-1;
					//bubble-up
					while( j >= 0 ) {
						compare = exon.get(ids[j]);
						swap = current[0]*current[1] < compare[0]*compare[1];
						if( swap ) {
							swapped = true;
							int help = ids[j];
							ids[j] = ids[j+1];
							ids[j+1] = help;
						} else {
							break;
						}
						j--;
					}
					
				}
				if( swapped ) {
					il.clear();
					for( int i = 0; i < ids.length; i++ ) {
						il.add(ids[i]);
					}
				}
			}
		}

		void precompute() {
			start = Integer.MAX_VALUE;
			end = Integer.MIN_VALUE;
			for( int i = 0; i < exon.size(); i++ ) {
				int[] current = exon.get(i);
				if ( start > current[1] ) {
					start = current[1];
				}
				if ( end < current[2] ) {
					end = current[2];
				}
			}
			Iterator<Transcript> it = transcript.values().iterator();
			while( it.hasNext() ) {
				Transcript t = it.next();
				if( t.add != null ) {
					for( int j = 0; j<t.add.size(); j++ ) {
						String[] s = t.add.get(j);
						start = Math.min(start, Integer.parseInt(s[3]) );
						end = Math.max(end, Integer.parseInt(s[4]));
					}
				}
			}
		}
		
		@Override
		public int compareTo(Gene o) {
			if( start == -1 ) {
				precompute();
			}
			if( o.start == -1 ) {
				o.precompute();
			}
			return Integer.compare( start+(end-start)/2, o.start+(o.end-o.start)/2 );
		}
	}
	
	private static boolean identical( IntList il1, IntList il2 ) {
		if( il1.length() != il2.length() ) {
			return false;
		}
		for( int i = 0; i < il1.length(); i++ ) {
			if( il1.get(i) != il2.get(i) ) {
				return false;
			}
		}
		return true;
	}	
	
	private void extract( boolean stopCodonEx, boolean fullLength, Ambiguity ambi, Protocol protocol, boolean verbose, String comment, int[] info,
			StringBuffer seq, HashMap<String, HashMap<String,Gene>> annot, HashMap<String,Character> code, boolean discardPreMatureStop, boolean longComment,
			String replacement, HashSet<String> geneIDs
			) throws Exception {
		if( comment == null ) {
			return;
		}
		//out.get(1).write("#"+comment+"\n");
		
		int idx = comment.indexOf(' ');
		String chr = idx>0 ? comment.substring(0,idx) : comment;
		HashMap<String,Gene> chrAnnot = annot.get(chr);
		if( chrAnnot == null ) {
			return;
		}
		//log
		//log.writeln(chr + "\t" + chrAnnot.size() + "\t" + seq.length() );
	
		//all += chrAnnot.size();
		ArrayList<Gene> genes = new ArrayList<Gene>( chrAnnot.values() );
		Collections.sort(genes);

		int max = 5000;
		boolean[] used = new boolean[max];
		boolean[] donS = new boolean[max];
		boolean[] accS = new boolean[max];
		String[] don = new String[max];
		String[] acc = new String[max];
		SafeOutputStream ident = out.get(6);
		HashMap<String, int[]> introns = new HashMap<String,int[]>();
		for( Gene gene: genes ) {
			if( gene.transcript.size()>0 ) {
				if( gene.id.split("\\s").length>1 ) {
					if( replacement!= null ) {
						String newID = gene.id.replaceAll("\\s", replacement);
						protocol.append("replaced gene-ID: " + gene.id + " by " + newID + "\n");
						gene.id=newID;
					} else {
						throw new IllegalArgumentException( "Gene ID ("+gene.id+") contains a whitespace character, please use parameter replace" );
					}
				}
				if( geneIDs.contains( gene.id ) ) {
					throw new IllegalArgumentException( "At least two genes with gene id: " + gene.id );
				} else {
					geneIDs.add(gene.id);
				}
				int[] val = null;
				boolean[][] splits = new boolean[gene.exon.size()][gene.exon.size()];
				for( int k = 0; k < splits.length; k++ ) {
					Arrays.fill( splits[k], false );
				}
				
				int strand = gene.strand;
				boolean forward=strand==1;
				
				gene.reduce( info, ident );
				part.clear();
				int i, j;
				
				String[] id = new String[gene.transcript.size()];
				gene.transcript.keySet().toArray(id);
				Arrays.sort(id);
				Arrays.fill( accS, false );
				Arrays.fill( donS, false );
				for( int k = 0; k < id.length; k++ ) {
					IntList il = gene.transcript.get( id[k] ).b;
					for( j = 0; j < il.length(); j++ ) {
						i = il.get(j);
						if( j != 0 ) {
							accS[i]=true;
						}
						if( j+1<il.length() ) {
							donS[i]=true;
						}
					}
				}
				
				for( i = 0; i < gene.exon.size(); i++ ) {
					val = gene.exon.get(i);
					
					int off1 = val[1]-1>=2 ? 2 : 0;
					int off2 = val[2]+2<=seq.length() ? 2 : 0;
					String p, s;
					try {
						p = seq.substring( val[1]-1-off1, val[2]+off2 );
						//check
						if( strand < 0 ) {
							p=Tools.rc(p);
						}
						
						if( strand > 0 ) {
							s=p.substring(off1,p.length()-off2);
							acc[i] = off1 > 0 ? p.substring(0,off1) : "";
							don[i] = off2 > 0 ? p.substring(p.length()-off2,p.length()) : "";
						} else {
							s=p.substring(off2,p.length()-off1);
							acc[i] = off2 > 0 ? p.substring(0,off2) : "";
							don[i] = off1 > 0 ? p.substring(p.length()-off1,p.length()) : "";
						}
					} catch( StringIndexOutOfBoundsException sioobe ) {
						s=null;
					}
					part.add(new Part(s,val));
				}
				
				Arrays.fill( used, false );
				for( int k = 0; k < id.length; k++ ) {
					IntList il = gene.transcript.get( id[k] ).b;
					boolean[] set = new boolean[il.length()];
					for( j = 0; j < il.length(); j++ ) {
						Part current = part.get( il.get(j) );
						set[j] = current.aa != null;
					}
					int prob = transcript( seq, stopCodonEx, chr, gene, id[k], -1, splits, fullLength, info, ambi, code, protocol, verbose, used, acc, don, discardPreMatureStop, longComment );
					
					//try to repair
					if( prob>=0 && rep) {
						int phase = -1, test;
						do {
							phase++;
							//clear CONDITIONALLY
							for( j = 0; j < il.length(); j++ ) {
								Part current = part.get( il.get(j) );
								if( !set[j] ) {
									current.offsetLeft = Part.NO_PHASE;
									current.aa = null;
								}
							}
							test = transcript( seq, stopCodonEx, chr, gene, id[k], phase, splits, fullLength, info, ambi, code, protocol, false, used, acc, don, discardPreMatureStop, longComment );
						} while( test >= 0 && phase <= 2 );
						if( test < 0 ) {
							if( verbose ) protocol.appendWarning(id[k] + "\trepaired with start phase " + phase + "\n" );
							repair++;
							prob=-1;
						} else {
							if( prob>= 0 ) {
								//clear CONDITIONALLY							
								for( j = 0; j < il.length(); j++ ) {
									Part current = part.get( il.get(j) );
									if( !set[j] ) {
										current.offsetLeft = Part.NO_PHASE;
										current.aa = null;
									}
								}
							}
						}
					}
					
					if( prob >= 0 ) {
						problem[prob]++;
						discarded.append( (discarded.length()>0?", ":"") + id[k] );
					} else {
						//intron length
						int last=-1;
						for( j = 0; j < il.length(); j++ ) {
							Part current = part.get( il.get(j) );
							if( last != -1 ) {
								int inLe = forward ? current.start-last : last - current.end;
								int[] num = intronL.get(inLe);
								if( num == null ) {
									num=new int[1];
									intronL.put(inLe, num);
								}
								num[0]++;
								anz++;
							}
							last=forward ? current.end : current.start;
						}
					}
				}
				
				//write cds parts, ...
				for( j = 0; j < gene.exon.size(); j++ ) {
					if( used[j] ) {
						Part p = part.get(j);
						if( p.aa.length() > 0 ) {
							out.get(0).write(">" + gene.id + "_" + j + "\n" + p.aa + "\n");
						}
						/*if( p.acc != null ) {
							out.get(4).write(">" + gene.id + "_" + j + "\n" + p.acc + "\n");
						}
						if( p.don != null ) {
							out.get(5).write(">" + gene.id + "_" + j + "\n" + p.don + "\n");
						}*/
						
						if( !out.get(5).doesNothing() ) {
							for( int k = 0; k < splits.length; k++ ) {
								if( splits[j][k] ) {
									int st, en;
									if( forward ) {
										st = gene.exon.get(j)[2]+1;
										en = gene.exon.get(k)[1];
									} else {
										st = gene.exon.get(k)[2]+1;
										en = gene.exon.get(j)[1];
									}
									String key = chr + "\tannotation\tintron\t" + st + "\t" + en + "\t\t"+ (forward?"+":"-") + "\t.\t.";
									int[] counts = introns.get(key); 
									if( counts == null ) {
										counts = new int[1];
										introns.put(key, counts);
									}
									counts[0]++;
									/*
									String intr = seq.substring(st-1, en-1);
									if( !forward ) {
										intr = Tools.rc(intr);
									}
									System.out.println( Arrays.toString( gene.exon.get(j) ) );
									System.out.println( Arrays.toString( gene.exon.get(k) ) );
									System.out.println(forward + "\t" + intr);
									*/
									//out.get(5).writeln(chr + "\tannotation\tintron\t" + st + "\t" + en + "\t.\t" + (forward?"+":"-") + "\t.\t." );
								}
							}
						}
					}
				}
			}
		}

		if( !out.get(5).doesNothing() ) {
			//write introns
			String[] keys = introns.keySet().toArray(new String[0]);
			Arrays.sort(keys);//sub optimal sorting due to string instead of integers for positions
			for( int i = 0;i < keys.length; i++ ) {
				int[] count = introns.get(keys[i]);
				out.get(5).writeln(keys[i].replaceFirst("\t\t", "\t" + count[0] + "\t") );
			}
		}
		
	}
	
	ArrayList<SafeOutputStream> out;
	ArrayList<Part> part = new ArrayList<Part>();
	StringBuffer dnaSeqBuff = new StringBuffer();
	IntList message = new IntList();
	HashMap<String, int[]> donor, acceptor;
	HashMap<Integer,int[]> count;
	int a=0, d=0;
	
	int transcript(StringBuffer seq, boolean stopCodonEx, String chr, Gene gene, String trans, int s, boolean[][] splits, boolean fullLength, int[] info, Ambiguity ambi, HashMap<String,Character> code, Protocol protocol,boolean verbose, boolean[] used, String[] acc, String[] don, boolean discardPreMatureStop, boolean longComment ) throws IOException {
		int j;
		dnaSeqBuff.delete(0, dnaSeqBuff.length());
		int currentProb=-1;
		
		Transcript tr = gene.transcript.get( trans ); 
		IntList il = tr.b;
		if( il.length() == 0 ) {
			protocol.append("No coding exon(s) for: " + trans + "\n");
			return -1;
		}
		int start = gene.strand>0 ? gene.exon.get(il.get(0))[1] : gene.exon.get(il.get(il.length()-1))[1];
		int end = gene.strand>0 ? gene.exon.get(il.get(il.length()-1))[2] : gene.exon.get(il.get(0))[2];
		
		int startPhase = (s>=0 && s<3) ? s : part.get(il.get(0)).offsetLeft;
		if( startPhase == Part.NO_PHASE ) {
			startPhase = 0;
		}
		int offset = 3-startPhase, pa = -1;
		message.clear();
		Part current = null;
		int minExon=Integer.MAX_VALUE, maxIntron=-1, lastPos=-1;
		for( j = 0; j < il.length(); j++ ) {
			pa = il.get(j);
			current = part.get(pa);
			if( lastPos>=0 ) {
				if( (gene.strand>0 && current.start <= lastPos)
					|| (gene.strand<0 && current.end >= lastPos) ) {
					if( verbose ) {
						int oldPa = il.get(j-1);
						Part old = part.get(oldPa);
						protocol.appendWarning(trans + "\tnot linear: "+ il.toString() + " strand=" + gene.strand + " " + oldPa + ":" + old.start + ".." + old.end + " " + pa + ":" + current.start + ".." + current.end + "\n" );
					}
					return 8;
				}
			}
			if( gene.strand>0 ) {
				lastPos = current.end;
			} else {
				lastPos = current.start;
			}
			if( current.dna == null ) {
				currentProb=1;
				break;
			}
			dnaSeqBuff.append( current.dna );
			
			//translate
			if( current.offsetLeft == Part.NO_PHASE ) {
				current.offsetLeft = (3-offset)%3;
			}
			
			if( current.aa == null ) {
				try {
					current.aa = Tools.translate(current.offsetLeft, current.dna, code, false, ambi);
				} catch( IllegalArgumentException iae ) {
					current.aa=null;
					currentProb=0;
					break;
				}							
				current.offsetRight = current.dna.length() - current.offsetLeft - 3*current.aa.length();
			} else {
				if( (fullLength || j>0) && current.offsetLeft != (3-offset)%3 ) {
					currentProb=2;
					break;
				}
			}
			
			if( current.aa!= null && current.aa.length()>0 && !current.aa.matches("[A-Za-z" + (discardPreMatureStop?"":"\\*")+ "]*" +(j+1==il.length() && !stopCodonEx?("\\*"+(fullLength?"{1}":"{0,1}")):"")) ) {//since > v1.3.1
				message.add(pa);
			}
			if( current.aa!= null && current.aa.length()<minExon ) {
				minExon=current.aa.length();
			}
			
			if( j>0 ) {
				int l;
				if( gene.strand == 1 ) {
					l = gene.exon.get( il.get(j) )[1] - gene.exon.get( il.get(j-1) )[2]-1;
				} else {
					l = gene.exon.get( il.get(j-1) )[1] - gene.exon.get( il.get(j) )[2]-1;
				}
				
				if( l > maxIntron ) {
					maxIntron=l;
				}
			}
			
			offset = current.offsetRight;
		}
		
		if( stopCodonEx && current != null ) {
			if( current.aa == null || current.aa.length()==0 ) {
				current.aa="*";
			} else if( current.aa.charAt(current.aa.length()-1)!='*' ) {
				current.aa +="*";
			}
		}

		String p=null;
		if( j == il.length() ) {
			if( ambi == Ambiguity.EXCEPTION && !dnaSeqBuff.toString().matches("[ACGT]*") ) {//to be compatible with older versions
				j=il.length()+1;
				currentProb=0;
			} else {
				try {
					p = Tools.translate(startPhase, dnaSeqBuff.toString(), code, false, ambi);
					if( stopCodonEx ) {
						p += "*";
					}
				} catch( IllegalArgumentException iae ) {
					j=il.length()+1;
					currentProb=0;
				}
			}
		}
		
		if( j == il.length() ) {
			if( p.length() == 0 ) {
				return 9;
			}
			
			int idx = p.indexOf('*');
			boolean preMature = 0<=idx && idx<p.length()-1;
			
			if( fullLength && startPhase!= 0 ) {
				if( verbose ) protocol.appendWarning(trans + "\tskip start phase not zero\n" );
				return 1;
			} else if( fullLength && p.charAt(0)!='M' ) {
				if( verbose ) writeWarning(protocol, trans, il.length(), "skip missing start",p, dnaSeqBuff );
				return 2;
			} else if( fullLength && p.charAt(p.length()-1)!='*' ){
				if( verbose ) writeWarning(protocol, trans, il.length(), "skip missing stop",p, dnaSeqBuff );
				return 3;
			} else if( preMature && discardPreMatureStop ) {
				if( verbose ) writeWarning(protocol, trans, il.length(), "skip premature stop",p, dnaSeqBuff );
				return 4;
			} else if( message.length() > 0 ) {
				if( verbose ) {
					String c = "";
					for( j = 0; j < il.length(); j++ ) {
						pa = il.get(j);
						current = part.get(pa);
						c+="cds-parts: " + pa + " (phase: " + current.offsetLeft + ")\nDNA: " + current.dna +"\nAA: " +current.aa +"\n";
					}
					
					protocol.appendWarning( trans + "\tskip wrong phase for coding part(s) = "+message+"\n" +
						"\nCDS: " + dnaSeqBuff.toString() +
						"\nprotein: " + p +
						"\n\nparts:\n" + c ); //since > v1.3.1
				}
				return 6;
			} else {				
				info[2]++;
				String comment = trans;
				if( longComment ) {
					comment += " gene=" + gene.id + " chr=" + chr + " strand=" + gene.strand + " interval=" + start + ".." + end;
					if( tr.attributes != null ) {
						String a = deleteAttribute( "Parent=", tr.attributes );
						a = deleteAttribute( "ID=", a );
						a=a.replaceAll(";", " ");
						comment += " " + a;
					}
				}
				out.get(3).write( ">" + comment + "\n" + dnaSeqBuff.toString() + "\n" );
				out.get(2).write( ">" + comment + "\n" + p + "\n" );
				String x = il.toString();
				SafeOutputStream sos = out.get(1);
				sos.write( gene.id + "\t" + trans + "\t" + x.substring(1,x.length()-1).replaceAll(" ", "") );
				String splitAA = "";
				int currentPos = 0;
				for( j = 0; j < il.length(); j++ ) {
					current = part.get(il.get(j));
					sos.write( (j==0?"\t":",") + current.offsetLeft );
					
					int pos = currentPos / 3;					
					int oldPos=currentPos;
					currentPos += current.dna.length();
					if( j>1 ) {
						splitAA += ",";
					}
					if( oldPos % 3 > 0 && pos < (currentPos/3) ) {
						splitAA += p.substring(pos, pos+1);
					}
				}
				sos.write( "\t" + chr + "\t" + gene.strand + "\t" + start + "\t" + end + "\t" + (p.charAt(0)=='M' && p.charAt(p.length()-1)=='*') + "\t" + (maxIntron>0?maxIntron:"NA") + "\t" + minExon + "\t" + splitAA + "\n" );
				for( j = 0; j < il.length(); j++ ) {
					used[il.get(j)] = true;
				}
															
				int[] c = count.get(il.length());
				if( c == null ) {
					c = new int[1];
					count.put(il.length(), c);
				}
				c[0]++;
				
				//genomic region
				if( !out.get(4).doesNothing() ) {
					StringBuffer genomicRegion = new StringBuffer();
					int off = 300, st;
					
					Part firstP = part.get(il.get(0));
					Part lastP = part.get(il.get(il.length()-1));
					boolean forward = gene.strand==1;
					if( forward ) {
						st = Math.max(firstP.start-off,1);
						genomicRegion = new StringBuffer( seq.substring(st-1, Math.min( seq.length(), lastP.end+off)).toLowerCase() );
					} else {
						st = Math.min(seq.length(), firstP.end+off);
						genomicRegion = new StringBuffer( Tools.rc(seq.substring(Math.max(lastP.start-off,1)-1, 
								st )).toLowerCase() );
					}
	
					for( j = 0; j < il.length(); j++ ) {
						Part t = part.get(il.get(j));
						int a = t.start-st, b = t.end-st;
						if( forward ) {
							genomicRegion.replace(a, b+1, genomicRegion.substring(a, b+1).toUpperCase() );
						} else {
							genomicRegion.replace(-b, -a+1, genomicRegion.substring(-b, -a+1).toUpperCase() );
						}
					}
					out.get(4).writeln(">" + comment );
					out.get(4).writeln(genomicRegion);
				}
				
				//splice sites
				int ignoreAAForSpliceSite=30, targetStart, targetEnd;
				int last=-1;
				for( j = 0; j < il.length(); j++ ) {
					pa = il.get(j);
					current = part.get(pa);
					if( last != -1 ) {
						splits[last][pa] = true;
					}
					last = pa;
					
					//find splice sites
					int[] exon = gene.exon.get(pa);
					if( gene.strand == 1 ) {
						targetStart = exon[1] + current.offsetLeft;
						targetEnd = exon[2] - current.offsetRight;
					} else {
						targetStart = exon[1] + current.offsetRight;
						targetEnd = exon[2] - current.offsetLeft;
					}
					
					int l = Math.abs(targetStart-1-targetEnd), add, t = 3*ignoreAAForSpliceSite;
					
					//add = the length which is used inside the exon to find a splice site
					if( l / 3 < t ) {
						add=l/3;
						add-=add%3; //;)
					} else {
						add=t;
					}
					
					if( current.dna.length()>0 ) {
						if( j > 0 && acc[pa].length() > 0 ) {
							int[] stat = acceptor.get(acc[pa]);
							if( stat == null ) {
								stat = new int[1];
								acceptor.put(acc[pa], stat);
							}
							stat[0]++;
							a++;
						}
						if( j+1<il.length() && don[pa].length() > 0 ) {
							int[] stat = donor.get(don[pa]);
							if( stat == null ) {
								stat = new int[1];
								donor.put(don[pa], stat);
							}
							stat[0]++;
							d++;
						}
					}
				}
				return -1;
			}
		} else {
			switch ( currentProb ) {
				case 0:
					if( verbose ) {
						if( j < il.length() ) {
							protocol.appendWarning(trans + "\tskip non-ACGT coding part "+j+"\n");
						} else {
							protocol.appendWarning(trans + "\tskip non-ACGT coding protein\n");
						}
					}
					return 0;
				case 1:
					if( verbose ) protocol.appendWarning(trans + "\tskip no DNA for coding part "+j+"\n");
					return 5;
				case 2:
					if( verbose ) protocol.appendWarning(trans + "\tskip conflicting phase for coding part "+j+"\n");
					return 7;
					
			}
			return 8;	
		}
	}	
	
	String deleteAttribute( String attribute, String all ) {
		int startIndex = all.indexOf(attribute);
		if( startIndex>=0 ) {
			int endIndex = all.indexOf(';', startIndex);
			if( endIndex < 0 ) {
				all=all.substring(0,startIndex);
			} else {
				all=all.substring(0,startIndex) + all.substring(endIndex+1);
			}
		}
		return all;
	}
	
	void writeWarning( Protocol protocol, String trans, int numExons, String reason, String p, CharSequence dna ) {
		protocol.appendWarning(trans + "\tnumExons=" + numExons + "\t" + reason + "\n"+p+"\n"+dna+"\n" );
	}
	
	static class Part {
		static final int NO_PHASE = -100000;
		
		String dna, aa;
		int offsetLeft, offsetRight, start, end;
		
		Part( String dna, int[] val ) {
			this.dna= dna;
			aa = null;
			start = val[1];
			end = val[2];
			offsetLeft = val[3];
			offsetRight = NO_PHASE;
		}
	}
	
	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getShortName(),
				new FileParameter( "annotation", "Reference annotation file (GFF or GTF), which contains gene models annotated in the reference genome", "gff,gff3,gtf,gff.gz,gff3.gz,gtf.gz", true, new FileExistsValidator(), true ),
				new FileParameter( "genome", "Reference genome file (FASTA)", "fasta,fa,fas,fna,fasta.gz,fa.gz,fas.gz,fna.gz",  true, new FileExistsValidator(), true ),

				new FileParameter( "genetic code", "optional user-specified genetic code", "tabular", false ),
					
				new SimpleParameter(DataType.BOOLEAN, Extractor.name[2], "whether the complete proteins sequences should returned as output", true, false ),
				new SimpleParameter(DataType.BOOLEAN, Extractor.name[3], "whether the complete CDSs should returned as output", true, false ),
				new SimpleParameter(DataType.BOOLEAN, Extractor.name[4], "whether the genomic regions should be returned (upper case = coding, lower case = non coding)", true, false ),
				new SimpleParameter(DataType.BOOLEAN, Extractor.name[5], "whether introns should be extracted from annotation, that might be used for test cases", true, false ),
				new SimpleParameter(DataType.BOOLEAN, Extractor.name[6], "if CDS is identical Extractor only used one transcript. This parameter allows to return a table that lists in the first column the used transcript and in the second column the discarded transcript. If no transcript is discarded, the list is empty.", true, false ),
				new SimpleParameter(DataType.BOOLEAN, "upcase IDs", "whether the IDs in the GFF should be upcased", true, false ),
				new SimpleParameter(DataType.BOOLEAN, "repair", "if a transcript annotation can not be parsed, the program will try to infer the phase of the CDS parts to repair the annotation", true, false ),
				
				/*
				new SelectionParameter(DataType.PARAMETERSET, new String[]{"no","yes"}, new ParameterSet[]{
						new SimpleParameterSet(),
						new SimpleParameterSet(
								new SimpleParameter( DataType.INT, "intronic", "The number of bp return from the intron side", true, new NumberValidator<Integer>(2,1000), 10 ),
								new SimpleParameter( DataType.INT, "exonic", "The number of bp return from the exon side", true, new NumberValidator<Integer>(0,1000), 8 )
								//negative?
						)
					}, "splice sites", "whether splice sites should be returned or not", true ),
				/**/	
				new FileParameter( "selected", "The path to list file, which allows to make only a predictions for the contained transcript ids. The first column should contain transcript IDs as given in the annotation. Remaining columns will be ignored.", "tabular,txt", maxSize>-1, new FileExistsValidator() ),
				new EnumParameter( Ambiguity.class, "This parameter defines how to deal with ambiguities in the DNA. There are 3 options: "
						+ "EXCEPTION, which will remove the corresponding transcript, "
						+ "AMBIGUOUS, which will use an X for the corresponding amino acid, and "
						+ "RANDOM, which will randomly select an amnio acid from the list of possibilities.", true, Ambiguity.EXCEPTION.toString() ),
				new SimpleParameter( DataType.BOOLEAN, "discard pre-mature stop", "if *true* transcripts with pre-mature stop codon are discarded as they often indicate misannotation", true, true ),
				new SimpleParameter( DataType.BOOLEAN, "stop-codon excluded from CDS", "A flag that states whether the reference annotation contains the stop codon in the CDS annotation or not", true, false ),
				new SimpleParameter( DataType.BOOLEAN, "full-length", "A flag which allows for choosing between only full-length and all (i.e., full-length and partial) transcripts", true, true ),
				new SimpleParameter( DataType.BOOLEAN, "long fasta comment", "whether a short (transcript ID) or a long (transcript ID, gene ID, chromosome, strand, interval) fasta comment should be written for proteins, CDSs, and genomic regions", true, false ),
				new SimpleParameter( DataType.BOOLEAN, "replace", "Replacement for whitespace characters in gene IDs", false ),
				new SimpleParameter( DataType.BOOLEAN, "verbose", "A flag which allows to output a wealth of additional information", true, false )
			);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public String getToolName() {
		return "Extractor";
	}

	@Override
	public String getShortName() {
		return getToolName();
	}

	@Override
	public String getDescription() {
		return "extracts parts of CDSs as annotated in a genome (assembly)";
	}

	@Override
	public String getHelpText() {
		return "This tool can be used to create input files for **GeMoMa**, i.e., it creates at least a fasta file containing the translated parts of the CDS and a tabular file containing the assignment of transcripts to genes and parts of CDS to transcripts."
				+ " In addition, **Extractor** can be used to create several additional files from the final prediction, e.g. proteins, CDSs, ... ."
				+ " Two inputs are mandatory: The genome as fasta or fasta.gz and the corresponding annotation as gff or gff.gz."
				+ " The gff file should be sorted."
				+ " If you like to set a user-specific genetic code, please use a tab-delimited file with two columns. The first column contains the amino acid in one letter code, the second a list of tripletts."
				+ MORE;
	}
	
	static final String EXAMPLE = " Here is an example\n\n"
			+ "+---+------------------------------+\n"
			+ "| I | ATT, ATC, ATA                |\n"
			+ "+---+------------------------------+\n"
			+ "| L | CTT, CTC, CTA, CTG, TTA, TTG |\n"
			+ "+---+------------------------------+\n"
			+ "| V | GTT, GTC, GTA, GTG           |\n"
			+ "+---+------------------------------+\n"
			+ "|...| ...                          |\n"
			+ "+---+------------------------------+";
	
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		ResultEntry[] re = new ResultEntry[2];
		for( int i = 0; i < re.length; i++ ) {
			re[i] = new ResultEntry(TextResult.class, type[i], name[i]);
		}
		return re;
	}
	
	@Override
	public ToolResult[] getTestCases( String path ) {
		try {
			return new ToolResult[]{new ToolResult(FileManager.readFile(path+File.separator+"tests/gemoma/xml/extractor-test.xml"))};
		} catch( Exception e ) {
			e.printStackTrace();
			return null;
		}
	}
}