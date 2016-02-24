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
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import projects.gemoma.Tools.Ambiguity;
import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.IntList;
import de.jstacs.utils.SafeOutputStream;

/**
 * This class extracts the information need to run GeMoMa.
 * 
 * @author Jens Keilwagen
 * 
 * @see GeMoMa
 */
public class Extractor implements JstacsTool {
	
	private int maxSize;
	
	public Extractor( int maxSize ) {
		this.maxSize = maxSize;
	}
	
	private static void getOut( String prefix, List<File> file, List<SafeOutputStream> out ) throws IOException {
		File f = prefix == null ? null : File.createTempFile(prefix,"tmp_GeMoMa");
		BufferedOutputStream b = (f == null) ? null : new BufferedOutputStream( new FileOutputStream( f ) );
		file.add(f);
		out.add(SafeOutputStream.getSafeOutputStream(b));
	}

	private static String[] name = {"cds-parts", "assignment", "proteins", "transcripts", "acceptor", "donor"};
	private static String[] type;
	static {
		type = new String[name.length];
		for( int i = 0; i < name.length; i++ ) {
			type[i] = (i!=1?"fasta":"tabular");
		};
	}
	
	@Override
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		progress.setIndeterminate();
		
		HashMap<String, String> selected = null;
		BufferedReader r;
		String line, comment=null;
		Parameter p = parameters.getParameterForName("selected"); 
		if( p.isSet() ) {
			selected = Tools.getSelection( p.getValue().toString(), maxSize, protocol );
			protocol.append("selected: " + selected.size() + "\t"+ selected+"\n");
		}
		
		HashMap<String, HashMap<String,Gene>> annot = readGFF( parameters.getParameterForName("annotation").getValue().toString(), selected );
		HashMap<Integer,int[]> count = new HashMap<Integer, int[]>();
		int[] problem = new int[4];
		Arrays.fill(problem, 0);
		int[] info = new int[3];
		Arrays.fill(info, 0);
		
		count.clear();
		InputStream in;// = parameters.getParameterForName("genetic code").getValue().toString();
		p = parameters.getParameterForName("genetic code");
		if( p.isSet() ) {
			in = new FileInputStream( p.getValue().toString() );
		} else {
			in = Extractor.class.getClassLoader().getResourceAsStream("projects/gemoma/test_data/genetic_code.txt" );
		}
		HashMap<String,Character> code = Tools.getCode( in );
		
		boolean verbose = (Boolean) parameters.getParameterForName("verbose").getValue();
		
		ArrayList<File> file = new ArrayList<File>();
		ArrayList<SafeOutputStream> out = new ArrayList<SafeOutputStream>();
		getOut( name[0], file, out );
		getOut( name[1], file, out );
		getOut( ((Boolean)parameters.getParameterForName(name[2]).getValue()) ? name[2] : null, file, out );
		getOut( ((Boolean)parameters.getParameterForName(name[3]).getValue()) ? name[3] : null, file, out );
		
		int intronic, exonic;
		exonic = intronic = 0;
		p = parameters.getParameterForName("splice sites");
		boolean splice = p!= null;
		if( splice ) {
			ParameterSet ps = (ParameterSet) p.getValue();
			splice = ps.getNumberOfParameters() > 0;
			if( splice ) {
				exonic = (Integer) ps.getParameterForName("exonic").getValue();
				intronic = (Integer) ps.getParameterForName("intronic").getValue();
				System.out.println(intronic + "\t" + exonic);
			}
		}
		getOut( splice ? name[4] : null, file, out );
		getOut( splice ? name[5] : null, file, out  );
	
		//read genome contig by contig
		r = new BufferedReader( new FileReader( parameters.getParameterForName("genome").getValue().toString() ) );
		
		StringBuffer seq = new StringBuffer();
		while( (line=r.readLine()) != null ) {
			if( line.startsWith(">") ) {
				//do
				extract( protocol, verbose, comment, problem, info, count, intronic, exonic, seq, annot, code, out );
				//clear
				comment = line.substring(1);
				seq.delete(0, seq.length());
			} else {
				//add
				seq.append(line.trim().toUpperCase() );
			}
		}
		//do
		extract( protocol, verbose, comment, problem, info, count, intronic, exonic, seq, annot, code, out );
		r.close();

		
		ArrayList<TextResult> res = new ArrayList<TextResult>();
		for( int i = 0; i < file.size(); i++ ) {
			File current = file.get(i);
			if( current != null ) {
				out.get(i).close();
				res.add( new TextResult(name[i], "Result", new FileParameter.FileRepresentation(current.getAbsolutePath()), type[i], getToolName(), null, true) );
				current.deleteOnExit();
			}
		}
		
		protocol.append( "genes\t" + info[0] +"\n");
		protocol.append( "removed transcripts\t" + info[1] +"\n");
		protocol.append( "transcripts\t" + info[2]+"\n\n");
		
		protocol.append(Arrays.toString(problem)+"\n\n");
		//log
		protocol.append("parts\t#\n");
		Integer[] array = new Integer[count.size()];
		count.keySet().toArray(array);
		Arrays.sort(array);
		for( int i = 0; i < array.length; i++ ) {
			protocol.append( array[i] + "\t" + count.get(array[i])[0] + "\n");
		}
		
		return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
	}

	private static String par = "Parent=";
	private static String geneTag = "gene";
	
	//gff has to be sorted 
	private static HashMap<String, HashMap<String,Gene>> readGFF( String input, HashMap<String,String> selected ) throws Exception {
		HashMap<String, HashMap<String,Gene>> annot = new HashMap<String, HashMap<String,Gene>>();
		HashMap<String,Gene> chr;
		Gene gene;
		BufferedReader r;
		String line, geneID = null, transcriptID, t;
		String[] split;
		int idx, h, end;
		/*//old
		r = new BufferedReader( new FileReader(input) );
		while( (line=r.readLine()) != null && line.startsWith("##") );
		do {
			if( line.length() == 0 ) continue; 
			idx = line.indexOf('\t')+1;
			idx = line.indexOf('\t',idx)+1;
			end = line.indexOf('\t',idx); 
			t = line.substring(idx,end);
			if( t.equals( geneTag ) || t.equals( "CDS" ) ) {
				//System.out.println(line);
				split = line.split("\t");
				if( t.equals( geneTag ) ) {
					idx = split[8].indexOf("ID=")+3;
					h = split[8].indexOf(';',idx);
					geneID = split[8].substring(idx, h>0?h:split[8].length() );
				} else {//CDS
					idx = split[8].indexOf(par)+par.length();
					h = split[8].indexOf(';',idx);
					transcriptID = split[8].substring(idx, h>0?h:split[8].length() );
					
					chr = annot.get(split[0]);
					if( chr == null ) {
						chr = new HashMap<String,Gene>();
						annot.put(split[0], chr);
					}
					gene = chr.get(geneID);
					if( gene == null ) {
						gene = new Gene();
						chr.put(geneID, gene);
					}
					gene.add( transcriptID, new int[]{
							split[6].charAt(0)=='+'?1:-1, //strand
							Integer.parseInt( split[3] ), //start
							Integer.parseInt( split[4] ) //end
					} );
					//System.out.println( geneID + "\t" + transcriptID + "\t" + gene.exon.size() );
				}
			}
		} while( (line=r.readLine()) != null );
		r.close();
		/**/
		
		//read genes
		r = new BufferedReader( new FileReader(input) );
		ArrayList<String> transcript = new ArrayList<String>();
		ArrayList<String> cds = new ArrayList<String>();
		while( (line=r.readLine()) != null && line.startsWith("##") );
		do {
			if( line.length() == 0 ) continue; 
			idx = line.indexOf('\t')+1;
			idx = line.indexOf('\t',idx)+1;
			end = line.indexOf('\t',idx); 
			t = line.substring(idx,end);
			if( t.equalsIgnoreCase( "CDS") ) {
				cds.add(line);
			} else if( t.equalsIgnoreCase( geneTag ) ) {
				//System.out.println(line);
				split = line.split("\t");
				idx = split[8].indexOf("ID=")+3;
				h = split[8].indexOf(';',idx);
				geneID = split[8].substring(idx, h>0?h:split[8].length() );
				chr = annot.get(split[0]);
				if( chr == null ) {
					chr = new HashMap<String,Gene>();
					annot.put(split[0], chr);
				}
				chr.put(geneID, new Gene());
			} else if( t.equalsIgnoreCase( "mRNA" ) || t.equalsIgnoreCase("transcript") ) {
				transcript.add(line);
			}
		} while( (line=r.readLine()) != null );		
		r.close();
		
		//read transcripts
		HashMap<String, Gene> trans = new HashMap<String, Gene>();
		for( int i = 0 ; i < transcript.size(); i++ ) {
			line = transcript.get(i);
			split = line.split("\t");
			
			idx = split[8].indexOf("ID=")+3;
			h = split[8].indexOf(';',idx);
			transcriptID = split[8].substring(idx, h>0?h:split[8].length() ).toUpperCase();
			if( selected == null || selected.containsKey(transcriptID) ) {
				idx = split[8].indexOf(par)+par.length();
				h = split[8].indexOf(';',idx);
				geneID = split[8].substring(idx, h>0?h:split[8].length() );
				
				gene = annot.get(split[0]).get(geneID);
				gene.add( transcriptID );
				
				trans.put(transcriptID, gene);
			}
		}
		
		//read cds
		for( int i = 0 ; i < cds.size(); i++ ) {
			line = cds.get(i);
			split = line.split("\t");
				
			idx = split[8].indexOf(par)+par.length();
			h = split[8].indexOf(';',idx);
			transcriptID = split[8].substring(idx, h>0?h:split[8].length() ).toUpperCase();
			
			if( selected == null || selected.containsKey(transcriptID) ) {
				gene = trans.get(transcriptID);
				gene.add( transcriptID, new int[]{
						split[6].charAt(0)=='+'?1:-1, //strand
						Integer.parseInt( split[3] ), //start
						Integer.parseInt( split[4] ) //end
				} );
			}
		}
		return annot;
	}
	
	private static class Gene {
		HashMap<String,IntList> transcript;
		ArrayList<int[]> exon;
		
		Gene() {
			transcript = new HashMap<String, IntList>();
			exon = new ArrayList<int[]>();
		}
		
		void add( String t ) {
			transcript.put(t,new IntList());
		}
	
		void add( String t, int[] border ) {
			IntList x = transcript.get(t);
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
		
		void reduce( String geneName, int[] info ) throws IOException {
			info[0]++;
			String[] s = new String[transcript.size()];
			transcript.keySet().toArray(s);
			Arrays.sort(s);
			boolean[] in = new boolean[s.length];
			Arrays.fill( in, true );
			for( int i = 0; i < s.length; i++ ) {
				IntList il = transcript.get( s[i] );
				for( int j = i+1; in[i] && j < s.length; j++ ) {
					if( in[j] && identical(transcript.get(s[j]),il) ) {
						in[j] = false;
						//out.append( geneName + "\t" + s[j] + "\t \n" );
						transcript.remove(s[j]);
						info[1]++;
					}
				}
			}
			sortExons();
		}
		
		void sortExons() {
			Iterator<Entry<String,IntList>> it = transcript.entrySet().iterator();
			Entry<String,IntList> e;
			IntList il;
			int[] ids, current, compare;
			boolean swapped = false, swap;
			while( it.hasNext() ) {
				e=it.next();
				il = e.getValue();
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
	
	private static void extract( Protocol protocol, boolean verbose, String comment, int[] problem, int[] info, HashMap<Integer,int[]> count,
			int intronic, int exonic, StringBuffer seq, HashMap<String, HashMap<String,Gene>> annot, HashMap<String,Character> code,
			ArrayList<SafeOutputStream> out) throws Exception {
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
		Set<Entry<String, Gene>> set = chrAnnot.entrySet();
		TreeSet<Entry<String, Gene>> t = new TreeSet<Entry<String, Gene>>( new EntryNameComparator<Gene>() );
		t.addAll(set);
		Iterator<Entry<String, Gene>> it = t.iterator();
		Entry<String, Gene> e;
		ArrayList<Part> part = new ArrayList<Part>();
		int[] val = null;
		StringBuffer dnaSeqBuff = new StringBuffer();
		Part current;
		boolean[] used = new boolean[5000];
		boolean[] donS = new boolean[5000];
		boolean[] accS = new boolean[5000];
		String[] don = new String[5000];
		String[] acc = new String[5000];
		while( it.hasNext() ) {
			e = it.next();
			Gene gene = e.getValue();
			if( gene.transcript.size()>0 ) {
				int strand = gene.exon.get(0)[0];
				
				gene.reduce( e.getKey(), info );
				part.clear();
				int i, j;
				
				String[] id = new String[gene.transcript.size()];
				gene.transcript.keySet().toArray(id);
				Arrays.sort(id);
				Arrays.fill( accS, false );
				Arrays.fill( donS, false );
				for( int k = 0; k < id.length; k++ ) {
					IntList il = gene.transcript.get( id[k] );
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
					
					int off1 = val[1]-1>=intronic ? intronic : 0;
					int off2 = val[2]+intronic<=seq.length() ? intronic : 0;
					String p = seq.substring( val[1]-1-off1, val[2]+off2 );
					//check
					if( strand < 0 ) {
						p=Tools.rc(p);
					}
					String s;
					if( strand > 0 ) {
						s=p.substring(off1,p.length()-off2);
					} else {
						s=p.substring(off2,p.length()-off1);
					}
					acc[i]=don[i]="";
					if( s.matches( "[ACGT]*") ) {
						//acceptors
						if( accS[i] && ((strand>0 & off1>0) || (strand<0 && off2>0)) ) {
							out.get(4).writeln(">" + e.getKey() + "_" + i );
							//System.out.println( p.length() + "\t" + Arrays.toString(val) + "\t" + e.getKey() + "\t" + i + "\t" + p);
							out.get(4).writeln(p.substring(0, intronic+exonic));
							acc[i] = p.substring(intronic-2,intronic);
						}
						//donors
						if( donS[i] && ((strand>0 & off2>0) || (strand<0 && off1>0)) ) {
							out.get(5).writeln(">" + e.getKey() + "_" + i );
							out.get(5).writeln(p.substring(p.length()-(intronic+exonic)));
							don[i] = p.substring(p.length()-intronic,p.length()-intronic+2);
						}
					} else {
						s=null;
					}
					part.add(new Part(s));
				}
				
				Arrays.fill( used, false );
				for( int k = 0; k < id.length; k++ ) {
					int gt = 0, gc = 0, ag = 0;
					dnaSeqBuff.delete(0, dnaSeqBuff.length());
					
					String trans = id[k];
					IntList il = gene.transcript.get( trans );
					int start = strand>0 ? gene.exon.get(il.get(0))[1] : gene.exon.get(il.get(il.length()-1))[1];
					int end = strand>0 ? gene.exon.get(il.get(il.length()-1))[2] : gene.exon.get(il.get(0))[2];
					int offset = 0;
					for( j = 0; j < il.length(); j++ ) {
						current = part.get(il.get(j));
						if( j!= 0 && acc[il.get(j)].equalsIgnoreCase("AG") ) {
							ag++;
						}
						if( j+1 < il.length() ) {
							if( don[il.get(j)].equalsIgnoreCase("GT") ) {
								gt++;
							} else if( don[il.get(j)].equalsIgnoreCase("GC") ) {
								gc++;
							}
						}
						if( current.dna == null ) {
							break;
						}
						dnaSeqBuff.append( current.dna );
						
						//translate
						if( current.offsetLeft < 0 ) {
							current.offsetLeft = (3-offset)%3;
							//System.out.println( trans + "\t" + il.get(j) +"\t" +current.dna );
							current.aa = Tools.translate(current.offsetLeft, current.dna, code, false, Ambiguity.EXCEPTION);
							current.offsetRight = current.dna.length() - current.offsetLeft - 3*current.aa.length();
						} else {
							if( current.offsetLeft != (3-offset)%3 ) {
								break;
							}
						}
						offset = current.offsetRight;
					}
					
					if( j == il.length() ) {
						String p=null;
						
						p = Tools.translate(0, dnaSeqBuff.toString(), code, false, Ambiguity.EXCEPTION);
						int anz = 0, index=-1, last = 0;
						while( (index=p.indexOf('*',index+1))>= 0 ) {
							anz++;
							last = index;
						}
						if( anz > 1 ) {
							if( verbose ) protocol.appendWarning(trans + "\tskip premature stop, " + p + "\n" + dnaSeqBuff+"\n");
							problem[2]++;
						} else if( last != p.length()-1 ){
							if( verbose ) protocol.appendWarning(trans + "\tskip missing stop\n" );
							problem[1]++;
						} else if( p.charAt(0)!='M' ) {
							if( verbose ) protocol.appendWarning(trans + "\tskip missing start\n" );
							problem[3]++;
						} else {
							info[2]++;
							out.get(3).write( ">" + trans + "\n" + dnaSeqBuff.toString() + "\n" );
							out.get(2).write( ">" + trans + "\n" + p + "\n" );
							p = il.toString();
							out.get(1).write( e.getKey() + "\t" + trans + "\t" + p.substring(1,p.length()-1) + "\t" + chr + "\t" + strand + "\t" + start + "\t" + end + "\n" );
							for( j = 0; j < il.length(); j++ ) {
								used[il.get(j)] = true;
							}
																		
							int[] c = count.get(il.length());
							if( c == null ) {
								c = new int[1];
								count.put(il.length(), c);
							}
							c[0]++;
							
							//System.out.println(trans + "\t" + (j-1) + "\t"+gt+"\t"+gc+"\t"+ag);
						}
					} else {
						if( verbose ) protocol.appendWarning(trans + "\tskip ACGT coding part "+i +"\n");
						problem[0]++;
					}
				}
				
				for( j = 0; j < gene.exon.size(); j++ ) {
					if( used[j] && part.get(j).aa.length() > 0 ) {
						out.get(0).write(">" + e.getKey() + "_" + j + "\n" + part.get(j).aa + "\n");
					}
				}
			}
		}
	}
	
	private static class EntryNameComparator<T> implements Comparator<Entry<String,T>> {
		public int compare(Entry<String,T> o1, Entry<String,T> o2) {
			return o1.getKey().compareTo( o2.getKey() );
		}
	}
	
	static class Part {
		String dna, aa;
		int offsetLeft, offsetRight;
		
		Part( String dna ) {
			this.dna= dna;
			aa = null;
			offsetLeft = offsetRight = -100000;
		}
	}
	
	public ParameterSet getToolParameters() {
		try{
			return new SimpleParameterSet(
				new FileParameter( "annotation", "Reference annotation file (GFF), which contains gene models annotated in the reference genome", "gff", true ),
				new FileParameter( "genome", "Reference genome file (FASTA)", "fasta",  true ),

				new FileParameter( "genetic code", "optional user-specified genetic code", "tabular", false ),
					
				new SimpleParameter(DataType.BOOLEAN, "proteins", "whether the complete proteins sequences should returned as output", true, false ),
				new SimpleParameter(DataType.BOOLEAN, "transcripts", "whether the complete transcripts sequences should returned as output", true, false ),

				/*
				new SelectionParameter(DataType.PARAMETERSET, new String[]{"no","yes"}, new ParameterSet[]{
						new SimpleParameterSet(),
						new SimpleParameterSet(
								new SimpleParameter( DataType.INT, "intronic", "The number of bp return from the intron side", true, new NumberValidator<Integer>(0,1000), 10 ),//TODO
								new SimpleParameter( DataType.INT, "exonic", "The number of bp return from the exon side", true, new NumberValidator<Integer>(0,1000), 8 )
								//negative?
						)
					}, "splice sites", "whether splice sites should be returned or not", true ),
				*/	
				new FileParameter( "selected", "The path to list file, which allows to make only a predictions for the contained transcript ids. The first column should contain transcript IDs as given in the annotation. Remaining columns will be ignored.", "tabular,txt", maxSize>-1 ),
					
				new SimpleParameter( DataType.BOOLEAN, "verbose", "A flag which allows to output wealth of additional information", true, false )
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
		return "**What it does**\n\nThis tools can be used to create input files for GeMoMa, i.e., it creates at least a fasta file containing the translated parts of the CDS and a tabular file containing the assignment of transcripts to genes and parts of CDS to transcripts\n\n"
				+ "**Data format**\n\nThe input format is fasta for the genome and gff for the annotation. The gff file should be sorted.\n\nIf you like to set a user-specific genetic code, please use a tab-delimited file with two columns. The first contains the amino acid in one letter code. The second a list of tripletts. Here is a sample example\n\n"
					+ "+---+------------------------------+\n"
					+ "| I | ATT, ATC, ATA                |\n"
					+ "+---+------------------------------+\n"
					+ "| L | CTT, CTC, CTA, CTG, TTA, TTG |\n"
					+ "+---+------------------------------+\n"
					+ "| V | GTT, GTC, GTA, GTG           |\n"
					+ "+---+------------------------------+\n"
					+ "|...| ...                          |\n"
					+ "+---+------------------------------+\n\n"
				+ "**References**\n\nFor more information please contact jens.keilwagen@jki.bund.de.";
	}
	
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		ResultEntry[] re = new ResultEntry[2];
		for( int i = 0; i < re.length; i++ ) {
			re[i] = new ResultEntry(TextResult.class, type[i], name[i]);
		}
		return re;
	}

	@Override
	public String getToolVersion() {
		return "1.1.3";
	}
}