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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import projects.gemoma.Tools.Ambiguity;
import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
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
	
	private static BufferedWriter intron;
	
	@Override
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		progress.setIndeterminate();
		intron = null;//new BufferedWriter(new FileWriter("intron.gff"));//TODO
		
		HashMap<String, String> selected = null;
		BufferedReader r;
		String line, comment=null;
		Parameter p = parameters.getParameterForName("selected"); 
		if( p.isSet() ) {
			selected = Tools.getSelection( p.getValue().toString(), maxSize, protocol );
			protocol.append("selected: " + selected.size() + "\t"+ selected+"\n");
		}
		
		HashMap<String, HashMap<String,Gene>> annot = readGFF( parameters.getParameterForName("annotation").getValue().toString(), selected, protocol );
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
		exonic = intronic = 2;
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
		HashMap<String, int[]> donor = new HashMap<String, int[]>();
		HashMap<String, int[]> acceptor = new HashMap<String, int[]>();
		while( (line=r.readLine()) != null ) {
			if( line.startsWith(">") ) {
				//do
				extract( protocol, verbose, comment, problem, info, count, intronic, exonic, seq, annot, code, out, donor, acceptor );
				//clear
				comment = line.substring(1);
				seq.delete(0, seq.length());
			} else {
				//add
				seq.append(line.trim().toUpperCase() );
			}
		}
		//do
		extract( protocol, verbose, comment, problem, info, count, intronic, exonic, seq, annot, code, out, donor, acceptor );
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
		
		protocol.append( "\ngenes\t" + info[0] +"\n");
		protocol.append( "identical CDS of same gene\t" + info[1] +"\n");
		protocol.append( "transcripts\t" + info[2]+"\n\n");
		
		protocol.append( "reasons for discarding transcripts:\n");
		protocol.append( "ambigious nucleotide\t" + problem[0] +"\n");
		protocol.append( "missing stop\t" + problem[1] +"\n");
		protocol.append( "premature stop\t" + problem[2]+"\n");
		protocol.append( "missing start\t" + problem[3]+"\n\n");
		
		//log
		protocol.append("coding exons\t#\n");
		Integer[] array = new Integer[count.size()];
		count.keySet().toArray(array);
		Arrays.sort(array);
		for( int i = 0; i < array.length; i++ ) {
			protocol.append( array[i] + "\t" + count.get(array[i])[0] + "\n");
		}
		
		Iterator<Entry<String,int[]>> it;
		Entry<String,int[]> e;
		if( acceptor.size() > 0 ) {
			protocol.append("\nacceptor\t#\n");
			it = acceptor.entrySet().iterator();
			while( it.hasNext() ) {
				e = it.next();
				protocol.append( e.getKey() + "\t" + e.getValue()[0] + "\n");
			}
		}
		if( donor.size() > 0 ) {
			protocol.append("\ndonor\t#\n");
			it = donor.entrySet().iterator();
			while( it.hasNext() ) {
				e = it.next();
				protocol.append( e.getKey() + "\t" + e.getValue()[0] + "\n");
			}
		}
		
		if( intron != null ) {
			intron.close();
		}
		return new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
	}

	private static String par = "Parent=";
	private static String geneTag = "gene";
	
	//gff has to be sorted 
	private static HashMap<String, HashMap<String,Gene>> readGFF( String input, HashMap<String,String> selected, Protocol protocol ) throws Exception {
		HashMap<String, HashMap<String,Gene>> annot = new HashMap<String, HashMap<String,Gene>>();
		HashMap<String,Gene> chr;
		Gene gene = null;
		BufferedReader r;
		String line, geneID = null, transcriptID, t;
		String[] split;
		int idx, h, end;
		
		//read genes
		r = new BufferedReader( new FileReader(input) );
		ArrayList<String> transcript = new ArrayList<String>();
		ArrayList<String> cds = new ArrayList<String>();
		while( (line=r.readLine()) != null ) {
			if( line.equalsIgnoreCase("##FASTA") ) break; //http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
			if( line.length() == 0 || line.startsWith("#") ) continue; 
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
				chr.put(geneID, new Gene(geneID,split[3],split[4],split[6]));
			} else if( t.equalsIgnoreCase( "mRNA" ) || t.equalsIgnoreCase("transcript") ) {
				transcript.add(line);
			}
		} while( (line=r.readLine()) != null );		
		r.close();
		
		//read transcripts
		HashMap<String, Gene> trans = new HashMap<String, Gene>();
		String[] parent;
		for( int i = 0 ; i < transcript.size(); i++ ) {
			line = transcript.get(i);
			split = line.split("\t");
			
			idx = split[8].indexOf("ID=")+3;
			h = split[8].indexOf(';',idx);
			transcriptID = split[8].substring(idx, h>0?h:split[8].length() ).toUpperCase();
			if( selected == null || selected.containsKey(transcriptID) ) {
				idx = split[8].indexOf(par)+par.length();
				h = split[8].indexOf(';',idx);
				parent = split[8].substring(idx, h>0?h:split[8].length() ).split(",");
				HashMap<String,Gene> x = annot.get(split[0]);
				int j = 0;
				gene = null;
				if( x != null) { 
					while( j < parent.length && (gene = x.get(parent[j])) == null ) {
						j++;
					}
				} else {
					protocol.appendWarning("Could not parse a feature \""+geneTag+"\" on sequence \"" + split[0] + "\": " + line + ".\n" );
				}
				if( gene != null ) {
					gene.add( transcriptID );
					trans.put(transcriptID, gene);
				}
			}
		}
		
		//read cds
		for( int i = 0 ; i < cds.size(); i++ ) {
			line = cds.get(i);
			split = line.split("\t");
				
			idx = split[8].indexOf(par)+par.length();
			h = split[8].indexOf(';',idx);
			parent = split[8].substring(idx, h>0?h:split[8].length() ).toUpperCase().split(",");
			int j = 0;
			while( j < parent.length && (gene = trans.get(parent[j]) ) == null ) {
				j++;
			}
			if( gene != null && (selected==null || selected.containsKey(parent[j])) ) {
				gene.add( parent[j], new int[]{
						split[6].charAt(0)=='+'?1:-1, //strand
						Integer.parseInt( split[3] ), //start
						Integer.parseInt( split[4] ) //end
				} );
			}
		}
		return annot;
	}
	
	private static class Gene implements Comparable<Gene>{
		HashMap<String,IntList> transcript;
		ArrayList<int[]> exon;
		int start, end;
		int strand;
		String id;
		
		Gene(String id, String start, String end, String strand) {
			transcript = new HashMap<String, IntList>();
			exon = new ArrayList<int[]>();
			this.id = id;
			this.start = Integer.parseInt(start);
			this.end = Integer.parseInt(end);
			this.strand = strand.charAt(0)=='+' ? 1: -1;
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

		@Override
		public int compareTo(Gene o) {
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
	
	private static void extract( Protocol protocol, boolean verbose, String comment, int[] problem, int[] info, HashMap<Integer,int[]> count,
			int intronic, int exonic, StringBuffer seq, HashMap<String, HashMap<String,Gene>> annot, HashMap<String,Character> code,
			ArrayList<SafeOutputStream> out, HashMap<String, int[]> donor, HashMap<String, int[]> acceptor) throws Exception {
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

		ArrayList<Part> part = new ArrayList<Part>();
		int[] val = null;
		StringBuffer dnaSeqBuff = new StringBuffer();
		Part current;
		boolean[] used = new boolean[5000];
		boolean[] donS = new boolean[5000];
		boolean[] accS = new boolean[5000];
		String[] don = new String[5000];
		String[] acc = new String[5000];
		StringBuffer aa = new StringBuffer();
		StringBuffer spliceSeq = new StringBuffer();
		for( Gene gene: genes ) {
			if( gene.transcript.size()>0 ) {
				boolean[][] splits = new boolean[gene.exon.size()][gene.exon.size()];
				for( int k = 0; k < splits.length; k++ ) {
					Arrays.fill( splits[k], false );
				}
				
				int strand = gene.strand;
				boolean forward=strand==1;
				
				gene.reduce( gene.id, info );
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
					if( !s.matches( "[ACGT]*") ) {
						s=null;
					}
					part.add(new Part(s));
				}
				
				Arrays.fill( used, false );
				for( int k = 0; k < id.length; k++ ) {
					//int gt = 0, gc = 0, ag = 0;
					dnaSeqBuff.delete(0, dnaSeqBuff.length());
					
					String trans = id[k];
					IntList il = gene.transcript.get( trans );
					if( il.length() == 0 ) {
						System.out.println("No coding exon(s) for: " + id[k] );
						continue;
					}
					int start = strand>0 ? gene.exon.get(il.get(0))[1] : gene.exon.get(il.get(il.length()-1))[1];
					int end = strand>0 ? gene.exon.get(il.get(il.length()-1))[2] : gene.exon.get(il.get(0))[2];
					int offset = 0;
					for( j = 0; j < il.length(); j++ ) {
						current = part.get(il.get(j));
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
							SafeOutputStream sos = out.get(1);
							sos.write( gene.id + "\t" + trans + "\t" + p.substring(1,p.length()-1) );
							for( j = 0; j < il.length(); j++ ) {
								sos.write( (j==0?"\t":",") + part.get(il.get(j)).offsetLeft );		
							}
							sos.write( "\t" + chr + "\t" + strand + "\t" + start + "\t" + end + "\n" );
							for( j = 0; j < il.length(); j++ ) {
								used[il.get(j)] = true;
							}
																		
							int[] c = count.get(il.length());
							if( c == null ) {
								c = new int[1];
								count.put(il.length(), c);
							}
							c[0]++;
							
							//splice sites
							int MAX_INTRON_LENGTH=15000, ignoreAAForSpliceSite=30, targetStart, targetEnd;
							int gt = 0, gc = 0, ag = 0;
							last=-1;
							for( j = 0; j < il.length(); j++ ) {
								current = part.get(il.get(j));
								if( last != -1 ) {
									splits[last][il.get(j)] = true;
								}
								last = il.get(j);
								
								//find splice sites
								int[] exon = gene.exon.get(il.get(j));
								if( forward ) {
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
								
//System.out.println(j + "\t" + Arrays.toString(exon) + "\t" + targetStart + "\t" + targetEnd + "\t" + current.offsetLeft + "\t" + current.offsetRight + "\t" + add);
								if( current.dna.length()>0 ) {

									if( j!=0 && current.acc==null ) {
										aa.delete(0, aa.length());
										spliceSeq.delete(0, spliceSeq.length());
										Tools.getMaximalExtension(seq, forward, false, -strand, forward ? targetStart+add-1 : (targetEnd-add), '*','*', spliceSeq, aa, current.aa.substring(0,add/3), exonic, intronic-1, MAX_INTRON_LENGTH, code );
										int d = aa.length()-add/3-(current.offsetLeft>0?1:0);
/*
System.out.print( getFiller(intronic-1, '.' ) );
System.out.print( matchAA2DNA(aa.subSequence(0, d),'>') );
System.out.println( getFiller( (3-current.offsetLeft)%3, '.') );
System.out.println( spliceSeq );
*/
										int border = intronic-1 + 3*d + ((3-current.offsetLeft)%3);
										current.acc = spliceSeq.substring(0, border).toLowerCase() + spliceSeq.substring(border).toUpperCase();
										
										out.get(4).writeln(">" + gene.id + "_" + il.get(j) + "\t" + border );
										out.get(4).writeln( current.acc );
										acc[il.get(j)] = spliceSeq.substring(border-2, border);
										
										int[] stat = acceptor.get(acc[il.get(j)]);
										if( stat == null ) {
											stat = new int[1];
											acceptor.put(acc[il.get(j)], stat);
										}
										stat[0]++;
									}
									if( j+1<il.length() && current.don==null ){
										aa.delete(0, aa.length());
										spliceSeq.delete(0, spliceSeq.length());
										
										Tools.getMaximalExtension(seq, forward, true, forward?1:-1, forward ? targetEnd-add : (targetStart+add-1), '*','*', spliceSeq, aa, current.aa.substring(current.aa.length() - add/3), exonic, intronic-1, MAX_INTRON_LENGTH, code );
/*								
System.out.print( getFiller(exonic, '.' ) );
System.out.print( matchAA2DNA(current.aa.substring(current.aa.length() - add/3), '>') );
System.out.println( getFiller( current.offsetRight, '.') );
System.out.println(spliceSeq);
*/	
										int border = exonic + add + current.offsetRight;
										current.don = spliceSeq.substring(0, border).toUpperCase() + spliceSeq.substring(border).toLowerCase();
										
										out.get(5).writeln(">" + gene.id + "_" + il.get(j) + "\t" + border );
										out.get(5).writeln( current.don );
										don[il.get(j)] = spliceSeq.substring(border, border+2);
										
										int[] stat = donor.get(don[il.get(j)]);
										if( stat == null ) {
											stat = new int[1];
											donor.put(don[il.get(j)], stat);
										}
										stat[0]++;
									}
									
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
								}
							}
							//TODO? System.out.println(trans + "\t" + (j-1) + "\t"+gt+"\t"+gc+"\t"+ag);
						}
					} else {
						if( verbose ) protocol.appendWarning(trans + "\tskip ACGT coding part "+i +"\n");
						problem[0]++;
					}
				}
				
				for( j = 0; j < gene.exon.size(); j++ ) {
					if( used[j]) {
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
						
						if( intron != null ) {
							for( int k = 0; k < splits.length; k++ ) {
								if( splits[j][k] ) {
									String intr;
									int st, en;
									if( forward ) {
										st = gene.exon.get(j)[2]+1;
										en = gene.exon.get(k)[1];
									} else {
										st = gene.exon.get(k)[2]+1;
										en = gene.exon.get(j)[1];
									}
									intr = seq.substring(st-1, en-1);
									if( !forward ) {
										intr = Tools.rc(intr);
									}
									//System.out.println( Arrays.toString( gene.exon.get(j) ) );
									//System.out.println( Arrays.toString( gene.exon.get(k) ) );
									//System.out.println(forward + "\t" + intr);
									intron.append(chr + "\tannotation\tintron\t" + st + "\t" + en + "\t.\t" + (forward?"+":"-") + "\t.\t." );
									intron.newLine();
								}
							}
						}
					}
				}
			}
		}
	}
	
	private static String getFiller( int anz, char c ) {
		char[] ch = new char[anz];
		Arrays.fill(ch, c);
		return new String(ch);
	}
	
	private static String matchAA2DNA( CharSequence up, char c ) {
		char[] ch = new char[3*up.length()];
		for( int i = 0; i < up.length(); i++ ) {
			ch[3*i] = up.charAt(i);
			ch[3*i+1] = ch[3*i+2] = c;
		}
		return new String(ch);
	}
	
	private static class EntryNameComparator<T> implements Comparator<Entry<String,T>> {
		public int compare(Entry<String,T> o1, Entry<String,T> o2) {
			return o1.getKey().compareTo( o2.getKey() );
		}
	}
	
	static class Part {
		String dna, aa, don, acc;
		int offsetLeft, offsetRight;
		
		Part( String dna ) {
			this.dna= dna;
			aa = null;
			don=null;
			acc=null;
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
								new SimpleParameter( DataType.INT, "intronic", "The number of bp return from the intron side", true, new NumberValidator<Integer>(2,1000), 10 ),
								new SimpleParameter( DataType.INT, "exonic", "The number of bp return from the exon side", true, new NumberValidator<Integer>(0,1000), 8 )
								//negative?
						)
					}, "splice sites", "whether splice sites should be returned or not", true ),
				/**/	
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