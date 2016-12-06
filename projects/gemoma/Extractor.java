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
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import de.jstacs.DataType;
import de.jstacs.parameters.EnumParameter;
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
import projects.gemoma.Tools.Ambiguity;

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

	private static String[] name = {"cds-parts", "assignment", "proteins", "cds"};
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
		int[] problem = new int[6];
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
		
		Ambiguity ambi = (Ambiguity) parameters.getParameterForName("Ambiguity").getValue();
		boolean fullLength = (Boolean) parameters.getParameterForName("full-length").getValue();
		boolean verbose = (Boolean) parameters.getParameterForName("verbose").getValue();
		
		ArrayList<File> file = new ArrayList<File>();
		ArrayList<SafeOutputStream> out = new ArrayList<SafeOutputStream>();
		getOut( name[0], file, out );
		getOut( name[1], file, out );
		getOut( ((Boolean)parameters.getParameterForName(name[2]).getValue()) ? name[2] : null, file, out );
		getOut( ((Boolean)parameters.getParameterForName(name[3]).getValue()) ? name[3] : null, file, out );
	
		out.get(1).writeln("#geneID\ttranscript\tcds-parts\tphases\tchr\tstrand\tstart\tend\tfull-length" );
		
		//read genome contig by contig
		r = new BufferedReader( new FileReader( parameters.getParameterForName("genome").getValue().toString() ) );
		
		StringBuffer seq = new StringBuffer();
		HashMap<String, int[]> donor = new HashMap<String, int[]>();
		HashMap<String, int[]> acceptor = new HashMap<String, int[]>();
		while( (line=r.readLine()) != null ) {
			if( line.startsWith(">") ) {
				//do
				extract( fullLength, ambi, protocol, verbose, comment, problem, info, count, seq, annot, code, out, donor, acceptor );
				//clear
				comment = line.substring(1);
				seq.delete(0, seq.length());
			} else {
				//add
				seq.append(line.trim().toUpperCase() );
			}
		}
		//do
		extract( fullLength, ambi, protocol, verbose, comment, problem, info, count, seq, annot, code, out, donor, acceptor );
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
		protocol.append( "missing start\t" + problem[3]+"\n");
		protocol.append( "no DNA\t" + problem[4]+"\n");
		protocol.append( "frame problems\t" + problem[5]+"\n\n");
		
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
	
	//gff has to be sorted 
	private static HashMap<String, HashMap<String,Gene>> readGFF( String input, HashMap<String,String> selected, Protocol protocol ) throws Exception {
		HashMap<String, HashMap<String,Gene>> annot = new HashMap<String, HashMap<String,Gene>>();
		HashMap<String,Gene> chr;
		Gene gene = null;
		BufferedReader r;
		String line, t;
		String[] split;
		int idx, h, end;
		
		//read transcripts
		r = new BufferedReader( new FileReader(input) );
		HashMap<String, Gene> trans = new HashMap<String, Gene>();
		ArrayList<String> cds = new ArrayList<String>();
		while( (line=r.readLine()) != null ) {
			if( line.equalsIgnoreCase("##FASTA") ) break; //http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
			if( line.length() == 0 || line.startsWith("#") ) continue; 
			idx = line.indexOf('\t')+1;
			idx = line.indexOf('\t',idx)+1;
			end = line.indexOf('\t',idx); 
			
			t = line.substring(idx,end);
			switch( t ) {
				case "CDS":
					cds.add(line);
					break;
				case "mRNA": case "transcript":
					split = line.split("\t");
					
					idx = split[8].indexOf("ID=")+3;
					h = split[8].indexOf(';',idx);
					String transcriptID = split[8].substring(idx, h>0?h:split[8].length() ).toUpperCase();
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
						
						chr = annot.get(split[0]);
						if( chr == null ) {
							chr = new HashMap<String,Gene>();
							annot.put(split[0], chr);
						}
						
						gene = chr.get(geneID);			
						if( gene == null ) {
							gene = new Gene(geneID,split[6]);
							chr.put(geneID, gene);
							
						}
						gene.add( transcriptID );
						trans.put(transcriptID, gene);
					}
					break;
			}
		}
		r.close();
		
		//read cds
		for( int i = 0 ; i < cds.size(); i++ ) {
			line = cds.get(i);
			split = line.split("\t");
				
			idx = split[8].indexOf(par)+par.length();
			h = split[8].indexOf(';',idx);
			String[] parent = split[8].substring(idx, h>0?h:split[8].length() ).toUpperCase().split(",");
			for( int j = 0; j < parent.length; j++ ) {
				gene = trans.get(parent[j]);
				if( gene != null && (selected==null || selected.containsKey(parent[j])) ) {
					gene.add( parent[j], new int[]{
							split[6].charAt(0)=='+'?1:-1, //strand
							Integer.parseInt( split[3] ), //start
							Integer.parseInt( split[4] ), //end
							split[7].charAt(0)=='.' ? -1 : Integer.parseInt(split[7]) //phase
					} );
				}
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
		
		Gene(String id, String strand) {
			transcript = new HashMap<String, IntList>();
			exon = new ArrayList<int[]>();
			this.id = id;
			this.strand = strand.charAt(0)=='+' ? 1: -1;
			start = end = -1;
		}
		
/*		Gene(String id, String start, String end, String strand) {
			this( id, strand );
			this.start = Integer.parseInt(start);
			this.end = Integer.parseInt(end);
		}*/
		
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
	
	private static void extract( boolean fullLength, Ambiguity ambi, Protocol protocol, boolean verbose, String comment, int[] problem, int[] info, HashMap<Integer,int[]> count,
			StringBuffer seq, HashMap<String, HashMap<String,Gene>> annot, HashMap<String,Character> code,
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
						s=null;//TODO
					}					
					part.add(new Part(s,val[3]));
				}
				
				Arrays.fill( used, false );
				int currentProb=-1;
				for( int k = 0; k < id.length; k++ ) {
					dnaSeqBuff.delete(0, dnaSeqBuff.length());
					
					String trans = id[k];
					IntList il = gene.transcript.get( trans );
					if( il.length() == 0 ) {
						System.out.println("No coding exon(s) for: " + id[k] );
						continue;
					}
					int start = strand>0 ? gene.exon.get(il.get(0))[1] : gene.exon.get(il.get(il.length()-1))[1];
					int end = strand>0 ? gene.exon.get(il.get(il.length()-1))[2] : gene.exon.get(il.get(0))[2];
					
					int startPhase = fullLength ? 0 : part.get(il.get(0)).offsetLeft;
					int offset = 3-startPhase, pa = -1;
					for( j = 0; j < il.length(); j++ ) {
						pa = il.get(j);
						current = part.get(pa);
						if( current.dna == null ) {
							currentProb=1;
							break;
						}
						dnaSeqBuff.append( current.dna );
						
						//translate
						if( current.offsetLeft < 0 ) {
							current.offsetLeft = (3-offset)%3;
						}
						
						if( current.aa == null ) {
							//System.out.println( trans + "\t" + pa +"\t" +current.dna );
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
						offset = current.offsetRight;
					}

					String p=null;
					if( j == il.length() ) {
						if( ambi == Ambiguity.EXCEPTION && !dnaSeqBuff.toString().matches("[ACGT]*") ) {//to be compatible with older versions
							j=il.length()+1;
							currentProb=0;
						} else {
							try {
								p = Tools.translate(startPhase, dnaSeqBuff.toString(), code, false, ambi);
							} catch( IllegalArgumentException iae ) {
								j=il.length()+1;
								currentProb=0;
							}
						}
					}
					
					if( j == il.length() ) {
						int anz = 0, index=-1, last = 0;
						while( (index=p.indexOf('*',index+1))>= 0 ) {
							anz++;
							last = index;
						}
						if( anz > 1 ) {
							if( verbose ) protocol.appendWarning(trans + "\tskip premature stop, " + p + "\n" + dnaSeqBuff+"\n");
							problem[2]++;
						} else if( fullLength && last != p.length()-1 ){
							if( verbose ) protocol.appendWarning(trans + "\tskip missing stop\n" );
							problem[1]++;
						} else if( fullLength && p.charAt(0)!='M' ) {
							if( verbose ) protocol.appendWarning(trans + "\tskip missing start\n" );
							problem[3]++;
						} else {
							info[2]++;
							out.get(3).write( ">" + trans + "\n" + dnaSeqBuff.toString() + "\n" );
							out.get(2).write( ">" + trans + "\n" + p + "\n" );
							String x = il.toString();
							SafeOutputStream sos = out.get(1);
							sos.write( gene.id + "\t" + trans + "\t" + x.substring(1,x.length()-1) );
							for( j = 0; j < il.length(); j++ ) {
								sos.write( (j==0?"\t":",") + part.get(il.get(j)).offsetLeft );		
							}
							sos.write( "\t" + chr + "\t" + strand + "\t" + start + "\t" + end + "\t" + (p.charAt(0)=='M' && p.charAt(p.length()-1)=='*') + "\n" );
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
							last=-1;
							for( j = 0; j < il.length(); j++ ) {
								pa = il.get(j);
								current = part.get(pa);
								if( last != -1 ) {
									splits[last][pa] = true;
								}
								last = pa;
								
								//find splice sites
								int[] exon = gene.exon.get(pa);
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
								
								if( current.dna.length()>0 ) {
									if( j > 0 && acc[pa].length() > 0 ) {
										int[] stat = acceptor.get(acc[pa]);
										if( stat == null ) {
											stat = new int[1];
											acceptor.put(acc[pa], stat);
										}
										stat[0]++;
									}
									if( j+1<il.length() && don[pa].length() > 0 ) {
										int[] stat = donor.get(don[pa]);
										if( stat == null ) {
											stat = new int[1];
											donor.put(don[pa], stat);
										}
										stat[0]++;
									}
								}
							}
						}
					} else {
						switch (currentProb ) {
							case 0:
								if( verbose ) {
									if( j < il.length() ) {
										protocol.appendWarning(trans + "\tskip non-ACGT coding part "+j+"\n");
									} else {
										protocol.appendWarning(trans + "\tskip non-ACGT coding protein\n");
									}
								}
								problem[0]++;
								break;
							case 1:
								if( verbose ) protocol.appendWarning(trans + "\tskip no DNA for coding part "+j+"\n");
								problem[4]++;
								break;
							case 2:
								if( verbose ) protocol.appendWarning(trans + "\tskip frame problems for coding part "+j+"\n");
								problem[5]++;
								break;
						}
						
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
	
	static class Part {
		String dna, aa;
		int offsetLeft, offsetRight;
		
		Part( String dna ) {
			this(dna, -100000);
		}
		
		Part( String dna, int phase ) {
			this.dna= dna;
			aa = null;
			offsetLeft = phase;
			offsetRight = -100000;
		}
	}
	
	public ParameterSet getToolParameters() {
		try{
			return new SimpleParameterSet(
				new FileParameter( "annotation", "Reference annotation file (GFF), which contains gene models annotated in the reference genome", "gff", true ),
				new FileParameter( "genome", "Reference genome file (FASTA)", "fasta",  true ),

				new FileParameter( "genetic code", "optional user-specified genetic code", "tabular", false ),
					
				new SimpleParameter(DataType.BOOLEAN, Extractor.name[2], "whether the complete proteins sequences should returned as output", true, false ),
				new SimpleParameter(DataType.BOOLEAN, Extractor.name[3], "whether the complete CDSs should returned as output", true, false ),
				
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
				new EnumParameter( Ambiguity.class, "This parameter defines how to deal with ambiguities in the DNA. There are 3 options: "
						+ "EXCEPTION, which will remove the corresponding transcript, "
						+ "AMBIGUOUS, which will use an X for the corresponding amino acid, and "
						+ "RANDOM, which will randomly select an amnio acid from the list of possibilities.", true, Ambiguity.EXCEPTION.toString() ),	
				new SimpleParameter( DataType.BOOLEAN, "full-length", "A flag which allows for choosing between only full-length and all (i.e., full-length and partial) transcripts", true, true ),
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
				+ "**References**\n\nFor more information please visit http://www.jstacs.de/index.php/GeMoMa or contact jens.keilwagen@julius-kuehn.de.\n";
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
		return "1.3";
	}
}