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
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Random;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

import javax.script.Bindings;
import javax.script.ScriptEngine;
import javax.script.ScriptException;

import de.jstacs.parameters.Parameter;
import de.jstacs.tools.Protocol;

/**
 * Some methods used for GeMoMa.
 * 
 * @author Jens Keilwagen
 */
public class Tools {
	
	public static final String GeMoMa_TEMP = "GeMoMa_temp" + File.separator;//TODO
	
	/**
	 * 
	 * @param parameter the parameter that defines the file for the {@link InputStream}
	 * @param alternative the alternative that is used if the parameter is not set
	 * 
	 * @return an {@link InputStream}
	 * 
	 * @throws FileNotFoundException the {@link File} cannot be found
	 * 
	 * @see {@link Parameter#isSet()}
	 */
	public static InputStream getInputStream( Parameter parameter, String alternative ) throws FileNotFoundException {
		InputStream in;
		if( parameter != null && parameter.isSet() ) {
			in = new FileInputStream( parameter.getValue().toString() );
		} else {
			in = Tools.class.getClassLoader().getResourceAsStream( alternative );
		}
		return in;
	}
	
	/**
	 * Create a temporary file with prefix &quot;GeMoMa&quot; and user-specified infix.
	 * The method {@link File#deleteOnExit()} is invoked.  
	 * 
	 * @param infix an infix for the file
	 * 
	 * @return a temporary file
	 * 
	 * @throws IOException
	 * 
	 * @see {@link #GeMoMa_TEMP}
	 * @see {@link File#deleteOnExit()}
	 */
	public static File createTempFile( String infix ) throws IOException {
			return createTempFile(infix, GeMoMa_TEMP);
	}
	
	public static File createTempFile( String infix, String temp ) throws IOException {
		File d = new File(temp);
		if( !d.exists() ) {
			d.mkdirs();
		}
		File f = File.createTempFile("GeMoMa-"+infix + "-", ".temp", d);
		f.deleteOnExit();
		return f;
	}
	
	public static File[] externalSort( String input, int size, int files, Protocol protocol, boolean cdsParts ) throws Exception {
		return externalSort(input, size, files, protocol, cdsParts, GeMoMa_TEMP);
	}
	
	public static File[] externalSort( String input, int size, int files, Protocol protocol, boolean cdsParts, String temp ) throws Exception {
		protocol.append("sorting the search results\n");
		BufferedReader r = new BufferedReader( new FileReader(input) );
		ArrayList<File> f = new ArrayList<File>();
		String[] line = new String[size];
		int anz = 0, a = 0, shortCut = 0;
		//read, split and sort parts;
		do {
			a=0;
			while( a < size && (line[a]=r.readLine()) != null ) {
				a++;
			}
			Arrays.sort(line,0,a);
			f.add( createTempFile("sort-"+anz, temp) );
			BufferedWriter w = new BufferedWriter( new FileWriter(f.get(anz)) );
			for( int i = 0; i < a; i++ ) {
				w.append(line[i]);
				w.newLine();
			}
			w.close();
			anz++;
		} while( a == size );
		r.close();
		protocol.append("files:\t" + anz + "\n" );
		//System.out.println("read\t" + t.getElapsedTime() );
		
		//merge
		line = new String[anz];
		BufferedReader[] parts = new BufferedReader[anz];
		for( int i = 0; i < anz; i++ ) {
			parts[i]  = new BufferedReader( new FileReader(f.get(i)) );
			line[i] = parts[i].readLine();
		}
		File[] sorted = new File[files];
		BufferedWriter[] w = new BufferedWriter[files];
		for( int s = 0; s < files; s++ ) {
			sorted[s] = createTempFile("sorted",temp);
			w[s] = new BufferedWriter( new FileWriter( sorted[s] ) );
		}
		int s = -1;
		String old = null;
		int flip = 0;
		while( true ) {
			int best = -1;
			for( int i = 0; i < anz; i++ ) {
				if( line[i] != null ) {
					if( best == -1 ) {
						best=i;
					} else {
						if( line[i].compareTo(line[best]) < 0 ) {
							best = i;
						}
					}
				}
			}
			if( best >= 0 ) {
				/*
				w.append(line[best]);
				w.newLine();
				line[best] = parts[best].readLine();
				/**/
				
				String start = line[best].substring(0, line[best].indexOf('\t')+1);
				if( old == null || !start.startsWith(old) ) {
					if( cdsParts ) {
						int idx = start.lastIndexOf('_');
						if( idx < 0 ) throw new IllegalArgumentException("You selected cdsParts=true, but the ID ("+start+")seems to be no CDS part.");
						old = start.substring(0,idx+1);
					} else {
						old = start;
					}
					s++;
					if( s == sorted.length ) {
						flip++;
						s=0;
					}
				}
				int b = 0;
				do {
					w[s].append(line[best]);
					w[s].newLine();
					b++;
				} while( (line[best] = parts[best].readLine()) != null && line[best].startsWith(start) );
				shortCut += (b-1);/**/ 
			} else {
				break;
			}
		}
		for( int i = 0; i < anz; i++ ) {
			parts[i].close();
		}
		for( int i = 0; i < files; i++ ) {
			w[i].close();
		}
		if( flip<1 ) {
			File[] small = new File[s+1];
			System.arraycopy(sorted, 0, small, 0, small.length);
			sorted=small;
		}
		protocol.append( "shortCuts:\t" + shortCut + "\n" );
		return sorted;
	}
	
	/**
	 * This method determines the maximal extension of a hit until the amino acids <code>c1</code> or <code>c2</code> occur (typically START (M), STOP (*))
	 * 
	 * @param chr the chromosome/contig
	 * @param end a switch for deciding where to insert
	 * @param forward the strand of the hit
	 * @param direction the direction of the search
	 * @param startPos the position for starting
	 * @param c1 typically START (M), STOP (*)
	 * @param c2 typically START (M), STOP (*)
	 * @param seq contains the DNA sequence in correct orientation
	 * @param translated contains the amino acid sequence in correct orientation
	 * @param backup backup sequence for ambiguous codons in the alignment
	 */
	public static void getMaximalExtension( CharSequence chr, boolean forward, boolean end, int direction, int startPos, char c1, char c2, StringBuffer seq, StringBuffer translated, String backup, int upstream, int downstream, int MAX_INTRON_LENGTH,  HashMap<String, Character> code ) {
		int s = startPos, e=-1;
		int i=0;
		CharSequence add;
		if( translated != null ) {
			translated.delete(0, translated.length() );
			seq.delete(0, seq.length() );
			
			//upstream of startPos
			if( direction==1 ) {
				s=Math.max(startPos-upstream,0);
				add = chr.subSequence(s,startPos);
			} else {
				s=Math.min(startPos+upstream,chr.length());
				add = chr.subSequence(startPos,s);
			}
			if( !forward ) {
				add = Tools.rc(add);
			}
			insert( end, seq, add );
			int l = add.length();
			while( l < upstream ) {
				insert( end, seq, "N" );
				l++;
			}
		}
		char c = '!';
		CharSequence codon = null;
		
		if( direction == -1 ) {
			startPos -= 3;
		}
		
		int problem = 0;
		boolean onlyStop = c1==c2 && c1=='*';
		boolean fromBackup = true;
		while( 3*i < MAX_INTRON_LENGTH //TODO  for runtime reasons (in N-stretches) 
				&& startPos >= 0 && startPos+3 <= chr.length() && ((fromBackup && onlyStop) || (c != c1 && c != c2)) ) {
			codon = chr.subSequence(startPos,startPos+3);
			if( !forward ) {
				codon = Tools.rc(codon);
			}
			fromBackup = backup != null && i < backup.length();
			if( fromBackup ) {
				c = backup.charAt( end ? i : (backup.length()-1-i) );
				problem = 0;
			} else {
				try{ 
					c = Tools.translate(codon, code, Ambiguity.EXCEPTION);
					problem = 0;
				} catch( Exception ex ) {
					problem++;
					if( problem == 3 ) {
						break;
					} else {
						c='X';
					}
				}
			}
			if( seq != null ) {
				insert( end, seq, codon );
				insert( end, translated, ""+c );
			}
			startPos=startPos+direction*3;
			i++;
		}
		
		if( direction == -1 ) {
			e = startPos+3;
		} else {
			e = startPos;
		}
		startPos=e;
		
		if( seq != null ) {
			/*old
			if( problem > 0 ) {
				if( startPos >= 0 && startPos+3 <= chr.length() ) {
					codon = chr.substring(startPos,startPos+3);
					if( !forward ) {
						codon = Tools.rc(codon);
					}
				} else {
					codon="NNN";
				}
			}
			insert( end, seq, codon != null? codon.substring(end?0:intronic, 1+(end?0:intronic)) : "N" );
			*/
			
			//downstream of the region
			if( direction==1 ) {
				//System.out.println( startPos + "\t" + Math.min(startPos+intronic-1,chr.length()) );
				e=Math.min(startPos+downstream,chr.length());
				add = chr.subSequence( startPos, e );
			} else {
				e=Math.max(startPos-downstream,0);
				add = chr.subSequence( e, startPos );
			}
			if( !forward ) {
				add = Tools.rc(add);
			}

			insert( end, seq, add );
			int l = add.length();
			while( l < downstream ) {
				insert( end, seq, "N" );
				l++;
			}
		}

/*
//test of correctness
		CharSequence h = direction==1 ? chr.substring(s, e) : chr.subSequence(e, s);
		if( !forward ){
			h = Tools.rc(h);
		}
		System.out.println(h);
		System.out.println(seq);
/**/
	}
	
	/**
	 * This method allows for inserting one amino acid or codon at the beginning or the end of a {@link StringBuffer}.
	 * 
	 * @param end a switch for deciding where to insert the {@link String}
	 * @param sb the {@link StringBuffer} representing the sequence
	 * @param s to be inserted
	 */
	static void insert( boolean end, StringBuffer sb, CharSequence s ) {
		if( end ) {
			sb.append( s );
		} else {
			sb.insert(0, s);
		}
	}
	
	
	public static HashMap<String,String> getSelection( String fName, int maxSize, Protocol protocol ) throws IOException {
		if( fName == null ) {
			return null;			
		}
		HashMap<String,String> selected = new HashMap<String, String>();
		BufferedReader r = new BufferedReader( new FileReader( fName ) );
		String line;
		while( (line=r.readLine()) != null && (maxSize<0 || selected.size() < maxSize) ) {
			int idx = line.indexOf('\t'), second = idx+1;
			if( idx < 0 ) {
				second = idx = line.length();
			}
			selected.put(line.substring(0,idx), line.substring(second));
		}
		if( maxSize >= 0 ) {
			protocol.appendWarning("Only used the first " + maxSize + " lines of selected gene IDs.");
		}
		r.close();
		return selected;
	}
	
	public static HashMap<String,Character> getCode( String fName ) throws Exception {
		return getCode( new FileInputStream(fName) );
	}
	
	public static HashMap<String,Character> getCode( InputStream input ) throws Exception {
		HashMap<String, Character> annot = new HashMap<String, Character>();
		BufferedReader r = new BufferedReader( new InputStreamReader(input) );
		String line;
		while( (line=r.readLine()) != null ) {
			char code = line.charAt(0);
			line = line.substring(2);
			String[] triplets = line.split(",");
			for( String t: triplets ) {
				t= t.trim();
				Character c = annot.get(t);
				if( c == null ) {
					annot.put(t, code);
				} else {
					throw new IllegalArgumentException( t + ": " + line );
				}
			}
		}
		r.close();
		return annot;
	}
	
	
	public static char[] rc( char[] in ) {
		char[] out = new char[in.length]; 
		for( int i = 0; i < in.length; i++ ) {
			switch( in[i] ) {
				case 'a': case 'A': out[in.length-1-i]='T'; break;
				case 'c': case 'C': out[in.length-1-i]='G'; break;
				case 'g': case 'G': out[in.length-1-i]='C'; break;
				case 't': case 'T': out[in.length-1-i]='A'; break;
				
				case 'n': case 'N': out[in.length-1-i]='N'; break;
				
				case 'y': case 'Y': out[in.length-1-i]='R'; break;
				case 'r': case 'R': out[in.length-1-i]='Y'; break;
				
				case 'w': case 'W': out[in.length-1-i]='W'; break;
				case 's': case 'S': out[in.length-1-i]='S'; break;
				
				case 'k': case 'K': out[in.length-1-i]='M'; break;
				case 'm': case 'M': out[in.length-1-i]='K'; break;
				
				case 'v': case 'V': out[in.length-1-i]='T'; break;
				case 'h': case 'H': out[in.length-1-i]='G'; break;
				case 'd': case 'D': out[in.length-1-i]='C'; break;
				case 'b': case 'B': out[in.length-1-i]='A'; break;

				case '-': out[in.length-1-i]='-'; break;

				
				default: throw new IllegalArgumentException("unknown character: '"+in[i]+"'");
			}
		}
		return out;
	}
	
	public static String rc( CharSequence in ) {
		char[] out = new char[in.length()]; 
		for( int i = 0; i < in.length(); i++ ) {
			switch( in.charAt(i) ) {
				case 'a': case 'A': out[in.length()-1-i]='T'; break;
				case 'c': case 'C': out[in.length()-1-i]='G'; break;
				case 'g': case 'G': out[in.length()-1-i]='C'; break;
				case 't': case 'T': out[in.length()-1-i]='A'; break;
				
				case 'n': case 'N': out[in.length()-1-i]='N'; break;
				
				case 'y': case 'Y': out[in.length()-1-i]='R'; break;
				case 'r': case 'R': out[in.length()-1-i]='Y'; break;
				
				case 'w': case 'W': out[in.length()-1-i]='W'; break;
				case 's': case 'S': out[in.length()-1-i]='S'; break;
				
				case 'k': case 'K': out[in.length()-1-i]='M'; break;
				case 'm': case 'M': out[in.length()-1-i]='K'; break;
				
				case 'v': case 'V': out[in.length()-1-i]='T'; break;
				case 'h': case 'H': out[in.length()-1-i]='G'; break;
				case 'd': case 'D': out[in.length()-1-i]='C'; break;
				case 'b': case 'B': out[in.length()-1-i]='A'; break;

				case '-': out[in.length()-1-i]='-'; break;

				
				default: throw new IllegalArgumentException("unknown character: '"+in.charAt(i)+"'");
			}
		}
		return new String( out );
	}
	
	public static Random r = new Random();
	public static HashMap<String,String[]> c = new HashMap<String, String[]>();
	static{
		c.put( "A", new String[]{"A"} );
		c.put( "C", new String[]{"C"} );
		c.put( "G", new String[]{"G"} );
		c.put( "T", new String[]{"T"} );
		
		c.put( "N", new String[]{"A","C","G","T"} );
		
		c.put( "Y", new String[]{"C","T"} );
		c.put( "R", new String[]{"A","G"} );
		
		c.put( "W", new String[]{"A","T"} );
		c.put( "S", new String[]{"C","G"} );
		
		c.put( "K", new String[]{"G","T"} );
		c.put( "M", new String[]{"A","C"} );
		
		c.put( "B", new String[]{"C","G","T"} );
		c.put( "D", new String[]{"A","G","T"} );
		c.put( "H", new String[]{"A","C","T"} );
		c.put( "V", new String[]{"A","C","G"} );
	}
	
	private static HashMap<Character, int[]> help = new HashMap<Character, int[]>();
	
	/**
	 * Enum 
	 * 
	 * @author Jens Keilwagen
	 */
	public static enum Ambiguity {
		EXCEPTION,
		AMBIGUOUS,
		RANDOM;
	}
	
	public static char translate( CharSequence triplett, HashMap<String, Character> code, Ambiguity ambiguity ) {
		Character as = code.get(triplett);
		if( as == null ) {
			String[][] current = new String[3][];
			int anz = 0;
			for( int i = 0; i < 3; i++ ) {
				current[i] = c.get(triplett.subSequence(i,i+1));
				if( current[i] == null ) {
					throw new IllegalArgumentException( "Check nucleotide: " + triplett.subSequence(i,i+1) );
				}
				anz = Math.max(anz, current[i].length);
			}
			if( anz == 1 ) {
				as = code.get(current[0][0]+current[1][0]+current[2][0]);
			} else {
				help.clear();
				anz = 0;
				for( int i = 0; i < current[0].length; i++ ) {
					for( int j = 0; j < current[1].length; j++ ) {
						for( int k = 0; k < current[2].length; k++ ) {
							as = code.get(current[0][i]+current[1][j]+current[2][k]);
							int[] count = help.get(as);
							if( count == null ) {
								help.put(as, new int[]{1});
							} else {
								count[0]++;
							}
							anz++;
						}
					}					
				}
				if( help.size() > 1 ) {
					switch( ambiguity ) {
						case AMBIGUOUS:
							as='X';
							break;
						case RANDOM: 
							int idx = r.nextInt(anz), v=0;
							Iterator<Entry<Character,int[]>> it = help.entrySet().iterator();
							Entry<Character,int[]> e = null;
							while( it.hasNext() && v <= idx ) {
								e = it.next();
								v += e.getValue()[0];
							}
							as = e.getKey();
							as = Character.toLowerCase(as);
							break;
						case EXCEPTION:
							throw new IllegalArgumentException( triplett.toString() );
					}
				} //else 'as' is already set correctly
			}
		}
		return as;
	}
	

	public static boolean getUnambigious(StringBuffer seq) {
		boolean changed = false;
		for( int i = 0; i < seq.length(); i++ ) {
			String base = seq.substring(i, i+1).toUpperCase();
			String[] possibilities = c.get(base);
			if( possibilities.length > 1 ) {
				seq.replace(i, i+1, possibilities[r.nextInt(possibilities.length)]);
				changed=true;
			}
		}
		return changed;
	}
	
	public static String translate( int offset, String seq, HashMap<String, Character> code, boolean check, Ambiguity ambiguity ) throws IllegalArgumentException {
		StringBuffer aaSeqBuff = new StringBuffer();
		char as;
		boolean last = false;
		for( int j = offset; j+3 <= seq.length(); j+=3 ) {
			as = translate( seq.substring(j, j+3), code, ambiguity );
			if( check ) {
				if( as == '*' ) {
					last = true;
				} else if( last ) {
					throw new IllegalArgumentException();// aaSeqBuff.toString() + as );
				}
			}
			aaSeqBuff.append( as );
		}
		return aaSeqBuff.toString();
	}
	
	public static HashMap<String,String[]> getAlias( String fName, int oldNameIdx, int newNameIdx, int countIdx ) throws Exception {
		HashMap<String, String[]> res = new HashMap<String, String[]>();
		getAlias(res, fName, null, oldNameIdx, newNameIdx, countIdx);
		return res;
	}
	
	public static void getAlias( HashMap<String,String[]> res, String fName, String prefix, int oldNameIdx, int newNameIdx, int countIdx ) throws Exception {
		if( fName == null || !(new File(fName).exists())) {
			return;
		}
		if( prefix == null ) {
			prefix="";
		} else {
			prefix = prefix.toUpperCase() + "_";
		}
		BufferedReader r = new BufferedReader( new FileReader(fName) );
		String line;
		String[] split;
		while( (line=r.readLine()) != null ) {
			if( line.length()==0 || line.charAt(0)=='#' ) continue;
			split = line.split("\t");
			split[newNameIdx] = split[newNameIdx].toUpperCase();
			res.put(prefix+split[oldNameIdx].toUpperCase(), countIdx < 0 ? new String[]{ prefix+split[newNameIdx] } : new String[]{ prefix+split[newNameIdx], split[countIdx].split(",").length+"" });
		}
		r.close();
	}
	
	public static BufferedReader openGzOrPlain( String fName ) throws IOException {
		InputStream ins = new FileInputStream(fName);
		try {
			GZIPInputStream gz = new GZIPInputStream(ins);
			ins = gz;
		} catch( EOFException | ZipException ze ) {
			ins.close();
			ins = new FileInputStream(fName);
		}
		return new BufferedReader( new InputStreamReader(ins) );
	}
	
	public static HashMap<String,String> getFasta( String fName, int initSize ) throws Exception {
		return getFasta(fName, initSize, "([^_^\t^,^;^=^\"]+)|([^\t^,^;^=^\"]+_\\d+)");
	}
	
	public static HashMap<String,String> getFasta( String fName, int initSize, String seqIdRegex ) throws Exception {
		return getFasta(fName, initSize, "\\s.*", "", seqIdRegex);
	}
	
	public static HashMap<String,String> getFasta( String fName, int initSize, String regex, String replace, String seqIDRegex ) throws Exception {
		HashMap<String,String> seqs = new HashMap<String, String>(initSize);
		if( fName!=null ) {
			BufferedReader r = Tools.openGzOrPlain( fName );//new BufferedReader( new FileReader( fName ) );
			StringBuffer seq = new StringBuffer();
			String comment=null, line;
			while( (line=r.readLine()) != null ) {
				if( line.startsWith(">") ) {
					if( comment != null ) {
						add( seqs, comment, seq );
					}
					//clear
					seq.delete(0, seq.length());
					
					comment = line.substring(1);
					if( regex != null ) {
						comment=comment.replaceAll(regex, replace);
					}
					if( !comment.matches(seqIDRegex) ) {
						throw new IllegalArgumentException("Sequence ID ("+comment+") in fasta comment line (" + line +") does not match the regular expression for sequence IDs (" + seqIDRegex +")" );
					}
				} else {
					//add
					seq.append(line.trim().toUpperCase() );
				}
			}
			add( seqs, comment, seq );
			r.close();
		}
		return seqs;
	}
	
	private static void add( HashMap<String, String> seqs, String comment, StringBuffer seq ) {
		String old = seqs.get(comment), newSeq = seq.toString();
		if( old==null ) {
			seqs.put( comment, newSeq );
		} else {
			if( old.equals(newSeq) ) {			
				System.err.println("WARNING: duplicated sequence ignored: " + comment);
			} else {
				throw new IllegalArgumentException("At least two sequences with the same ID but different sequence: " + comment );
			}
		}
		
	}
	
	/**
	 * This method modifies a filter String that can be used in a {@link ScriptEngine}.
	 *  
	 * @param filter the user-specified filter String
	 * 
	 * @return the modified filter String
	 * 
	 * @see ScriptEngine
	 * @see #filter(ScriptEngine, String, HashMap)
	 */
	public static String prepareFilter( String filter ) {
		if( filter == null ) {
			filter = "";
		} else {
			filter = filter.trim();
		}
		filter = filter.replaceAll( " or ", " || " );
		filter = filter.replaceAll( " and ", " && " );
		return filter;
	}
	
	/**
	 * This method evaluates an <code>expression</code>.
	 * 
	 * @param engine the {@link ScriptEngine} to be used for filtering
	 * @param expression to be evaluated
	 * @param hash the hash containing the attributes in key-value-pairs
	 * 
	 * @return the evaluated expression as String
	 * 
	 * @throws ScriptException if the expression could not be evaluated properly
	 * 
	 * @see ScriptEngine
	 */
	public static String eval( ScriptEngine engine, String expression, HashMap<String,String> hash ) throws ScriptException {
		Bindings b = engine.createBindings();
		Iterator<Entry<String,String>> it = hash.entrySet().iterator();
		while( it.hasNext() ) {
			Entry<String,String> e = it.next();
			String key = e.getKey();
			if( expression.indexOf(key)>=0 ) {
				String val = e.getValue();
				b.put(key, val.equals("NA")?""+Double.NaN:val);
			}
		}			
		
		return engine.eval(expression, b).toString();
	}
	
	/**
	 * This method evaluates a filter String and returns a boolean.
	 * 
	 * @param engine the {@link ScriptEngine} to be used for filtering
	 * @param filter the filter String
	 * @param hash the hash containing the attributes in key-value-pairs
	 * @return the result of the evaluated filter String
	 * 
	 * @throws ScriptException if the filter String could not be evaluated properly
	 * 
	 * @see #prepareFilter(String)
	 * @see #eval(ScriptEngine, String, HashMap)
	 */
	public static boolean filter( ScriptEngine engine, String filter, HashMap<String,String> hash ) throws ScriptException {
		if( hash == null ) {
			return false;
		} else if( filter.length()== 0 ) {
			return true;
		} else {			
			String s = eval(engine, filter, hash);
			Boolean bool = Boolean.parseBoolean( s );
			return bool;
		}
	}
	
}
