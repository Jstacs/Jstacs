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
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Random;

import de.jstacs.tools.Protocol;

/**
 * Some methods used for GeMoMa.
 * 
 * @author Jens Keilwagen
 */
public class Tools {
	
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
		while( 3*i < MAX_INTRON_LENGTH //TODO  for runtime reasons (in N-stretches) 
				&& startPos >= 0 && startPos+3 <= chr.length() && c != c1 && c != c2 ) {
			codon = chr.subSequence(startPos,startPos+3);
			if( !forward ) {
				codon = Tools.rc(codon);
			}
			if( backup != null && i < backup.length() ) {
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
			selected.put(line.substring(0,idx).toUpperCase(), line.substring(second));
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
		if( fName == null || !(new File(fName).exists())) {
			return null;
		}
		HashMap<String, String[]> res = new HashMap<String, String[]>();
		BufferedReader r = new BufferedReader( new FileReader(fName) );
		String line;
		String[] split;
		while( (line=r.readLine()) != null ) {
			split = line.split("\t");
			split[newNameIdx] = split[newNameIdx].toUpperCase();
			res.put(split[oldNameIdx].toUpperCase(), countIdx < 0 ? new String[]{ split[newNameIdx] } : new String[]{ split[newNameIdx], split[countIdx].split(",").length+"" });
		}
		r.close();
		return res;
	}
	
	public static HashMap<String,String> getFasta( String fName, int initSize, char c ) throws Exception {
		HashMap<String,String> seqs = null;
		if( fName!=null ) {
			seqs = new HashMap<String, String>(initSize);
			BufferedReader r = new BufferedReader( new FileReader( fName ) );
			StringBuffer seq = new StringBuffer();
			String comment=null, line;
			while( (line=r.readLine()) != null ) {
				if( line.startsWith(">") ) {
					if( comment != null ) {
						add( seqs, comment, seq );
					}
					//clear
					seq.delete(0, seq.length());
					int idx = line.indexOf(c);
					comment = line.substring(1,idx<0?line.length():idx);				
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
		seqs.put( comment, seq.toString() );
	}
}