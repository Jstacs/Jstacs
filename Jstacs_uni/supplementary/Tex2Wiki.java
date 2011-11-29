/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package supplementary;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Map.Entry;

import de.jstacs.io.FileManager;

/**
 * This class is a VERY SIMPLE parser for Tex-Files to wiki.
 * 
 * @author Jens Keilwagen
 */
public class Tex2Wiki {

	private static final String HOME = "./supplementary/";
	
	private static final String START = "\\begin{document}";
	private static final String END = "\\end{document}";
	
	private static final String CMD = "\\newcommand{";
	private static final String INPUT = "\\input{";
	
	private static final String OPEN = "{";
	private static final String CLOSE = "}";


	private static boolean isOkay( char c ) {
		return ('a' <= c && c <='z') || ('A' <= c && c <='Z');
	}
	
	private static String getReplacement( String template, ArrayList<String> params ) {
		int anz = params.size();
		while( anz > 0 ) {
			template = template.replace( "#" + anz, params.get( --anz ) );			
		}
		return template;
	}
	
	private static int findOpeningTag( String tag, int startIdx, StringBuffer s ) {
		return s.indexOf( tag, startIdx );
	}
	
	private static int findClosingTag( String startTag, String endTag, int startIdx, StringBuffer s ) throws Exception {
		int anz = 1, idx, h;
		do {
			//System.out.println( "a\t" + startIdx + " " + anz + " " + s.substring( startIdx, Math.min( s.length(), startIdx+25 ) ).replace( "\n", " " ) );
			idx = s.indexOf( endTag, startIdx );
			anz--;
			while( (h = s.indexOf( startTag, startIdx )) < idx && h > 0 ) {
				anz++;
				//System.out.println( "xxx " + startIdx + " " + h + " " + anz + " " + s.substring( h, Math.min( s.length(), h+25 ) ).replace( "\n", " " ) );
				startIdx = h+1;
			}
			//System.out.println( startIdx + " " + anz );
			startIdx = idx+1;
		} while( idx > 0 && anz > 0 );
		if( idx < 0 ) {
			throw new Exception( "could not parse: did not find closing \"" + endTag + "\"" );
		}
		return idx;
	}
	
	private static void addCommands( StringBuffer s, Hashtable<String, Entry<Integer,String>> hash ) throws Exception {
		int idx1, idx2, idx3, idx4, startIdx = 0, params;
		String cmd, rep;
		while( (idx1=s.indexOf( CMD, startIdx )) >= 0 ) {
			idx2 = s.indexOf( "}", idx1 );
			idx3 = findOpeningTag( "{", idx2+1, s );
			idx4 = findClosingTag( "{", "}", idx3+1, s );
			
			cmd = s.substring( idx1+CMD.length(), idx2 );
			rep = s.substring( idx3+1,idx4 );
			
			if( idx2+1 == idx3 ){
				params = 0;
			} else {
				if( s.charAt( idx2+1 ) == '[' && s.charAt( idx3-1 ) == ']' ) {
					params = Integer.parseInt( s.substring( idx2+2, idx3-1 ) );
				} else {
					throw new Exception( "could not parse command" );
				}
			}
			hash.put( cmd, new MyEntry( params, rep ) );
			
			startIdx = idx4+1;
		}
		
		//recursive -> input
		startIdx = 0;
		StringBuffer sb;
		while( (idx1=s.indexOf( INPUT, startIdx )) >= 0 ) {
			idx2 = s.indexOf( "}", idx1 );
			rep = s.substring( idx1+INPUT.length(), idx2 );
			sb = FileManager.readFile( new File( HOME + rep ) );
			addCommands( sb, hash );			
			startIdx = idx2+1;
		}
	}
	
	/**
	 * This is the main of the program.
	 * 
	 * @param args the arguments
	 * 
	 * @throws Exception if something went wrong
	 */
	public static void main( String[] args ) throws Exception {

		Hashtable<String,Entry<Integer,String>> hash = new Hashtable<String,Entry<Integer,String>>(50);
		hash.put( "\\section", new MyEntry( 1, "= #1 =") );
		hash.put( "\\subsection", new MyEntry( 1, "== #1 ==") );
		hash.put( "\\emph", new MyEntry( 1, "''#1''") );
		hash.put( "\\textit", new MyEntry( 1, "''#1''") );
		hash.put( "\\textbf", new MyEntry( 1, "'''#1'''") );
		hash.put( "\\textsc", new MyEntry( 1, "<span style=\"font-variant: small-caps;\">#1</span>") );//TODO
		hash.put( "\\href", new MyEntry( 2, "[#1 #2]") );
		
		StringBuffer s = FileManager.readFile( new File( HOME + "Mixtures.tex" ) ); //TODO args[0]
		int start = s.indexOf( START ) + START.length();
		int end = s.indexOf( END );
		
		addCommands( new StringBuffer( s.substring( 0, start ).replace( "\\\\", "\n" ) ), hash );
		
		StringBuffer wiki = new StringBuffer( s.substring( start, end ).replace( "\\\\", "\n" ) ); 
		
		Entry<Integer,String> e;
		String cmd;
		int idx1 = 0, idx2, idx3, idx4, anz;
		ArrayList<String> list = new ArrayList<String>();
		while( (idx1 = wiki.indexOf( "\\", idx1 ) ) >= 0 ){
			//extract command
			idx2 = idx1+1;
			while( isOkay( wiki.charAt( idx2 ) ) ) {
				idx2++;
			}
			cmd = wiki.substring( idx1, idx2 );
			e = hash.get( cmd );
			
			//replace
			if( e == null ) {
				idx1++; //skip
				//throw new Exception( "Could not parse command: " + cmd );//TODO
			} else {
				anz = e.getKey();
				list.clear();
				idx4 = idx2;
				for( int i = 0; i < anz; i++ ) {
					idx3 = findOpeningTag( OPEN, idx4, wiki );
					idx4 = findClosingTag( OPEN, CLOSE, idx3+1, wiki );
					list.add( wiki.substring( idx3+1, idx4 ) );
					idx4++;
				}
				if( wiki.charAt( idx1-1 ) == '{' && wiki.charAt( idx4 ) == '}' ) {
					idx1--;
					idx4++;
				}
				wiki.replace( idx1, idx4, getReplacement( e.getValue(), list ) );
			}
		}
		System.out.println( wiki );	
	}

	private static class MyEntry implements Entry<Integer,String> {

		private int key;
		private String value;
		
		public MyEntry( int key, String val ) {
			this.key = key;
			this.value = val;
		}
		
		public Integer getKey() {
			return key;
		}

		public String getValue() {
			return value;
		}

		public String setValue( String value ) {
			String old = this.value;
			this.value = value;
			return old;
		}	
	}	
}
