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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;

import de.jstacs.io.FileManager;
import de.jstacs.utils.SafeOutputStream;

/**
 * This class is a VERY SIMPLE parser for Tex-Files to wiki.
 * 
 * @author Jens Keilwagen
 */
public class Tex2Wiki {

	private static final String HOME = "./supplementary/cookbook/";
	
	private static final String START = "\\begin{document}";
	private static final String END = "\\end{document}";
	
	private static final String CMD = "\\newcommand{";
	private static final String INPUT = "\\input{";
	
	private static final String OPEN = "{";
	private static final String CLOSE = "}";


	private static boolean isOkay( char c ) {
		return ('a' <= c && c <='z') || ('A' <= c && c <='Z');
	}
	
	private static int findOpeningTag( String tag, int startIdx, StringBuffer s ) {
		return s.indexOf( tag, startIdx );
	}
	
	private static HashMap<Character, Character> tag = new HashMap<Character, Character>();
	
	static {
		tag.put( '{', '}' );
		tag.put( '[', ']' );
		tag.put( '(', ')' );
	}
	
	private static char getEndTag( char c ) {
		Character end = tag.get( c );
		if( end == null ) {
			return c;
		}
		return end;
	}
	
	private static int findClosingTag( int startIdx, StringBuffer s ) throws Exception {
		char startTag = s.charAt(startIdx), endTag = getEndTag( startTag );
		String eTag = ""+endTag, sTag = "" + startTag;
		startIdx++;
		return findClosingTag( startIdx, sTag, eTag, s );
	}
	
	private static int findClosingTag( int startIdx, String sTag, String eTag, StringBuffer s ) throws Exception {
		int hh = startIdx;
		int anz = 1, idx, h;
		do {
			//System.out.println( "a\t" + startIdx + " " + anz + " " + s.substring( startIdx, Math.min( s.length(), startIdx+25 ) ).replace( "\n", " " ) );
			idx = s.indexOf( eTag, startIdx );
			anz--;
			while( (h = s.indexOf( sTag, startIdx )) < idx && h > 0 ) {
				anz++;
				//System.out.println( "xxx " + startIdx + " " + h + " " + anz + " " + s.substring( h, Math.min( s.length(), h+25 ) ).replace( "\n", " " ) );
				startIdx = h+1;
			}
			//System.out.println( startIdx + " " + anz );
			startIdx = idx+1;
		} while( idx > 0 && anz > 0 );
		if( idx < 0 ) {
			System.out.println( "CLOSE?\t" + s.substring( hh ));
			throw new Exception( "could not parse: did not find closing \"" + eTag + "\"" );
		}
		return idx;
	}
	
	private static void add( StringBuffer s, HashMap<String, Replacement> hash ) throws Exception {
		int idx1, idx2, idx3, idx4, startIdx = 0, params;
		String cmd, rep;
		while( (idx1=s.indexOf( CMD, startIdx )) >= 0 ) {
			idx2 = s.indexOf( "}", idx1 );
			idx3 = findOpeningTag( OPEN, idx2+1, s );
			idx4 = findClosingTag( idx3, s );
			
			cmd = s.substring( idx1+CMD.length(), idx2 ).trim();
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
			hash.put( cmd, new SimpleReplacement( params, rep ) );
			
			startIdx = idx4+1;
		}
		
		//recursive -> input
		startIdx = 0;
		StringBuffer sb;
		while( (idx1=s.indexOf( INPUT, startIdx )) >= 0 ) {
			idx2 = s.indexOf( "}", idx1 );
			rep = s.substring( idx1+INPUT.length(), idx2 );
			sb = FileManager.readFile( new File( HOME + rep ) );
			add( sb, hash );			
			startIdx = idx2+1;
		}
	}
	
	private static HashMap<String,Replacement> hash = new HashMap<String,Replacement>(50);
	
	/**
	 * This is the main of the program.
	 * 
	 * @param args the arguments
	 * 
	 * @throws Exception if something went wrong
	 */
	public static void main( String[] args ) throws Exception {
		hash.put( "\\newcommand", new NewCommandReplacement() );
		hash.put( "\\renewcommand", new RenewCommandReplacement() );
		
		hash.put( "\\section", new SimpleReplacement( 1, "= #1 =") );
		hash.put( "\\subsection", new SimpleReplacement( 1, "== #1 ==") );
		hash.put( "\\subsubsection", new SimpleReplacement( 1, "=== #1 ===") );
		hash.put( "\\emph", new SimpleReplacement( 1, "''#1''") );
		hash.put( "\\textit", new SimpleReplacement( 1, "''#1''") );
		hash.put( "\\textbf", new SimpleReplacement( 1, "'''#1'''") );
		hash.put( "\\textsc", new SimpleReplacement( 1, "<span style=\"font-variant: small-caps;\">#1</span>") );//TODO
		hash.put( "\\href", new SimpleReplacement( 2, "[#1 #2]") );
		hash.put( "\\url", new SimpleReplacement( 1, "[#1 #1]") );
		hash.put( "\\caption", new SimpleReplacement( 1, "#1" ) );
		hash.put( "\\lstinline", new SimpleReplacement(1, "<code>#1</code>") );
		hash.put( "\\label", new SimpleReplacement( 1, "\n<span id=\"#1\"> </span>" ) );
		hash.put( "\\ref", new SimpleReplacement( 1, "[[##1 (link)]]" ) );//TODO
		hash.put( "\\lstset", new SimpleReplacement( 1, "" ));
		hash.put( "\\textcolor", new SimpleReplacement( 2, "<font color=#1>#2</font>" ) );
		
		hash.put( "\\newcounter", new NewCounterReplacement() );
		hash.put( "\\setcounter", new SetCounterReplacement() );
		hash.put( "\\stepcounter", new StepCounterReplacement() );
		hash.put( "\\addtocounter", new AddToCounterReplacement() );
		hash.put( "\\code", new CodeReplacement() );
		
		hash.put( "\\begin", new EnvironmentReplacement() );

		
		createWiki( false, "defs", System.out );
		
		createWiki( true, "preface", System.out );
		createWiki( true, "data", System.out );
		createWiki( true, "infrastructure", System.out );
		createWiki( true, "sequenceScore", System.out );
		createWiki( true, "classifier", System.out );
		createWiki( true, "optimization", System.out );
		createWiki( true, "goodies", System.out );
		createWiki( true, "recipes", System.out );
	}
	
	private static ArrayList<String> list = new ArrayList<String>();
	
	private static void createWiki( boolean create, String file, OutputStream out ) throws Exception {
		SafeOutputStream sos = SafeOutputStream.getSafeOutputStream( out );
		System.out.println("--------------------------------------------------");
		System.out.println( file );
		StringBuffer s = FileManager.readFile( new File( HOME + file + ".tex" ) );
		
		//remove comments
		int start = 0, end;
		while( (start = s.indexOf("%", 0)) >= 0 ) {
			if( s.charAt(start-1) != '\\' ) {
				end = s.indexOf("\n", start+1);
				if( end < 0 ) {
					end = s.length();
				}
				s.delete(start, end);
			} else {
				start++;
			}
		}
		
		start = s.indexOf( START );
		end = s.indexOf( END );
		StringBuffer wiki;
		wiki = new StringBuffer( s.toString()
				//naive replacements!
				.replaceAll( "\\\\\\\\", "\n" )
				.replaceAll( "~", " " )
				.replaceAll( "\\\\&", "&" )
				.replaceAll( "''", "&quot;" )
				.replaceAll( "``", "&quot;" )
		);
		if( start >= 0 && end >= 0 ) {
			wiki = new StringBuffer( wiki.substring( start+START.length(), end ) ); 
		} else if( start < 0 && end < 0 ) {
		} else {
			throw new Exception("could not localize part for parsing");
		}
		
		Replacement e;
		String cmd;
		int idx1 = 0, idx2;
		while( (idx1 = wiki.indexOf( "\\", idx1 ) ) >= 0 ){
			//extract command
			idx2 = idx1+1;
			while( isOkay( wiki.charAt( idx2 ) ) ) {
				idx2++;
			}
			cmd = wiki.substring( idx1, idx2 );
			e = hash.get( cmd );
			
			try {
				//replace
				if( e == null ) {
					idx1++; //skip
					sos.writeln( "PROBLEM\t" + cmd );
				} else {
					list.clear();
					idx1 = e.replace( wiki, idx1, idx2 );
				}
			} catch( Exception ex ) {
				System.out.println( "ERROR\t" + cmd + "\t" + ex.getMessage() );
				ex.printStackTrace();
				idx1++;
			}
		}
		
		idx1 = 0;
		boolean open = true;
		while( (idx1 = wiki.indexOf( "$", idx1 ) ) >= 0 ){
			
			if( idx1 > 0 && wiki.charAt(idx1-1) != '\\' ) {
				wiki.replace( idx1, idx1+1, open ? "<math>" : "</math>" );
				idx1 += open?6:7;
				open = !open;
			} else {
				idx1++;
			}
		}
		
		if( create ) {
			FileManager.writeFile( new File( HOME + file + ".wiki" ), wiki );
		}
	}
	
	private static int fillParams( StringBuffer wiki, int idx2, int anz ) throws Exception {
		int idx4 = idx2, idx3;
		for( int i = 0; i < anz; i++ ) {
			idx3 = idx4;
			idx4 = findClosingTag( idx3, wiki );
			list.add( wiki.substring( idx3+1, idx4 ) );
			idx4++;
		}
		return idx4;
	}
	

	private static interface Replacement {
		
		/**
		 * The method performs the replacement of a specific TeX command to a specific wiki command.
		 *  
		 * @param wiki the {@link StringBuffer} to be modified (from tex to wiki)
		 * @param start the index where the command starts
		 * @param startParams the index where the parameters of the command start
		 * 
		 * @return the index where to further replace, if recursive replacements are necessary you have to return <code>start</code>
		 * 
		 * @throws Exception if something went wrong
		 */
		public abstract int replace( StringBuffer wiki, int start, int startParams ) throws Exception;
	}
	
	private static class SimpleReplacement implements Replacement {

		int anz;
		String template;
		
		public SimpleReplacement( int anz, String template ) {
			this.anz = anz;
			this.template = template;
		}
		
		public int replace( StringBuffer wiki, int start, int startParams ) throws Exception {
			int end = fillParams( wiki, startParams, anz );			
			String result = template;
			int anz = list.size();
			//System.out.println( "BEFORE\t\"" + wiki.substring(start, end) + "\"" );
			//System.out.println( list );
			while( anz > 0 ) {
				//System.out.println( "AFTER\t" + (anz == list.size() ? "" : list.get( anz )) + "\t\"" + result + "\"" );
				String h = "#" + anz;
				--anz;
				while( result.indexOf(h) >= 0 ) {
					result = result.replace( h, list.get( anz ) );
				}
			}
			//System.out.println( "AFTER\t" + (anz == list.size() ? "" : list.get( anz )) + "\t\"" + result + "\"" );
			wiki.replace( start, end, result );
			return start;
		}
	}
	
	private static class NewCommandReplacement implements Replacement {		
		public int replace( StringBuffer wiki, int start, int startParams ) throws Exception {
			int end = fillParams( wiki, startParams, 1 );
			delete();
			int anz = wiki.charAt(end) == '[' ? 2: 1;
			int end2 = fillParams( wiki, end, anz );
			if( anz == 2 ) {
				anz = Integer.parseInt( list.get(1) );
			} else {
				anz = 0;
			}
			hash.put( list.get(0), new SimpleReplacement( anz, list.get(list.size()-1) ) );

			wiki.delete( start, end2 );
			return start;
		}
		
		private void delete() {}
	}
	
	private static class RenewCommandReplacement extends NewCommandReplacement {		
		private void delete() {
			hash.remove( list.get(0) );
		}
	}
	
	private static HashMap<String,int[]> counter = new HashMap<String, int[]>();
	
	private static class NewCounterReplacement implements Replacement {
		public int replace( StringBuffer wiki, int start, int startParams ) throws Exception {
			int end = fillParams( wiki, startParams, 1 );
			//System.out.println( list.get(0) + "\t" + 0 );	
			counter.put( list.get(0), new int[]{ 0 } );
			wiki.delete( start, end );
			return start;
		}
	}
	
	private static class SetCounterReplacement implements Replacement {
		public int replace( StringBuffer wiki, int start, int startParams ) throws Exception {
			int end = fillParams( wiki, startParams, 2 );
			int[] c = counter.get( list.get(0) );
			//System.out.print( list.get(0) + "\t" + c[0] );
			int v;
			try{
				v = Integer.parseInt(list.get(1));
			} catch( Exception e ) {
				String s = list.get(1);
				v = counter.get( s.substring(8,s.length()-1) )[0];
			}
			if( c == null ) {
				counter.put( list.get(0), new int[]{ v } );
			} else {
				c[0] = v;
			}
			//System.out.println( " -> " + counter.get( list.get(0) )[0] );
			wiki.delete( start, end );
			return start;
		}
	}
	
	private static class StepCounterReplacement implements Replacement {
		public int replace( StringBuffer wiki, int start, int startParams ) throws Exception {
			int end = fillParams( wiki, startParams, 1 );
			//System.out.print( list.get(0) + "\t" + counter.get( list.get(0) )[0] );
			counter.get( list.get(0) )[0]++;
			//System.out.println( " -> " + counter.get( list.get(0) )[0] );
			wiki.delete( start, end );
			return start;
		}
	}
	
	private static class AddToCounterReplacement implements Replacement {
		public int replace( StringBuffer wiki, int start, int startParams ) throws Exception {
			int end = fillParams( wiki, startParams, 2 );
			int[] c = counter.get( list.get(0) );
			//System.out.print( list.get(0) + "\t" + c[0] );
			c[0] += Integer.parseInt(list.get(1));
			//System.out.println( " -> " + c[0] );
			wiki.delete( start, end );	
			return start;
		}
	}

	private static class CodeReplacement implements Replacement {
		public int replace( StringBuffer wiki, int start, int startParams ) throws Exception {
			int end = fillParams( wiki, startParams, 1 );
			
			int off = counter.get( "off" )[0];
			StringBuffer _new = new StringBuffer();
			_new.append( "<source lang=\"java5\" enclose=\"div\">\n" );
			String s = "\\codefile";
			do {
				s = ((SimpleReplacement) hash.get( s )).template;
			} while( s.charAt(0) == '\\' );
			
			BufferedReader r = new BufferedReader( new FileReader( HOME + s ) );
			for( int i = 1; i < off; i++ ) {
				r.readLine();
			}
			int anz = Integer.parseInt( list.get(0) );
			String prefix = null, line;
			for( int i = 0; i <= anz; i++ ) {
				line = r.readLine();
				if( i == 0 ) {
					int p = 0;
					while( p < line.length() && line.charAt(p) == '\t' ) {
						p++;
					}
					prefix = line.substring(0,p);
				} 
				
				if( line.startsWith(prefix) ) {
					line = line.substring(prefix.length());
				} else {
					if( line.length() > 0 ) {
						System.out.println( "WARNING " + this.getClass().getSimpleName() + "\t" + (off+i) + "\t" + s  + "\t\"" + line + "\"" );
					}
				}
				_new.append( line );
				_new.append( "\n" );
			}
			r.close();
			_new.append( "</source>\n" );
			
			wiki.replace( start, end, _new.toString() );
			return start+_new.length();
		}
	}
	
	private static class EnvironmentReplacement implements Replacement {
		public int replace( StringBuffer wiki, int start, int startParams ) throws Exception {
			int end = fillParams( wiki, startParams, 1 );
			int end2 = findClosingTag( end, "\\begin{"+list.get(0)+"}", "\\end{"+list.get(0)+"}", wiki );
			
			int res = start;
			
			StringBuffer _new = new StringBuffer();
			String s = list.get(0), h;
			if( s.equals("itemize") ) {
				h = wiki.substring( end, end2 );
				h = h.replaceAll( "[ \\t]*\\\\item", "*" );
				_new.append( h );
			} else if( s.equals( "figure" ) ){
				h = wiki.substring( end, end2 );
				System.out.println( h );
				_new.append("[[File:TODO|thumb|");//TODO
				_new.append( h );
				_new.append( "]]\n" );
			}else if( s.equals( "align*" ) || s.equals( "align" ) || 
					s.equals( "equation*" ) || s.equals( "equation" ) ||
					s.equals( "eqnarray*" ) || s.equals( "eqnarray" )){
				String s2 = wiki.substring( end, end2 );
				String s3 = s;
				if(s.endsWith( "*" )){
					s3 = s.substring( 0, s.length()-1 );
				}
				_new.append( "\n<math>" );
				_new.append( "\\begin{" );
				res += _new.length();
				_new.append( s3 );
				_new.append( "}\n" );
				_new.append( s2 );
				_new.append( "\n\\end{" );
				_new.append( s3 );
				_new.append( "}" );
				_new.append( "</math>\n" );
			}else if( s.equals( "lstlisting" ) ){
				String s2 = wiki.substring( end,end2 );
				String lang = "java5";
				if(s2.startsWith( "[" )){
					int end3 = s2.indexOf( "]" );
					if(s2.startsWith( "[language=" )){
						lang = s2.substring( 10, end3 );
					}
					s2 = s2.substring( end3+1 );
				}
				_new.append( "<source lang=\""+lang+"\" enclose=\"div\">" );
				_new.append( s2 );
				_new.append( "</source>" );
			}else{
				throw new Exception( s );
			}

			wiki.replace( start, end2+s.length()+2+4, _new.toString() );
			return res;
		}
	}
}