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
package supplementary;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;

import de.jstacs.io.RegExFilenameFilter;
import de.jstacs.io.RegExFilenameFilter.Directory;

public class PreambleChecker {

	private static final String[] preamble = {
		"/*",
		" * This file is part of Jstacs.",
		" *",
		" * Jstacs is free software: you can redistribute it and/or modify it under the",
		" * terms of the GNU General Public License as published by the Free Software",
		" * Foundation, either version 3 of the License, or (at your option) any later",
		" * version.",
		" *",
		" * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY",
		" * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR",
		" * A PARTICULAR PURPOSE. See the GNU General Public License for more details.",
		" *",
		" * You should have received a copy of the GNU General Public License along with",
		" * Jstacs. If not, see <http://www.gnu.org/licenses/>.",
		" *",
		" * For more information on Jstacs, visit http://www.jstacs.de",
		" */"
	};
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		File dir = new File( "./" );
		
		System.out.println( "java: " + check( new RegExFilenameFilter("java",Directory.ALLOWED,true,".*\\.java"), dir, null, null ) );
		System.out.println( "html: " + check( new RegExFilenameFilter("html",Directory.ALLOWED,true,".*\\.html"), dir, "<!--", "--!>", "<head>" ) );
	}
	
	private static int check( FilenameFilter filter, File dir, String addStart, String addEnd, String... before ) throws Exception {
		File[] f = dir.listFiles(filter);
		int anz = 0;
		for( File g : f ) {
			if( g.isDirectory() ) {
				anz += check( filter, g, addStart, addEnd, before );
			} else {
				BufferedReader r = new BufferedReader( new FileReader( g ) );
				String line;
				int i = 0, k = 0;
				while( (line = r.readLine()) != null && i < before.length ) {
					if( line.startsWith( before[i] ) ) {
						i++;
					}
					k += line.length() + 1;
				}
				boolean missingPreamble;
				if( addStart == null ) {
					missingPreamble = line != null && !line.startsWith("/*");
				} else {
					missingPreamble = line != null && !line.startsWith( addStart );
				}
				r.close();
				if( missingPreamble ) {
					//copy preamble
					File h = new File( g.getAbsolutePath() + ".preamble" );
					BufferedInputStream bis = new BufferedInputStream( new FileInputStream(g) );
					BufferedOutputStream bos = new BufferedOutputStream( new FileOutputStream(h) );
					int j = 0;
					while( j < k ) {
						bos.write(bis.read());
						j++;
					}
					if( addStart != null ) {
						write( bos, addStart );
					}
					for( String s : preamble ) {
						write( bos, s );
					}
					write( bos, addEnd );
					while( (i=bis.read()) > 0 ) {
						bos.write(i);
					}
					bis.close();
					bos.close();
					
					//rename
					boolean b;
					b = g.delete();
					System.out.print( b + "\t" );
					b = h.renameTo( g );
					System.out.print( b + "\t" );
					System.out.println( g );
					if( !b ) {
						h.delete();
					}
					
					anz++;
				}
			}
		}
		return anz;
	}
	
	private static void write( BufferedOutputStream bos, String s ) throws IOException {
		if( s != null ) {
			bos.write( s.getBytes() );
			bos.write( '\n' );
		}
	}
}
