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

import de.jstacs.io.RegExFilenameFilter;

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
		RegExFilenameFilter filter = new RegExFilenameFilter("java",true,true,".*\\.java");
		File dir = new File( "./" );
		
		check( filter, dir );
	}
	
	private static void check( FilenameFilter filter, File dir ) throws Exception {
		File[] f = dir.listFiles(filter);
		for( File g : f ) {
			if( g.isDirectory() ) {
				check( filter, g );
			} else {
				BufferedReader r = new BufferedReader( new FileReader( g ) );
				String line = r.readLine();
				r.close();
				if( !line.startsWith("/*") ) {
					//copy preamble
					File h = new File( g.getAbsolutePath() + ".preamble" );
					BufferedInputStream bis = new BufferedInputStream( new FileInputStream(g) );
					BufferedOutputStream bos = new BufferedOutputStream( new FileOutputStream(h) );
					for( String s : preamble ) {
						bos.write( s.getBytes() );
						bos.write( '\n' );
					}
					int i;
					while( (i=bis.read()) > 0 ) {
						bos.write(i);
					}
					bis.close();
					bos.close();
					/*
					BufferedWriter w = new BufferedWriter( new FileWriter( h ) );
					for( String s : preamble ) {
						w.write( s );
						w.newLine();
					}
					w.write( line );
					w.newLine();
					while( (line=r.readLine()) != null ) {
						w.write( line );
						w.newLine();
					}
					w.close();
					r.close();
					*/
					
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
				}
			}
		}
	}
}
