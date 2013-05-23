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

package de.jstacs.io;

import java.io.File;
import java.io.FilenameFilter;
import java.util.regex.Pattern;

import javax.swing.filechooser.FileFilter;

/**
 * A simple filter on {@link File}s that accepts {@link File}s with a specific regular expression
 * in the filename.
 * 
 * @author Jens Keilwagen
 * 
 * @see String#matches(String)
 * @see java.util.regex.Pattern
 */
public class RegExFilenameFilter extends FileFilter implements java.io.FileFilter, FilenameFilter {

	private String desc;

	private Pattern[] regex;

	private Directory dir;
	private boolean ignoreCase;
	
	/**
	 * A switch to decide whether the results of a {@link RegExFilenameFilter} are forbidden, allowed or required to be a directory.
	 * 
	 * @author Jens Keilwagen
	 */
	public static enum Directory {
		/**
		 * Directories will not show up when filtering.
		 */
		FORBIDDEN,
		/**
		 * Directories can not show up when filtering.
		 */
		ALLOWED,
		/**
		 * Only directories will not show up when filtering.
		 */
		REQUIRED
	}

	/**
	 * Creates a new {@link RegExFilenameFilter} with given
	 * regular expressions <code>regex</code> to be found in the
	 * file name of the {@link File}s to be filtered.
	 * 
	 * @param desc
	 *            the description of the filter, e.g. &quot;text-files
	 *            (*.txt)&quot;
	 * @param dir
	 *            a switch whether the results are forbidden, allowed or required to be a directory
	 * @param ignoreCase
	 *            indicates whether to ignore the case of the file names or not
	 * @param regex
	 *            an array of regular expressions; at least one regular expression has to match
	 *            the file name
	 */
	public RegExFilenameFilter( String desc, Directory dir, boolean ignoreCase, String... regex ) {
		this.desc = desc;
		this.dir = dir;
		this.ignoreCase = ignoreCase;
		this.regex = new Pattern[regex.length];
		for( int i = 0; i < regex.length; i++ ) {
			this.regex[i] = Pattern.compile( ignoreCase ? regex[i].toLowerCase() : regex[i] );
		}
	}

	/* (non-Javadoc)
	 * @see javax.swing.filechooser.FileFilter#accept(java.io.File)
	 */
	@Override
	public boolean accept( File arg0 ) {
		boolean isDir = arg0.isDirectory();
		if( isDir ) {
			return dir != Directory.FORBIDDEN;
		} else {
			if( dir == Directory.REQUIRED ) {
				return false;
			} else {
				String name = arg0.getName();
				if( ignoreCase ) {
					name = name.toLowerCase();
				}
				int i = 0;
				while( i < regex.length && !regex[i].matcher( name ).matches() ) {
						i++;
				}
				return i < regex.length;
			}
		}
	}

	/* (non-Javadoc)
	 * @see java.io.FilenameFilter#accept(java.io.File, java.lang.String)
	 */
	public boolean accept( File dir, String name ) {
		return accept( new File(dir,name) );
	}
	

	/* (non-Javadoc)
	 * @see javax.swing.filechooser.FileFilter#getDescription()
	 */
	@Override
	public String getDescription() {
		return desc;
	}
	
	/**
	 * Returns a representation of all used regular expressions.
	 * 
	 * @return a representation of all used regular expressions
	 */
	public String getRegEx() {
		StringBuffer sb = new StringBuffer();
		for( Pattern s : regex ) {
			sb.append( s.pattern() + "\n" );
		}
		return sb.toString();
	}
	
	public String toString() {
		return getRegEx();
	}
}
