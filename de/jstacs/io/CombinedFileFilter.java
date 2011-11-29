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
import java.io.FileFilter;

/**
 * This class allows to combine several {@link FileFilter}s.
 * 
 * @author Jens Keilwagen
 *
 * @see RegExFilenameFilter
 * @see DateFileFilter
 */
public class CombinedFileFilter implements FileFilter {

	private int anz;
	private FileFilter[] filter;
	
	/**
	 * Creates an instance that accepts a {@link File} if at least <code>minAccepted</code> filters accept the {@link File}.
	 * 
	 * @param minAccepted the minimal number of filters that has to accept a {@link File}
	 * @param filter the filters internally used
	 */
	public CombinedFileFilter( int minAccepted, FileFilter... filter ) {
		if( minAccepted < 0 || minAccepted > filter.length ) {
			throw new IllegalArgumentException( "Check the value for the parameter anz" );
		}
		this.anz = minAccepted;
		this.filter = filter.clone();
	}
	
	public boolean accept(File pathname) {
		int accepted = 0, i = 0;
		for( ; i < filter.length; i++ ) {
			accepted += filter[i].accept(pathname) ? 1 : 0;
		}
		return accepted >= anz;
	}
	
	public String toString() {
		String desc = "filter (to be accepted: " + anz + "):\n";
		for( int i = 0; i < filter.length; i++ ) {
			desc += filter[i].toString() + "\n";
		}
		return desc;
	}
}
