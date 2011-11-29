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
import java.util.Date;
import java.util.GregorianCalendar;

/**
 * This class implements a {@link FileFilter} that accepts {@link File}s that were modified after the date that is given in the constructor.
 * 
 * @author Jens Keilwagen
 */
public class DateFileFilter implements FileFilter {

	private long time;
	private String desc;
	
	/**
	 * Creates an instance that accepts {@link File}s that were modified after the given year, month, ... .
	 * 
	 * @param year the year
	 * @param month the month
	 * @param dayOfMonth the day of the month
	 * @param hrs the hours
	 * @param min the minutes
	 * @param sec the seconds
	 * 
	 * @see GregorianCalendar#GregorianCalendar(int, int, int, int, int, int)
	 */
	public DateFileFilter( int year, int month, int dayOfMonth, int hrs, int min, int sec ) {
		this( new GregorianCalendar(year, month, dayOfMonth, hrs, min, sec).getTime() );
	}
	
	/**
	 * Creates an instance that accepts {@link File}s that were modified after <code>d</code>.
	 * 
	 * @param d the date used to decide whether a {@link File} will be accepted or not
	 */
	public DateFileFilter( Date d ) {
		time = d.getTime();
		desc = d.toString() + " (" + time +")";
	}
	
	public boolean accept(File f) {
		return f.isDirectory() || f.lastModified() >= time;
	}

	public String toString() {
		return desc;
	}
}

