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

package de.jstacs.utils;

import java.io.IOException;
import java.io.OutputStream;

/**
 * This interface is the framework for stopping the time of anything.
 * 
 * @author Jens Keilwagen
 */
public abstract class Time {

	/**
	 * Creates a new time object and starts the clock.
	 */
	public Time() {
		reset();
	}

	/**
	 * Returns the elapsed time since invoking the constructor.
	 * 
	 * @return the elapsed time in seconds
	 */
	public abstract double getElapsedTime();

	/**
	 * Restarts the time stopping.
	 */
	public abstract void reset();
	
	/**
	 * This method tries to return a {@link UserTime} instance, if not possible (due to native code) it returns a {@link RealTime} instance.
	 * 
	 * @param out a stream that allows to write a warning if a {@link RealTime} instance is returned; can be <code>null</code>
	 * 
	 * @return a {@link UserTime} or {@link RealTime} instance
	 * 
	 * @throws IOException forwarded from {@link OutputStream#write(byte[])} 
	 */
	public static Time getTimeInstance( OutputStream out ) throws IOException {
		Time t = null;
		try {
			t = new UserTime();
		} catch ( Error err ) {
			if( out != null ) {
				out.write( ("Warning: Could not load UserTime. Using RealTime instead.\n").getBytes() );
			}
			t = new RealTime();
		}
		return t;
	}

}
