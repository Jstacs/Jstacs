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

/**
 * This is an implementation of {@link Time} that uses a native method.
 * 
 * <br>
 * <br>
 * 
 * If you like to use this class you must set the VM option
 * <code>-Djava.library.path</code> to the directory where the native library
 * resides.<br>
 * The native library is
 * <ul>
 * <li> <code>UserTime.dll</code> for Windows,</li>
 * <li> <code>libUserTime.so</code> for Linux,</li>
 * <li> <code>libUserTime.jnilib</code> for Mac OS X.</li>
 * </ul>
 * 
 * If you want to compile the native library for your system, compile
 * <code>de_jstacs_utils_UserTime.c</code> as a dynamic library including the
 * path to <code>jni.h</code> on your system (
 * <code>-I/path/to/directory_containing_jni.h</code>).
 * 
 * @author Jens Keilwagen
 */
public class UserTime extends Time {

	static {
		System.loadLibrary( "UserTime" );
	}

	float start;

	float ticks;

	/* (non-Javadoc)
	 * @see de.jstacs.utils.Time#getElapsedTime()
	 */
	@Override
	public double getElapsedTime() {
		return ( getUserTime() - start ) / ticks;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.Time#reset()
	 */
	@Override
	public void reset() {
		ticks = getTicks();
		start = getUserTime();
	}

	private native float getUserTime();

	private native long getTicks();
}
