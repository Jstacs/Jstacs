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

package de.jstacs.results.savers;

import java.io.File;

import de.jstacs.results.Result;

/**
 * Interface for saving {@link Result}s to a file.
 * Specific implementations for different {@link Result} types required.
 * Each implementation must register using {@link ResultSaverLibrary#register(Class, ResultSaver)} to become available
 * in other classes.
 * 
 * @author Jan Grau
 *
 * @param <T> the type of the {@link Result}s that can be saved
 */
public interface ResultSaver<T extends Result> {

	/**
	 * Returns <code>true<code> if this {@link ResultSaver} is for storing atomic {@link Result}s.
	 * @return if this {@link ResultSaver} is for storing atomic {@link Result}s.
	 */
	public boolean isAtomic();
	
	/**
	 * Returns the file extensions (in descending preference) for storing the given {@link Result}
	 * @param result the result
	 * @return the file extension(s)
	 */
	public String[] getFileExtensions(T result);
	
	/**
	 * Writes the output (i.e., the result contents) to the supplied file.
	 * @param result the result
	 * @param path the output file
	 * @return if writing was successful
	 */
	public boolean writeOutput(T result, File path);
	
	/**
	 * Appends the output (i.e., the result contents) to the supplied {@link StringBuffer}
	 * @param result the result
	 * @param buf the buffer
	 * @return if appending was successful
	 */
	public boolean writeOutput(T result, StringBuffer buf);
	
}
