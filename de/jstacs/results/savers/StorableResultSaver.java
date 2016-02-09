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
import java.io.IOException;

import de.jstacs.io.FileManager;
import de.jstacs.results.StorableResult;

/**
 * Implements a {@link ResultSaver} for {@link StorableResult}. {@link StorableResult} are stored to files as their XML representation
 * using {@link de.jstacs.Storable#toXML()} of the contained {@link de.jstacs.Storable} object. This object is obtained from the {@link StorableResult#getResultInstance()}
 * method.
 * 
 * @author Jan Grau
 */
public class StorableResultSaver implements ResultSaver<StorableResult> {


	/**
	 * Registers this {@link ResultSaver} in the {@link ResultSaverLibrary}
	 */
	public static void register(){
		ResultSaverLibrary.register( StorableResult.class, new StorableResultSaver() );
	}
	
	private StorableResultSaver() { }

	@Override
	public boolean isAtomic() {
		return true;
	}

	@Override
	public String[] getFileExtensions(StorableResult result) {
		return new String[]{"xml"};
	}

	@Override
	public boolean writeOutput(StorableResult result, File path) {
		try{
			FileManager.writeFile(path, result.getResultInstance().toXML());
			return true;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
	}

	@Override
	public boolean writeOutput(StorableResult result, StringBuffer buf) {
		buf.append(result.getResultInstance().toXML());
		return true;
	}
	
	

	
	
}
