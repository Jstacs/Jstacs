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
import java.io.PrintWriter;
import java.io.StringWriter;

import de.jstacs.results.ListResult;

/**
 * {@link ResultSaver} for {@link ListResult}s.
 * Contents are saved to a text file or {@link StringBuffer} using the {@link ListResult#print(PrintWriter)} method.
 * 
 * @author Jan Grau
 *
 */
public class ListResultSaver implements ResultSaver<ListResult> {

	/**
	 * Registers this {@link ResultSaver} in the {@link ResultSaverLibrary}
	 */
	public static void register(){
		ResultSaverLibrary.register( ListResult.class, new ListResultSaver() );
	}
	
	private ListResultSaver(){}
	
	@Override
	public boolean isAtomic() {
		return true;
	}

	@Override
	public String[] getFileExtensions( ListResult result ) {
		return new String[]{"tsv"};
	}

	@Override
	public boolean writeOutput( ListResult result, File path ) {
		try{
			PrintWriter wr = new PrintWriter( path );
			result.print( wr );
			wr.close();
			return true;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
	}

	@Override
	public boolean writeOutput( ListResult result, StringBuffer buf ) {
		StringWriter sw = new StringWriter();
		result.print( new PrintWriter( sw ) );
		buf.append(sw.toString());
		return true;
	}

	
	
	
	
}
