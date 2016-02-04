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
import java.util.HashSet;

import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;
import de.jstacs.tools.ToolResult;

/**
 * {@link ResultSaver} for a {@link ResultSetResult}.
 * File names for the {@link Result}s in the contained {@link ResultSet} are defined using the
 * {@link Result#getName()} method, replacing whitespace with underscores. If a file with the corresponding name
 * already exists, a number is appended and the existing file is not overwritten.
 *  
 * @author Jan Grau
 *
 */
public class ResultSetResultSaver implements ResultSaver<ResultSetResult> {

	/**
	 * Registers this {@link ResultSaver} in the {@link ResultSaverLibrary}
	 */
	public static void register(){
		ResultSaverLibrary.register( ResultSetResult.class, new ResultSetResultSaver() );
		ResultSaverLibrary.register( ToolResult.class, new ResultSetResultSaver() );
	}
	
	private ResultSetResultSaver() {
	}

	@Override
	public boolean isAtomic() {
		return false;
	}

	@Override
	public String[] getFileExtensions( ResultSetResult result ) {
		return null;
	}

	@Override
	public boolean writeOutput( ResultSetResult result, File dir ) {
		if(!dir.exists()){
			dir.mkdirs();
		}
		
		if(!dir.isDirectory()){
			return false;
		}
		ResultSet set = result.getRawResult()[0];
		HashSet<String> names = new HashSet<String>();
		boolean wroteAll = true;
		for(int i=0;i<set.getNumberOfResults();i++){
			Result res = set.getResultAt( i );
			ResultSaver saver = ResultSaverLibrary.getSaver( res.getClass() );
			if(saver != null){
				String filename = res.getName().replaceAll( "[\\s\\:\\/]", "_" );
				String temp = filename;
				if(saver.isAtomic()){
					temp = temp + "." + saver.getFileExtensions( res )[0];
				}
				int j=1;
				while(names.contains( temp ) || new File(dir.getAbsolutePath()+File.separator+temp).exists()){//TODO correct?
					temp = filename+"_"+j;
					if(saver.isAtomic()){
						temp = temp + "." + saver.getFileExtensions( res )[0];
					}
					j++;
				}
				filename = temp;
				names.add( filename );
				
				
				wroteAll &= saver.writeOutput( res, new File(dir.getAbsolutePath()+File.separator+filename) );
				
			}
		}
		return wroteAll;
	}

	@Override
	public boolean writeOutput( ResultSetResult result, StringBuffer buf ) {
		throw new RuntimeException( "Not possible" );
	}

	
	
}
