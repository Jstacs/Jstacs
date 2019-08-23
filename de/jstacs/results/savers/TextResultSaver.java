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
import java.nio.file.Path;
import java.nio.file.Paths;

import de.jstacs.io.FileManager;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.results.TextResult;

/**
 * {@link ResultSaver} for a {@link TextResult}.
 * The contents of the {@link TextResult} are saved to disk using the {@link FileRepresentation#getContent()} method.
 * If the file defined by the contained {@link FileRepresentation} already exists, it is simply copied to the new location.
 * 
 * 
 * @author Jan Grau, Jens Keilwagen
 *
 */
public class TextResultSaver implements ResultSaver<TextResult> {

	/**
	 * Registers this {@link ResultSaver} in the {@link ResultSaverLibrary}
	 */
	public static void register(){
		ResultSaverLibrary.register( TextResult.class, new TextResultSaver() );
	}

	private TextResultSaver() {
	}



	@Override
	public String[] getFileExtensions( TextResult result ) {
		String[] mimes = result.getMime().split( "\\," );
		return mimes;
	}

	@Override
	public boolean writeOutput( TextResult result, File path ) {
		try{
			FileRepresentation rep = result.getValue();
			String relPath = getRelative( path.getAbsolutePath() );
			if(rep.getFilename() != null && (new File(rep.getFilename())).exists()){
				FileManager.copy(rep.getFilename(), relPath );
				new File( rep.getFilename() ).delete();////TODO move (Jan likes to move files);
			}else{
				PrintWriter wr = new PrintWriter( path );
				wr.println( rep.getContent() );
				wr.close();
			}
			rep.setFilename( relPath );
			return true;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
	}
	
	private String getRelative( String path ) {
		Path home = Paths.get(new File("").getAbsolutePath());
		Path absolute = Paths.get(path);
		try {
			return home.relativize(absolute).toString();
		} catch( IllegalArgumentException e ) {
			return path;
		}
	}

	@Override
	public boolean isAtomic() {
		return true;
	}

	@Override
	public boolean writeOutput( TextResult result, StringBuffer buf ) {
		buf.append( result.getValue().getContent() );
		return true;
	}

}
